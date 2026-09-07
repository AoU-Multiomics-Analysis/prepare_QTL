version 1.0

import "tasks/tensorqtl_trans_compat.wdl" as trans

workflow PhenotypePCQTL {
    input {
        File covariates
        File plink_pgen
        File plink_pvar
        File plink_psam
        String placeholder_chromosome = "chrY"
        Float maf_threshold = 0.05
        Float pval_threshold = 0.00001
        Int num_threads = 4
        Int memory_gb = 64
        Int disk_gb = 100
        Int preemptible_attempts = 1
        String tensorqtl_docker = "gcr.io/broad-cga-francois-gtex/tensorqtl@sha256:f6efb9e592eb32c46cb75070be2769b34381d60cbb2709d2885771324abfe32a"
    }

    call PreparePCPhenotypes {
        input:
            covariates = covariates,
            placeholder_chromosome = placeholder_chromosome,
            docker_image = tensorqtl_docker
    }

    call trans.TensorQTLTrans {
        input:
            plink_pgen = plink_pgen,
            plink_pvar = plink_pvar,
            plink_psam = plink_psam,
            phenotype_bed = PreparePCPhenotypes.phenotype_bed,
            covariates = PreparePCPhenotypes.remaining_covariates,
            maf_threshold = maf_threshold,
            pval_threshold = pval_threshold,
            num_threads = num_threads,
            memory_gb = memory_gb,
            disk_gb = disk_gb,
            preemptible_attempts = preemptible_attempts,
            docker_image = tensorqtl_docker
    }

    output {
        File phenotype_pc_bed = PreparePCPhenotypes.phenotype_bed
        File remaining_covariates = PreparePCPhenotypes.remaining_covariates
        File pc_selection = PreparePCPhenotypes.pc_selection
        File trans_qtl_pairs = TensorQTLTrans.trans_qtl_pairs
        File preparation_log = PreparePCPhenotypes.log
        File trans_log = TensorQTLTrans.log
    }

    parameter_meta {
        covariates: "TensorQTL TSV: covariate IDs in the first column, samples in the header."
        placeholder_chromosome: "Default chrY. Nearby chrY variants can be filtered. Use a label absent from PVAR (e.g. PC_SCAN) to avoid coordinate-based exclusion."
        plink_pvar: "Uncompressed PLINK 2 PVAR with #CHROM header. Files can have different basenames and source directories."
        pval_threshold: "Nominal P-value cutoff for saving sparse pairs; not a PC exclusion rule or FDR threshold."
    }
}

task PreparePCPhenotypes {
    input {
        File covariates
        String placeholder_chromosome
        String docker_image
    }

    command <<<
        set -euo pipefail
        echo '[prepare] Starting phenotype PC preparation'
        cat > prepare_pc_scan.py <<'PYTHON_SOURCE'
# BEGIN EMBEDDED PYTHON
"""Split an existing tensorQTL covariate matrix; do not recompute PCA."""
import argparse
import csv
import math
from pathlib import Path
import re
import sys


def readable(value):
    if '://' in value:
        raise ValueError('localization error: unresolved URI: ' + value)
    p = Path(value)
    if not p.is_file():
        raise ValueError('localization error: input is not a readable file: ' + value)
    with p.open():
        pass
    return p


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--covariates', required=True)
    parser.add_argument('--chromosome', default='chrY')
    args = parser.parse_args()
    matrix = readable(args.covariates)
    if not re.fullmatch(r'[A-Za-z0-9_.-]+', args.chromosome):
        raise ValueError('placeholder chromosome must contain only letters, digits, _, . or -')
    print('[prepare] Reading covariate matrix', flush=True)
    with matrix.open(newline='') as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader, [])
        samples = header[1:]
        if len(samples) < 3 or any(not x.strip() for x in samples) or len(set(samples)) != len(samples):
            raise ValueError('need at least three unique, nonempty sample IDs in the header')
        rows = {}
        for line, row in enumerate(reader, 2):
            if len(row) != len(header) or not row[0].strip():
                raise ValueError(f'malformed matrix row {line}')
            if row[0] in rows:
                raise ValueError('duplicate covariate ID: ' + row[0])
            try:
                values = [float(x) for x in row[1:]]
            except ValueError:
                raise ValueError('nonnumeric values in covariate: ' + row[0])
            if not all(math.isfinite(x) for x in values):
                raise ValueError('nonfinite values in covariate: ' + row[0])
            if len(set(values)) < 2:
                raise ValueError('constant covariate/phenotype PC: ' + row[0])
            rows[row[0]] = row[1:]
    ids = [name for name in rows if re.fullmatch(r'PC[0-9]+', name)]
    if not ids:
        raise ValueError('no molecular PC rows found; expected case-sensitive names PC1, PC2, etc.')
    selected = set(ids)
    remaining = [x for x in rows if x not in selected]
    if not remaining:
        raise ValueError('no remaining covariates; keep genetic PCs and other adjustment covariates')
    if len(samples) <= len(remaining) + 2:
        raise ValueError('insufficient samples for remaining covariates plus intercept and genotype')
    with Path('pc_scan.phenotypes.bed').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['#chr', 'start', 'end', 'phenotype_id'] + samples)
        for name in ids:
            writer.writerow([args.chromosome, 0, 1, name] + rows[name])
    with Path('pc_scan.covariates.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        for name in remaining:
            writer.writerow([name] + rows[name])
    with Path('pc_scan.selection.tsv').open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['covariate_id', 'role'])
        writer.writerows((name, 'phenotype' if name in selected else 'covariate') for name in rows)
    print(f'[prepare] Wrote {len(ids)} PC phenotypes, {len(remaining)} covariates, {len(samples)} samples', flush=True)
    print(f'[prepare] Placeholder {args.chromosome}:0-1; nearby variants on this chromosome can be excluded by tensorQTL', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        sys.exit(1)
# END EMBEDDED PYTHON
PYTHON_SOURCE
        python3 prepare_pc_scan.py \
            --covariates '~{sub(covariates, "'", "'\"'\"'")}' \
            --chromosome '~{sub(placeholder_chromosome, "'", "'\"'\"'")}' 2>&1 | tee preparation.log
    >>>

    output {
        File phenotype_bed = "pc_scan.phenotypes.bed"
        File remaining_covariates = "pc_scan.covariates.tsv"
        File pc_selection = "pc_scan.selection.tsv"
        File log = "preparation.log"
    }

    runtime {
        docker: docker_image
        cpu: 1
        memory: "2 GB"
        disks: "local-disk 10 SSD"
        bootDiskSizeGb: 25
        preemptible: 1
    }
}

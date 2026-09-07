version 1.0

# WDL 1.0 adaptation of tensorQTL_trans task from:
# https://github.com/AoU-Multiomics-Analysis/tensorQTL_trans/blob/ab5aa7befaae3c6d0b5f20178cb587f10b4da160/tensorQTL_trans.wdl
# Original task author: Francois Aguet. See README for compatibility changes.
task TensorQTLTrans {
    input {
        File plink_pgen
        File plink_pvar
        File plink_psam
        File phenotype_bed
        File covariates
        Float maf_threshold
        Float pval_threshold
        Int num_threads
        Int memory_gb
        Int disk_gb
        Int preemptible_attempts
        String docker_image
    }

    command <<<
        set -euo pipefail
        echo '[trans] Starting localized tensorQTL task'
        export OMP_NUM_THREADS=~{num_threads}
        export MKL_NUM_THREADS=~{num_threads}
        cat > run_trans_pc.py <<'PYTHON_SOURCE'
# BEGIN EMBEDDED PYTHON
"""Stage localized PLINK files and run the upstream tensorQTL trans command."""
import argparse
import csv
from pathlib import Path
import subprocess
import sys


def readable(value):
    if '://' in value:
        raise ValueError('localization error: unresolved URI: ' + value)
    p = Path(value)
    if not p.is_file():
        raise ValueError('localization error: input is not a readable file: ' + value)
    with p.open('rb'):
        pass
    return p.resolve()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['plink-pgen', 'plink-pvar', 'plink-psam', 'phenotype-bed', 'covariates']:
        parser.add_argument('--' + name, required=True)
    parser.add_argument('--maf-threshold', type=float, default=0.05)
    parser.add_argument('--pval-threshold', type=float, default=1e-5)
    args = parser.parse_args()
    files = {name: readable(getattr(args, name)) for name in ['plink_pgen', 'plink_pvar', 'plink_psam', 'phenotype_bed', 'covariates']}
    if not 0 <= args.maf_threshold <= 0.5 or not 0 < args.pval_threshold <= 1:
        raise ValueError('MAF threshold must be in [0,0.5] and P-value threshold in (0,1]')
    print('[trans] Checking localized inputs and sample IDs', flush=True)
    with files['phenotype_bed'].open() as handle:
        samples = next(csv.reader(handle, delimiter='\t'))[4:]
    with files['covariates'].open() as handle:
        cov_samples = next(csv.reader(handle, delimiter='\t'))[1:]
    if not samples or samples != cov_samples or len(set(samples)) != len(samples):
        raise ValueError('phenotype and covariate sample IDs/order must match and be unique')
    with files['plink_psam'].open() as handle:
        header = None
        genotype_samples = []
        for line in handle:
            if line.startswith('##') or not line.strip():
                continue
            fields = line.split()
            if header is None:
                header = [x.lstrip('#') for x in fields]
                if 'IID' not in header:
                    raise ValueError('PLINK PSAM requires an IID header')
                col = header.index('IID')
            else:
                if len(fields) != len(header):
                    raise ValueError('malformed PLINK PSAM row')
                genotype_samples.append(fields[col])
    if len(set(genotype_samples)) != len(genotype_samples):
        raise ValueError('duplicate PLINK IID values')
    missing = set(samples) - set(genotype_samples)
    if missing:
        raise ValueError(f'{len(missing)} phenotype samples are missing from PLINK PSAM')
    with files['plink_pvar'].open() as handle:
        chromosome_col = None
        has_placeholder = False
        with files['phenotype_bed'].open() as bed:
            reader = csv.reader(bed, delimiter='\t'); next(reader)
            placeholders = {(row[0][3:] if row[0].startswith('chr') else row[0]) for row in reader}
        for line in handle:
            if line.startswith('##') or not line.strip():
                continue
            fields = line.split()
            if chromosome_col is None:
                header = [x.lstrip('#') for x in fields]
                if 'CHROM' not in header:
                    raise ValueError('PVAR must be uncompressed text with a #CHROM header')
                chromosome_col = header.index('CHROM')
            elif (fields[chromosome_col][3:] if fields[chromosome_col].startswith('chr') else fields[chromosome_col]) in placeholders:
                has_placeholder = True
                break
    if chromosome_col is None:
        raise ValueError('empty PVAR file')
    if has_placeholder:
        print('[trans] WARNING: variants exist on the placeholder chromosome; tensorQTL can remove pairs within 5 Mb', flush=True)
    stage = Path('staged_genotypes')
    stage.mkdir()
    for ext in ['pgen', 'pvar', 'psam']:
        (stage / ('input.' + ext)).symlink_to(files['plink_' + ext])
    command = [sys.executable, '-m', 'tensorqtl', str(stage / 'input'), str(files['phenotype_bed']), 'pc_scan',
               '--mode', 'trans', '--covariates', str(files['covariates']),
               '--maf_threshold', str(args.maf_threshold), '--pval_threshold', str(args.pval_threshold)]
    print('[trans] Starting tensorQTL trans (sparse output)', flush=True)
    subprocess.run(command, check=True)
    if not Path('pc_scan.trans_qtl_pairs.parquet').is_file():
        raise ValueError('tensorQTL did not produce the expected sparse association file')
    print('[trans] Finished tensorQTL trans', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        sys.exit(1)
# END EMBEDDED PYTHON
PYTHON_SOURCE
        python3 run_trans_pc.py \
            --plink-pgen '~{sub(plink_pgen, "'", "'\"'\"'")}' \
            --plink-pvar '~{sub(plink_pvar, "'", "'\"'\"'")}' \
            --plink-psam '~{sub(plink_psam, "'", "'\"'\"'")}' \
            --phenotype-bed '~{sub(phenotype_bed, "'", "'\"'\"'")}' \
            --covariates '~{sub(covariates, "'", "'\"'\"'")}' \
            --maf-threshold ~{maf_threshold} \
            --pval-threshold ~{pval_threshold} 2>&1 | tee trans.log
    >>>

    output {
        File trans_qtl_pairs = "pc_scan.trans_qtl_pairs.parquet"
        File log = "trans.log"
    }

    runtime {
        docker: docker_image
        cpu: num_threads
        memory: "~{memory_gb} GB"
        disks: "local-disk ~{disk_gb} SSD"
        bootDiskSizeGb: 25
        preemptible: preemptible_attempts
    }
}

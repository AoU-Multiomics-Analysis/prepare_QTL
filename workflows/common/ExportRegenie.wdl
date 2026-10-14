version 1.0

# Python and REGENIE are already available in this tested image.
task ExportRegenieFiles {
    input {
        File Bed
        File? Covariates
        String CovariateFormat = "sample"
        String OutputPrefix = "qtl"
        Int PhenotypesPerFile = 100
        Int DiskGB = 100
        String DockerImage = "ghcr.io/aou-multiomics-analysis/regenie@sha256:4232a95edf04477deab443032cf8cff01a8fc3cd200e6a9505721230a2795c43"
    }
    command <<<
        set -euo pipefail
        {
        echo "stage=export_regenie start_time=$(date -u +%Y-%m-%dT%H:%M:%SZ) dimensions=pending"
        cat > export_regenie.py <<'REGENIE_EXPORT_SOURCE'
# BEGIN REGENIE EXPORT
"""Export prepared QTL BEDs as bounded REGENIE phenotype batches."""
import argparse
import csv
import gzip
import hashlib
import math
from pathlib import Path
import re
import sqlite3
import sys
import tempfile


MISSING = {'', '.', 'na', 'nan', '-nan'}
SAFE_NAME = re.compile(r'[A-Za-z0-9_][A-Za-z0-9_.-]*\Z')


def readable(value):
    if '://' in value:
        raise ValueError('localization error: unresolved URI: ' + value)
    path = Path(value)
    if not path.is_file():
        raise ValueError('localization error: input is not a readable file: ' + value)
    return path


def rows(path):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', encoding='utf-8', newline='') as stream:
        # Input values are literal TSV fields, not quoted CSV fields.
        yield from csv.reader(stream, delimiter='\t', quoting=csv.QUOTE_NONE)


def unique_ids(values, label):
    if not values or any(not v or re.search(r'\s', v) or v.lower() in MISSING for v in values):
        raise ValueError(label + ' contain blank, missing, or whitespace IDs')
    if len(values) != len(set(values)):
        raise ValueError(label + ' contain duplicate IDs')


def covariates(path, orientation):
    table = rows(path)
    header = next(table, [])
    if not header or header[0] not in {'ID', 'sample_id', 'IID'}:
        raise ValueError('covariates must start with ID, sample_id, or IID')
    data = list(table)
    if any(len(row) != len(header) for row in data):
        raise ValueError('covariates have an inconsistent number of fields')
    if orientation == 'tensor':
        samples = header[1:]
        names = [r[0] for r in data]
        values = list(zip(*(r[1:] for r in data))) if data else [() for _ in samples]
    else:
        samples = [r[0] for r in data]
        names = header[1:]
        values = [r[1:] for r in data]
    unique_ids(samples, 'covariate samples')
    if len(names) != len(set(names)) or any(not SAFE_NAME.fullmatch(n) or n in {'FID', 'IID'} for n in names):
        raise ValueError('covariate names must be unique REGENIE names, excluding FID and IID')
    return names, dict(zip(samples, values))


def covariate_value(value):
    # R's t(data.frame) can add numeric padding in MergeCovariates.R.
    value = value.strip()
    if value.lower() in MISSING or re.search(r'\s', value):
        raise ValueError('covariates contain missing values or whitespace')
    try:
        numeric = float(value)
    except ValueError:
        return value  # Categorical predictor; the REGENIE run must select it.
    if not math.isfinite(numeric):
        raise ValueError('covariates contain nonfinite values')
    return value


def phenotype_value(value):
    if value.lower() in MISSING:
        return 'NA'
    try:
        numeric = float(value)
    except ValueError as error:
        raise ValueError('phenotypes contain a nonnumeric value') from error
    if not math.isfinite(numeric):
        raise ValueError('phenotypes contain an infinite value')
    # Preserve all numeric text, including decimal precision and exponent notation.
    return value


def writer(path, header):
    stream = path.open('w', encoding='utf-8', newline='')
    output = csv.writer(stream, delimiter='\t', lineterminator='\n', quoting=csv.QUOTE_NONE)
    output.writerow(header)
    return stream, output


def export(args):
    if args.phenotypes_per_file < 1:
        raise ValueError('phenotypes-per-file must be positive')
    if not SAFE_NAME.fullmatch(args.output_prefix) or args.output_prefix in {'.', '..'}:
        raise ValueError('output-prefix must be a plain filename with letters, digits, _, . or -')
    bed = rows(readable(args.bed))
    header = next(bed, [])
    if len(header) < 5 or header[0] not in {'#chr', '#chrom', '#CHROM', 'chr', 'chrom'}:
        raise ValueError('BED must have four metadata columns and sample columns')
    samples = header[4:]
    unique_ids(samples, 'BED samples')
    names, cov = covariates(readable(args.covariates), args.covariate_format) if args.covariates else ([], None)
    indices = [i for i, sample in enumerate(samples) if cov is None or sample in cov]
    if not indices:
        raise ValueError('no samples overlap the BED and covariates')
    selected = [samples[i] for i in indices]
    output = Path('outputs')
    output.mkdir(exist_ok=True)
    prefix = args.output_prefix
    if list(output.glob(prefix + '.*')):
        raise ValueError('output files already exist for this prefix')
    excluded_file, excluded = writer(output / (prefix + '.excluded_samples.tsv'), ['IID', 'reason'])
    with excluded_file:
        for sample in samples:
            if cov is not None and sample not in cov:
                excluded.writerow([sample, 'missing_covariate_sample'])
    if names:
        cov_file, cov_out = writer(output / (prefix + '.covariates.tsv'), ['FID', 'IID'] + names)
        with cov_file:
            for sample in selected:
                cov_out.writerow(['0', sample] + [covariate_value(v) for v in cov[sample]])
    trait_file, traits = writer(output / (prefix + '.traits.tsv'),
                               ['phenotype_id', 'source_id', 'chromosome', 'start', 'end', 'phenotype_file'])
    batch_file, batches = writer(output / (prefix + '.batches.tsv'),
                                 ['phenotype_file', 'phenotype_columns', 'n_traits', 'n_samples'])
    total = batch_count = 0
    block_names, block_values = [], []

    def flush():
        nonlocal batch_count
        if not block_names:
            return
        batch_count += 1
        filename = f'{prefix}.phenotypes.{batch_count:06d}.tsv'
        stream, out = writer(output / filename, ['FID', 'IID'] + block_names)
        with stream:
            for sample, values in zip(selected, zip(*block_values)):
                out.writerow(['0', sample, *values])
        batches.writerow([filename, ','.join(block_names), len(block_names), len(selected)])
        print(f'Exported batch={batch_count} traits={len(block_names)} samples={len(selected)}', flush=True)
        block_names.clear()
        block_values.clear()

    with trait_file, batch_file, tempfile.TemporaryDirectory(prefix='regenie_ids_') as directory:
        with sqlite3.connect(str(Path(directory) / 'ids.sqlite')) as db:
            db.execute('CREATE TABLE ids (source TEXT PRIMARY KEY, alias TEXT UNIQUE)')
            for line_number, row in enumerate(bed, 2):
                if len(row) != len(header):
                    raise ValueError(f'BED row {line_number} has an inconsistent number of fields')
                chromosome, start, end, source_id = row[:4]
                if not source_id or any(c in source_id for c in '\r\n\t'):
                    raise ValueError('BED has a blank or invalid phenotype ID')
                alias = 'pheno_' + hashlib.sha256(source_id.encode('utf-8')).hexdigest()[:24]
                try:
                    db.execute('INSERT INTO ids VALUES (?, ?)', (source_id, alias))
                except sqlite3.IntegrityError as error:
                    raise ValueError('duplicate phenotype ID or alias collision: ' + source_id) from error
                # Validate all phenotype cells, including covariate-excluded samples.
                values = [phenotype_value(v) for v in row[4:]]
                block_names.append(alias)
                block_values.append([values[i] for i in indices])
                filename = f'{prefix}.phenotypes.{batch_count + 1:06d}.tsv'
                traits.writerow([alias, source_id, chromosome, start, end, filename])
                total += 1
                if len(block_names) == args.phenotypes_per_file:
                    flush()
            flush()
    if not total:
        raise ValueError('BED contains no phenotypes')
    print(f'Complete traits={total} samples={len(selected)} excluded_samples={len(samples)-len(selected)} batches={batch_count}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bed', required=True)
    parser.add_argument('--covariates')
    parser.add_argument('--covariate-format', choices=['sample', 'tensor'], default='sample')
    parser.add_argument('--output-prefix', default='qtl')
    parser.add_argument('--phenotypes-per-file', type=int, default=100)
    args = parser.parse_args()
    try:
        export(args)
    except (ValueError, OSError, csv.Error, sqlite3.Error) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
# END REGENIE EXPORT
REGENIE_EXPORT_SOURCE
        args=(--bed '~{sub(Bed, "'", "'\"'\"'")}'
              --output-prefix '~{sub(OutputPrefix, "'", "'\"'\"'")}'
              --phenotypes-per-file ~{PhenotypesPerFile})
        ~{if defined(Covariates) then "args+=(--covariates '" + sub(select_first([Covariates]), "'", "'\"'\"'") + "' --covariate-format '" + sub(CovariateFormat, "'", "'\"'\"'") + "')" else ""}
        python3 export_regenie.py "${args[@]}"
        echo "stage=export_regenie completion_time=$(date -u +%Y-%m-%dT%H:%M:%SZ) outputs=outputs/"
        } 2>&1 | tee export.log
    >>>
    runtime {
        docker: DockerImage
        memory: "4G"
        cpu: 1
        disks: "local-disk ~{DiskGB} HDD"
        preemptible: 0
        maxRetries: 1
    }
    output {
        Array[File] Phenotypes = glob("outputs/*.phenotypes.*.tsv")
        File? CovariatesOut = "outputs/~{OutputPrefix}.covariates.tsv"
        File TraitMap = "outputs/~{OutputPrefix}.traits.tsv"
        File Batches = "outputs/~{OutputPrefix}.batches.tsv"
        File ExcludedSamples = "outputs/~{OutputPrefix}.excluded_samples.tsv"
        File Log = "export.log"
    }
}

workflow ExportRegenie {
    input {
        File Bed
        File? Covariates
        String CovariateFormat = "sample"
        String OutputPrefix = "qtl"
        Int PhenotypesPerFile = 100
        Int DiskGB = 100
    }
    call ExportRegenieFiles {
        input:
            Bed = Bed,
            Covariates = Covariates,
            CovariateFormat = CovariateFormat,
            OutputPrefix = OutputPrefix,
            PhenotypesPerFile = PhenotypesPerFile,
            DiskGB = DiskGB
    }
    output {
        Array[File] Phenotypes = ExportRegenieFiles.Phenotypes
        File? CovariatesOut = ExportRegenieFiles.CovariatesOut
        File TraitMap = ExportRegenieFiles.TraitMap
        File Batches = ExportRegenieFiles.Batches
        File ExcludedSamples = ExportRegenieFiles.ExcludedSamples
        File Log = ExportRegenieFiles.Log
    }
}

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

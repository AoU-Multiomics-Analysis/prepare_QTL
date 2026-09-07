"""Normalize already selected trans associations without recomputing significance."""
import argparse
import csv
import gzip
import json
import re
from pathlib import Path
from .common import ASSOCIATION_FIELDS, assignments, chromosome, finite, lines, local_file_list, log, rows, sizes, token

def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['association-list', 'ancestry', 'chrom-sizes', 'chromosomes', 'genome-build']:
        p.add_argument('--' + name, required=True)
    p.add_argument('--pvalue-threshold', type=float, required=True)
    for name, default in [('variant-column','variant_id'),('phenotype-column','phenotype_id'),('pvalue-column','pval'),('chrom-column','chrom'),('position-column','pos'),('modality-column','modality'),('dataset-column','dataset_id'),('cell-type-column','cell_type')]:
        p.add_argument('--'+name, default=default)
    p.add_argument('--default-modality', default='unspecified')
    args = p.parse_args()
    files = local_file_list(args.association_list)
    if not files:
        raise ValueError('supply at least one association table')
    group_map = assignments(args.ancestry)
    lengths = sizes(args.chrom_sizes)
    chroms = [chromosome(x) for x in lines(args.chromosomes)]
    if not chroms or len(chroms) != len(set(chroms)) or set(chroms) - set(lengths):
        raise ValueError('genotype chromosomes must be unique, nonempty and have chromosome sizes')
    if not 0 < args.pvalue_threshold <= 1:
        raise ValueError('pvalue-threshold must be in (0,1]')
    seen = {}
    counts = {'input_rows':0, 'retained_rows':0}
    log('Normalizing association tables')
    with gzip.open('associations.tsv.gz', 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=ASSOCIATION_FIELDS, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        for file_index, file in enumerate(files, 1):
            for row in rows(file):
                counts['input_rows'] += 1
                pv = finite(row[args.pvalue_column], 'p-value')
                if not 0 <= pv <= 1:
                    raise ValueError('p-value outside [0,1]')
                if pv > args.pvalue_threshold:
                    continue
                variant = token(row[args.variant_column], 'variant ID')
                phenotype = token(row[args.phenotype_column], 'phenotype ID')
                has_chrom = args.chrom_column in row
                has_pos = args.position_column in row
                if has_chrom != has_pos:
                    raise ValueError('supply both chromosome and position columns, or neither')
                if has_chrom:
                    chrom = chromosome(row[args.chrom_column])
                    pos = int(str(row[args.position_column]))
                else:
                    match = re.fullmatch(r'([^:_]+)[:_](\d+)[:_]([^:_]+)[:_]([^:_]+)', variant)
                    if not match:
                        raise ValueError('variant ID needs explicit chromosome/position columns: ' + variant)
                    chrom, pos = chromosome(match[1]), int(match[2])
                if chrom not in chroms or not 1 <= pos <= lengths[chrom]:
                    raise ValueError('association chromosome lacks genotype shard, or position is outside chromosome: ' + variant)
                if variant in seen and seen[variant] != (chrom,pos):
                    raise ValueError('variant ID maps to multiple positions: ' + variant)
                seen[variant] = (chrom,pos)
                out = dict(variant_id=variant, chrom=chrom, pos=pos, phenotype_id=phenotype, p_value=pv,
                           modality=token(row.get(args.modality_column,args.default_modality),'modality'),
                           dataset_id=token(row.get(args.dataset_column,'input_'+str(file_index)),'dataset ID'),
                           cell_type=token(row.get(args.cell_type_column,'unspecified'),'cell type'))
                writer.writerow(out)
                counts['retained_rows'] += 1
    Path('ancestries.txt').write_text(''.join(x+'\n' for x in sorted(set(group_map.values()))))
    Path('preparation.json').write_text(json.dumps(dict(counts, unique_seeds=len(seen), ancestry_samples={g:sum(v==g for v in group_map.values()) for g in sorted(set(group_map.values()))}, parameters=vars(args)), indent=2)+'\n')
    log(f"Retained {counts['retained_rows']} association rows and {len(seen)} unique starting variants")

if __name__ == '__main__':
    main()

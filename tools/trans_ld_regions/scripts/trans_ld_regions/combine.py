"""Combine ancestry evidence and preserve phenotype links in conservative regions."""
import argparse
from collections import defaultdict
import csv
import gzip
import json
from pathlib import Path
from .common import ASSOCIATION_FIELDS, chrom_key, interval, local_file_list, log, rows, seeds, sizes, write_table

def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['associations','interval-list','chrom-sizes','genome-build']:
        p.add_argument('--'+name, required=True)
    p.add_argument('--padding-bp', type=int, default=100000)
    p.add_argument('--fallback-bp', type=int, default=1000000)
    p.add_argument('--min-region-bp', type=int, default=2000000)
    args = p.parse_args()
    if min(args.padding_bp,args.fallback_bp) < 0 or args.min_region_bp < 1:
        raise ValueError('invalid interval width parameters')
    lengths = sizes(args.chrom_sizes)
    starts = seeds(args.associations)
    evidence = defaultdict(list)
    all_records = []
    for file in local_file_list(args.interval_list):
        for row in rows(file):
            key = row['variant_id']
            if key not in starts or (row['chrom'],int(row['pos'])) != starts[key]:
                raise ValueError('LD evidence does not match association: ' + key)
            if any(r['ancestry'] == row['ancestry'] for r in evidence[key]):
                raise ValueError('duplicate seed/ancestry evidence: '+key)
            evidence[key].append(row)
            all_records.append(row)
    intervals = []
    seed_records = []
    for variant,(chrom,pos) in starts.items():
        if not evidence[variant]:
            raise ValueError('missing LD task output for seed: '+variant)
        assessed = [r for r in evidence[variant] if r['status'] in ('ok','search_limit')]
        flags = set()
        if assessed:
            lo = min(int(r['start']) for r in assessed)
            hi = max(int(r['end']) for r in assessed)
            left = sorted(r['ancestry'] for r in assessed if int(r['start']) == lo)
            right = sorted(r['ancestry'] for r in assessed if int(r['end']) == hi)
            if any(r['status']=='search_limit' for r in assessed): flags.add('search_limit')
            if len(assessed) != len(evidence[variant]): flags.add('some_groups_unassessed')
        else:
            lo, hi = pos-args.fallback_bp, pos+args.fallback_bp
            left, right = [], []
            flags.add('fallback')
        lo, hi = interval(lo-1-args.padding_bp, hi+args.padding_bp, lengths[chrom], args.min_region_bp)
        if hi-lo < args.min_region_bp: flags.add('chromosome_shorter_than_minimum')
        record = dict(variant_id=variant, chrom=chrom, pos=pos, start=lo, end=hi,
                      left_ancestries=','.join(left), right_ancestries=','.join(right),
                      flags=','.join(sorted(flags)), n_assessed_groups=len(assessed))
        seed_records.append(record)
        intervals.append((chrom,lo,hi,variant,flags))
    merged = []
    for chrom,lo,hi,variant,flags in sorted(intervals,key=lambda x:(chrom_key(x[0]),x[1],x[2],x[3])):
        if merged and merged[-1]['chrom']==chrom and lo < merged[-1]['end']:
            merged[-1]['end'] = max(hi,merged[-1]['end'])
            merged[-1]['seeds'].append(variant)
            merged[-1]['flag_set'].update(flags)
        else:
            merged.append(dict(chrom=chrom,start=lo,end=hi,seeds=[variant],flag_set=set(flags)))
    membership = {}
    for region in merged:
        region['region_id'] = f"chr{region['chrom']}_{region['start']}_{region['end']}"
        region['genome_build'] = args.genome_build
        region['width_bp'] = region['end']-region['start']
        region['n_seeds'] = len(region['seeds'])
        region['flags'] = ','.join(sorted(region['flag_set']))
        for variant in region['seeds']: membership[variant] = region
    fields = ['region_id','chrom','start','end','width_bp','n_seeds','genome_build','flags']
    write_table('regions.tsv',fields,merged)
    with open('regions.bed','w') as h:
        for r in merged: h.write(f"chr{r['chrom']}\t{r['start']}\t{r['end']}\t{r['region_id']}\n")
    seed_fields = ['variant_id','chrom','pos','start','end','left_ancestries','right_ancestries','flags','n_assessed_groups','region_id']
    for r in seed_records: r['region_id'] = membership[r['variant_id']]['region_id']
    write_table('seed_regions.tsv',seed_fields,seed_records)
    # Keep every ancestry diagnostic field, including boundary pair support.
    evidence_fields = list(dict.fromkeys(k for r in all_records for k in r)) or ['variant_id','ancestry','status']
    write_table('ancestry_intervals.tsv.gz',evidence_fields,all_records)
    connections = {}
    association_count = 0
    with gzip.open('region_associations.tsv.gz','wt',newline='') as h:
        writer = csv.DictWriter(h,fieldnames=['region_id','start','end']+ASSOCIATION_FIELDS,delimiter='\t',lineterminator='\n')
        writer.writeheader()
        for row in rows(args.associations):
            r = membership[row['variant_id']]
            writer.writerow(dict(row,region_id=r['region_id'],start=r['start'],end=r['end']))
            association_count += 1
            key = (r['region_id'],row['dataset_id'],row['cell_type'],row['modality'],row['phenotype_id'])
            if key not in connections or float(row['p_value']) < float(connections[key]['p_value']):
                connections[key] = dict(window_id=r['region_id'],chrom='chr'+r['chrom'],start=r['start'],end=r['end'],modality=row['modality'],molecular_trait_id=row['phenotype_id'],p_value=row['p_value'],dataset_id=row['dataset_id'],cell_type=row['cell_type'])
    write_table('trans_window_associations.tsv.gz',['window_id','chrom','start','end','modality','molecular_trait_id','p_value','dataset_id','cell_type'],connections.values())
    summary = dict(parameters=vars(args), n_regions=len(merged), n_seeds=len(starts), n_associations=association_count,
                   n_fallback_seeds=sum('fallback' in r['flags'] for r in seed_records),
                   n_search_limit_seeds=sum('search_limit' in r['flags'] for r in seed_records))
    Path('region_qc.json').write_text(json.dumps(summary,indent=2)+'\n')
    log(f'Wrote {len(merged)} regions; all {association_count} association rows retained')

if __name__ == '__main__':
    main()

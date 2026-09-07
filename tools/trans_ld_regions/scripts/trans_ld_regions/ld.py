"""Ancestry-specific pairwise LD with adaptive search and explicit unresolved seeds."""
import argparse
import csv
from decimal import Decimal
import gzip
import json
import math
from pathlib import Path
import subprocess
from .common import assignments, chromosome, log, readable, seeds, sizes, write_table

FIELDS = ['variant_id','chrom','pos','ancestry','start','end','status','n_samples','n_assigned_samples','n_nonfounders_excluded','seed_maf','seed_mac','seed_allele_observations','search_bp','n_ld_partners','left_variant','left_r2','left_mac','right_variant','right_r2','right_mac']

def variant_rows(path, fmt):
    if str(path).endswith('.zst'):
        import zstandard
        stream = zstandard.open(path,'rt')
    elif str(path).endswith('.gz'):
        stream = gzip.open(path,'rt')
    else:
        stream = open(path)
    with stream as h:
        if fmt == 'bed':
            for line in h:
                f=line.split()
                if len(f)!=6: raise ValueError('BIM must have six columns')
                yield f[1],f[0],int(f[3]),',' not in f[4] and ',' not in f[5]
        else:
            header=None
            for line in h:
                if line.startswith('##'): continue
                if header is None:
                    header=line.rstrip().lstrip('#').split()
                    if not {'CHROM','POS','ID','REF','ALT'} <= set(header): raise ValueError('invalid PVAR header')
                    continue
                row=dict(zip(header,line.split()))
                yield row['ID'],row['CHROM'],int(row['POS']),',' not in row['ALT']

def sample_rows(path, fmt):
    with open(path) as h:
        if fmt=='bed':
            for line in h:
                f=line.split()
                if len(f)!=6: raise ValueError('FAM must have six columns')
                yield f[1], f[2]=='0' and f[3]=='0'
        else:
            header=None
            for line in h:
                if line.startswith('##'): continue
                if header is None:
                    header=line.rstrip().lstrip('#').split()
                    if 'IID' not in header: raise ValueError('PSAM requires IID')
                    continue
                row=dict(zip(header,line.split()))
                yield row['IID'],row.get('PAT','0')=='0' and row.get('MAT','0')=='0'

def plink(command):
    log('Running PLINK: '+ ' '.join(command[1:]))
    subprocess.run(command,check=True)

def frequencies(path, minimum):
    counts={}
    with open(path) as h:
        reader=csv.DictReader(h,delimiter='\t')
        for row in reader:
            alt=Decimal(row['ALT_CTS'])
            obs=Decimal(row['OBS_CT'])
            if not alt.is_finite() or not obs.is_finite() or obs<0 or not 0<=alt<=obs:
                raise ValueError('invalid allele counts: '+row['ID'])
            mac=min(alt,obs-alt)
            maf=mac/obs if obs else Decimal(0)
            counts[row['ID']]=dict(maf=float(maf),mac=float(mac),obs=float(obs),eligible=obs>0 and maf>Decimal(str(minimum)))
    return counts

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['genotype','variants','samples','chromosome','ancestry','group','associations','chrom-sizes']:
        p.add_argument('--'+name,required=True)
    p.add_argument('--format',choices=['bed','pgen'],required=True)
    p.add_argument('--maf-threshold',type=float,default=.01)
    p.add_argument('--r2-threshold',type=float,default=.1)
    p.add_argument('--initial-search-bp',type=int,default=1000000)
    p.add_argument('--max-search-bp',type=int,default=10000000)
    p.add_argument('--edge-bp',type=int,default=100000)
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--memory-mb',type=int,default=14000)
    p.add_argument('--plink',default='plink2')
    args=p.parse_args()
    # Validate every required file before reading metadata or starting PLINK.
    files={name:readable(getattr(args,name)) for name in ['genotype','variants','samples','ancestry','associations','chrom_sizes']}
    if not 0<=args.maf_threshold<.5 or not 0<args.r2_threshold<=1:
        raise ValueError('invalid MAF or r2 threshold')
    if not 0<args.edge_bp<args.initial_search_bp<=args.max_search_bp or args.threads<1 or args.memory_mb<640:
        raise ValueError('invalid search or resource settings')
    chrom=chromosome(args.chromosome)
    length=sizes(files['chrom_sizes'])[chrom]
    seed_map=seeds(files['associations'],chrom)
    group_map=assignments(files['ancestry'])
    if args.group not in group_map.values(): raise ValueError('unknown ancestry group')
    selected=[]; seen_samples=set(); nonfounders=0
    for iid,founder in sample_rows(files['samples'],args.format):
        if iid in seen_samples: raise ValueError('IID must be unique in genotype sample file: '+iid)
        seen_samples.add(iid)
        if group_map.get(iid)==args.group:
            if founder: selected.append(iid)
            else: nonfounders+=1
    positions={}; biallelic=set(); seen_variants=set()
    for variant,vc,pos,bi in variant_rows(files['variants'],args.format):
        vc = vc.removeprefix('chr')
        if vc == '23': vc = 'X'
        if vc != chrom: continue
        if variant=='.' or variant in seen_variants: raise ValueError('missing or duplicate genotype variant ID: '+variant)
        if not 1<=pos<=length: raise ValueError('genotype position outside chromosome: '+variant)
        seen_variants.add(variant); positions[variant]=pos
        if bi: biallelic.add(variant)
    if not positions:
        raise ValueError('genotype shard has no variants on declared chromosome: '+chrom)
    for variant,(_,pos) in seed_map.items():
        if variant in positions and pos!=positions[variant]:
            raise ValueError('association/genotype position mismatch: '+variant)
    records={v:dict(variant_id=v,chrom=chrom,pos=pos,ancestry=args.group,start=pos,end=pos,status='missing_variant',n_samples=len(selected),n_assigned_samples=sum(g==args.group for g in group_map.values()),n_nonfounders_excluded=nonfounders,search_bp=0,n_ld_partners=0) for v,(_,pos) in seed_map.items()}
    Path('keep.txt').write_text('#IID\n'+''.join(i+'\n' for i in selected))
    Path('biallelic.txt').write_text(''.join(v+'\n' for v in sorted(biallelic)))
    Path('eligible_variants.txt').write_text('')
    count_map={}
    if args.format=='bed':
        base=[args.plink,'--bed',str(files['genotype']),'--bim',str(files['variants']),'--fam',str(files['samples'])]
    else:
        base=[args.plink,'--pgen',str(files['genotype']),'--pvar',str(files['variants']),'--psam',str(files['samples'])]
    base+=['--keep','keep.txt','--chr',chrom,'--threads',str(args.threads),'--memory',str(args.memory_mb)]
    version=subprocess.run([args.plink,'--version'],check=True,text=True,capture_output=True).stdout.strip()
    log(f'{args.group} chromosome {chrom}: {len(selected)} founders, {len(seed_map)} starting variants')
    if seed_map and len(selected)>=50 and biallelic:
        plink(base+['--extract','biallelic.txt','--freq','counts','--out','frequency'])
        count_map=frequencies('frequency.acount',args.maf_threshold)
        eligible={v for v,c in count_map.items() if c['eligible']}
        Path('eligible_variants.txt').write_text(''.join(v+'\n' for v in sorted(eligible)))
        for v,r in records.items():
            if v not in positions: continue
            if v not in biallelic:
                r['status']='multiallelic'; continue
            if v not in count_map: raise ValueError('PLINK omitted allele counts for '+v)
            c=count_map[v]
            r.update(seed_maf=c['maf'],seed_mac=c['mac'],seed_allele_observations=c['obs'])
            r['status']='ok' if v in eligible else 'maf_filtered'
        active={v for v,r in records.items() if r['status']=='ok'}
        # One record per seed; partners are never promoted into new starting seeds.
        partners={v:set() for v in active}
        radius=args.initial_search_bp
        while active:
            Path('active_seeds.txt').write_text(''.join(v+'\n' for v in sorted(active)))
            for v in active: records[v]['search_bp']=radius
            if len(eligible)>1:
                prefix='ld_'+str(radius)
                plink(base+['--extract','eligible_variants.txt','--r2-unphased','--ld-snp-list','active_seeds.txt','--ld-window',str(len(eligible)+1),'--ld-window-kb',str(radius/1000),'--ld-window-r2',str(args.r2_threshold),'--out',prefix])
                with open(prefix+'.vcor') as h:
                    reader=csv.DictReader(h,delimiter='\t')
                    for pair in reader:
                        r2=float(pair['UNPHASED_R2'])
                        if not math.isfinite(r2) or r2<args.r2_threshold: continue
                        # PLINK orders pairs by genomic position, not by the seed list.
                        for v,other in [(pair['ID_A'],pair['ID_B']),(pair['ID_B'],pair['ID_A'])]:
                            if v not in active: continue
                            if other not in eligible: raise ValueError('LD partner failed eligibility checks')
                            r=records[v]; pos=positions[other]
                            if abs(pos-r['pos'])>radius: continue
                            partners[v].add(other)
                            for side,compare in [('left',pos<r['start']),('right',pos>r['end'])]:
                                if compare:
                                    r['start' if side=='left' else 'end']=pos
                                    r[side+'_variant']=other; r[side+'_r2']=r2; r[side+'_mac']=count_map[other]['mac']
            expand=set()
            for v in active:
                r=records[v]; r['n_ld_partners']=len(partners[v])
                near_left=r['pos']-radius>1 and r['pos']-r['start']>=radius-args.edge_bp
                near_right=r['pos']+radius<length and r['end']-r['pos']>=radius-args.edge_bp
                if near_left or near_right:
                    if radius==args.max_search_bp:
                        r['status']='search_limit'
                        r['start']=max(1,r['pos']-radius); r['end']=min(length,r['pos']+radius)
                    else: expand.add(v)
            active=expand
            radius=min(radius*2,args.max_search_bp)
    else:
        for v,r in records.items():
            if v in positions:
                r['status']='insufficient_samples' if len(selected)<50 else 'multiallelic'
    for record in records.values():
        if record['status'] == 'ok' and record['n_ld_partners'] == 0:
            record['status'] = 'no_reported_ld'
    write_table('intervals.tsv',FIELDS,records.values())
    write_table('allele_counts.tsv.gz',['variant_id','maf','mac','obs','eligible'],(dict(variant_id=v,**c) for v,c in count_map.items()))
    Path('ld_qc.json').write_text(json.dumps(dict(parameters=vars(args),plink_version=version,n_founders=len(selected),n_nonfounders_excluded=nonfounders,n_assigned_missing_from_genotypes=sum(g==args.group and iid not in seen_samples for iid,g in group_map.items()),n_variants=len(positions),n_eligible=sum(c['eligible'] for c in count_map.values()),status_counts={s:sum(r['status']==s for r in records.values()) for s in sorted({r['status'] for r in records.values()})}),indent=2)+'\n')
    log('Wrote ancestry intervals and allele-count diagnostics')

if __name__=='__main__':
    main()

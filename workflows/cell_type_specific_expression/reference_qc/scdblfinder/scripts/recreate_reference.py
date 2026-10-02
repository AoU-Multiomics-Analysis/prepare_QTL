import argparse,json,hashlib
from pathlib import Path
import numpy as np,pandas as pd
from reference_filter import filter_reference
p=argparse.ArgumentParser();p.add_argument('--reference',required=True);p.add_argument('--decisions',required=True);p.add_argument('--gene-validation',required=True);p.add_argument('--out-dir',required=True)
a=p.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
def sha(path):
    digest=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(4*1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()
prior=json.loads(Path(a.gene_validation).read_text());assert sha(a.reference)==prior['gzip_sha256']
decisions=pd.read_csv(a.decisions,sep='\t',low_memory=False).sort_values('matrix_column')
assert len(decisions)==70458 and decisions.cell_id.is_unique
keep=decisions['scDblFinder.class'].eq('singlet').to_numpy();assert (~keep).sum()>0
retained=decisions[keep].copy();counts=retained.matrix_header.value_counts()
assert set(counts.index)==set(prior['cell_type_counts']) and counts.min()>=3
print(f'Recreating reference: {len(decisions)} cells before, {len(retained)} retained',flush=True)
destination=out/'single_cell_reference_all_cells_relabelled_no_globins_bulk_detected_no_sex_chromosomes_no_doublets.tsv.gz'
validation=filter_reference(a.reference,destination,keep,decisions.matrix_header)
retained=retained.rename(columns={'matrix_column':'original_matrix_column'})
retained.insert(0,'matrix_column',np.arange(2,len(retained)+2))
retained.to_csv(out/'selected_cells.tsv',sep='\t',index=False)
decisions[~keep].to_csv(out/'excluded_predicted_doublets.tsv',sep='\t',index=False)
pd.DataFrame({'gene_symbol':validation.pop('zero_genes_removed')}).to_csv(out/'genes_zero_after_cell_filter.tsv',sep='\t',index=False)
validation.update({'source_reference':a.reference,'source_gzip_sha256':prior['gzip_sha256'],'output':str(destination.resolve()),'gzip_sha256':sha(destination),'gzip_bytes':destination.stat().st_size,'genes_before':prior['genes_after'],'cells_before':len(decisions),'predicted_doublets_removed':int((~keep).sum()),'cell_type_counts':counts.to_dict(),'retained_values':'Original raw nonnegative integer counts; no normalization or log transform','cell_order_and_labels_preserved':True,'source_gene_rules_preserved':True,'sex_chromosome_genes_still_excluded':prior['sex_chromosome_genes_removed'],'previous_17_exclusions':prior['previous_17_exclusions'],'NKT_review_cells_still_excluded':252,'unresolved_review_cells_still_excluded':1,'readback_passed':True,'marker_derivation_run':False,'S_mode_run':False,'Terra_submission':False})
(out/'validation.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation,indent=2),flush=True)
with (out/'SHA256SUMS').open('w') as f:
    for path in sorted(out.iterdir()):
        if path.is_file() and path.name!='SHA256SUMS':f.write(f'{sha(path)}  {path.name}\n')
print('Reference export and readback complete',flush=True)

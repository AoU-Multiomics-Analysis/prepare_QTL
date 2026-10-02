import argparse,json
from pathlib import Path
import pandas as pd,numpy as np
from cell_decisions import reference_decisions
p=argparse.ArgumentParser();p.add_argument('--calls',nargs='+',required=True);p.add_argument('--selected',required=True);p.add_argument('--flagged',required=True);p.add_argument('--focused',required=True);p.add_argument('--out-dir',required=True)
a=p.parse_args();out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
calls=pd.concat([pd.read_csv(x,sep='\t') for x in a.calls],ignore_index=True)
assert len(calls)==71101 and calls.cell_id.is_unique and calls['10X_run'].nunique()==22
assert np.isfinite(calls['scDblFinder.score']).all()
calls['predicted_doublet']=calls['scDblFinder.class'].eq('doublet')
calls.to_csv(out/'all_cell_calls.tsv.gz',sep='\t',index=False)
selected=pd.read_csv(a.selected,sep='\t',low_memory=False);flagged=pd.read_csv(a.flagged,sep='\t')
ref=reference_decisions(selected,calls)
ref.to_csv(out/'reference_cell_decisions.tsv',sep='\t',index=False)
flagged_calls=flagged.merge(calls[['cell_id','scDblFinder.score','scDblFinder.class','predicted_doublet']],on='cell_id',how='left',validate='one_to_one')
assert len(flagged_calls)==226 and flagged_calls['scDblFinder.class'].notna().all()
flagged_calls.to_csv(out/'flagged_CD8_226_calls.tsv',sep='\t',index=False)
for keys,name,frame in [(['10X_run'],'capture',calls),(['donor_id'],'donor',calls),(['matrix_header'],'reference_cell_type',ref),(['donor_id','matrix_header'],'reference_donor_cell_type',ref)]:
    table=frame.groupby(keys,dropna=False).agg(cells=('cell_id','size'),doublets=('predicted_doublet','sum')).reset_index()
    table['singlets']=table.cells-table.doublets;table['doublet_fraction']=table.doublets/table.cells
    table.to_csv(out/(name+'_summary.tsv'),sep='\t',index=False)
focused=pd.read_csv(a.focused,sep='\t');cd8=focused[focused.matrix_header.eq('CD8_T')].merge(calls[['cell_id','predicted_doublet']],on='cell_id',validate='one_to_one')
expression=[]
for donor in ['all_donors']+sorted(cd8.donor_id.unique()):
    frame=cd8 if donor=='all_donors' else cd8[cd8.donor_id.eq(donor)]
    for group,subset in [('all_CD8',frame),('predicted_singlet_CD8',frame[~frame.predicted_doublet]),('predicted_doublet_CD8',frame[frame.predicted_doublet])]:
        for gene in ['S100A8','S100A9','S100A12','CD3E','CD8B','LYZ','LST1','CSF3R','FCGR3B','CEACAM8']:
            expression.append({'donor':donor,'group':group,'gene':gene,'cells':len(subset),'mean_CPM':subset[gene+'_CPM'].mean(),'detected_fraction':subset[gene+'_UMI'].gt(0).mean()})
pd.DataFrame(expression).to_csv(out/'CD8_marker_comparison.tsv',sep='\t',index=False)
summary={'scope':'All 22 eligible 10x 3-prime v3 captures; TSP10 remains excluded','method':'scDblFinder 1.14.0','cells_scored':len(calls),'donors':int(calls.donor_id.nunique()),'captures':22,'predicted_doublets_all_scored':int(calls.predicted_doublet.sum()),'reference_cells_before':len(ref),'reference_doublets_removed':int(ref.predicted_doublet.sum()),'reference_cells_retained':int((~ref.predicted_doublet).sum()),'flagged_CD8_total':226,'flagged_CD8_doublets':int(flagged_calls.predicted_doublet.sum()),'flagged_CD8_singlets':int((~flagged_calls.predicted_doublet).sum()),'all_CD8_before':len(cd8),'CD8_doublets_removed':int(cd8.predicted_doublet.sum()),'SoupX_run':False,'CellBender_run':False,'Terra_run':False}
(out/'run_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2),flush=True)

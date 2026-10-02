import argparse, gzip, hashlib, json
from pathlib import Path
import h5py, numpy as np, pandas as pd
from scipy import sparse, io

p=argparse.ArgumentParser()
p.add_argument('--h5ad',required=True);p.add_argument('--selected',required=True);p.add_argument('--flagged',required=True);p.add_argument('--out-dir',required=True);p.add_argument('--capture-plan',required=True)
a=p.parse_args(); plan=pd.read_csv(a.capture_plan,sep='\t'); out=Path(a.out_dir);out.mkdir(parents=True,exist_ok=True)
def column(g,k):
    v=g[k]
    if isinstance(v,h5py.Group):
        cats=v['categories'].asstr()[:];codes=v['codes'][:]
        return np.array([cats[c] if c>=0 else '' for c in codes])
    return v.asstr()[:] if v.dtype.kind in 'OSU' else v[:]
selected=pd.read_csv(a.selected,sep='\t',low_memory=False); flagged=pd.read_csv(a.flagged,sep='\t')
assert selected.cell_id.is_unique and flagged.cell_id.is_unique and len(flagged)==226
with h5py.File(a.h5ad,'r') as f:
    obs=f['obs']; ids=column(obs,'_index')
    meta=pd.DataFrame({'cell_id':ids,'atlas_row':np.arange(len(ids))})
    for key in ['donor_id','assay','10X_run','cell_type','n_genes_by_counts','pct_counts_mt','total_counts']:
        meta[key]=column(obs,key)
    meta=meta[(meta.assay=="10x 3' v3")&(meta.donor_id!='TSP10')].copy()
    assert meta.cell_id.is_unique and len(meta)==71101
    assert len(selected)==70458 and len(flagged)==226
    assert set(selected.cell_id).issubset(set(meta.cell_id))
    assert set(meta['10X_run'])==set(plan.capture) and len(plan)==22
    assert set(flagged.cell_id).issubset(set(selected.cell_id))
    assert (meta.set_index('cell_id').loc[selected.cell_id,'atlas_row'].to_numpy()==selected.atlas_row.to_numpy()).all()
    meta['in_reference']=meta.cell_id.isin(selected.cell_id)
    meta['flagged_CD8_226']=meta.cell_id.isin(flagged.cell_id)
    meta['reference_label']=meta.cell_id.map(selected.set_index('cell_id').matrix_header).fillna('context_only')
    meta.to_csv(out/'all_cells.tsv',sep='\t',index=False)
    features=column(f['raw/var'],'_index'); assert len(set(features))==len(features)
    ptr=f['raw/X/indptr'][:]; raw=f['raw/X']; manifest=[]
    for index,(capture,cells) in enumerate(meta.groupby('10X_run',sort=True)):
        expected=plan.loc[plan.capture.eq(capture)].iloc[0]
        assert len(cells)==int(expected.expected_cells)
        folder=out/capture;folder.mkdir(exist_ok=True)
        rows=cells.atlas_row.to_numpy();blocks=[]
        print(f'Exporting capture {index+1}/22: {capture}, {len(rows)} cells',flush=True)
        for start in range(0,len(rows),300):
            chunk=rows[start:start+300];lo=int(chunk[0]);hi=int(chunk[-1])+1
            left=int(ptr[lo]);right=int(ptr[hi])
            values=raw['data'][left:right]
            if not np.isfinite(values).all() or (values<0).any() or (values!=np.floor(values)).any() or (values>np.iinfo(np.int32).max).any():
                raise ValueError('Raw counts must be finite nonnegative integers')
            block=sparse.csr_matrix((values.astype(np.int32),raw['indices'][left:right],ptr[lo:hi+1]-left),shape=(hi-lo,len(features)))
            blocks.append(block[chunk-lo])
        counts=sparse.vstack(blocks,format='csr').T.tocsr()
        expressed=np.asarray(counts.sum(axis=1)).ravel()>0
        counts=counts[expressed]
        assert (np.asarray(counts.sum(axis=0)).ravel()>0).all()
        with gzip.open(folder/'counts.mtx.gz','wb',compresslevel=1) as target:
            io.mmwrite(target,counts,field='integer')
        pd.DataFrame({'feature_id':features[expressed]}).to_csv(folder/'features.tsv',sep='\t',index=False)
        cells.to_csv(folder/'cells.tsv',sep='\t',index=False)
        manifest.append({'capture':capture,'cells':len(cells),'features':counts.shape[0],'seed':int(expected.seed)})
        del counts,blocks,block
    assert len(manifest)==22
    pd.DataFrame(manifest).to_csv(out/'captures.tsv',sep='\t',index=False)
(out/'provenance.json').write_text(json.dumps({'source_h5ad':a.h5ad,'raw_shape':[len(ids),len(features)],'cells_scored':len(meta),'reference_cells':len(selected),'flagged_CD8':len(flagged),'capture_key':'10X_run','full_raw_counts':True},indent=2))
print('Count export complete',flush=True)

import numpy as np

def reference_decisions(selected,calls):
    if not selected.cell_id.is_unique or not calls.cell_id.is_unique:
        raise ValueError('Cell IDs must be unique')
    if not set(selected.cell_id).issubset(set(calls.cell_id)):
        raise ValueError('Some reference cells have no doublet call')
    if not calls['scDblFinder.class'].isin(['singlet','doublet']).all() or not np.isfinite(calls['scDblFinder.score']).all():
        raise ValueError('Invalid or incomplete doublet calls')
    out=selected.merge(calls[['cell_id','scDblFinder.score','scDblFinder.class']],on='cell_id',how='left',validate='one_to_one').sort_values('matrix_column').reset_index(drop=True)
    if not np.array_equal(out.matrix_column.to_numpy(),np.arange(2,len(out)+2)):
        raise ValueError('Reference matrix columns must be complete and in order')
    out['predicted_doublet']=out['scDblFinder.class'].eq('doublet')
    return out

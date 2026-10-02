import unittest
import pandas as pd
from scripts.cell_decisions import reference_decisions
class CellDecisionsTest(unittest.TestCase):
    def setUp(self):
        self.selected=pd.DataFrame({'cell_id':['a','b','c'],'matrix_column':[2,3,4],'matrix_header':['CD8_T','NK','B']})
        self.calls=pd.DataFrame({'cell_id':['c','a','b','context'],'scDblFinder.class':['singlet','doublet','singlet','doublet'],'scDblFinder.score':[.01,.95,.03,.9]})
    def test_restores_reference_order_and_excludes_context(self):
        out=reference_decisions(self.selected,self.calls)
        self.assertEqual(out.cell_id.tolist(),['a','b','c'])
        self.assertEqual(out.predicted_doublet.tolist(),[True,False,False])
    def test_rejects_missing_calls(self):
        with self.assertRaises(ValueError): reference_decisions(self.selected,self.calls[self.calls.cell_id!='b'])
    def test_rejects_duplicate_calls(self):
        with self.assertRaises(ValueError): reference_decisions(self.selected,pd.concat([self.calls,self.calls.iloc[[0]]]))
if __name__=='__main__': unittest.main()

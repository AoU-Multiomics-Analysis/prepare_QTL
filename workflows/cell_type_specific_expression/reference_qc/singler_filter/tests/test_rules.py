import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from prepare import keep_cell
class Rules(unittest.TestCase):
    def row(self,label,main,fine):return dict(matrix_header=label,Monaco_main_comparison=main,Monaco_fine_comparison=fine)
    def test_cd8_requires_explicit_fine(self):
        self.assertTrue(keep_cell(self.row('CD8_T','disagreement','compatible')))
        self.assertFalse(keep_cell(self.row('CD8_T','broadly_compatible','T_subtype_unresolved')))
        self.assertFalse(keep_cell(self.row('CD8_T','broadly_compatible','uncertain')))
    def test_plasma_both_uncertain_only(self):
        self.assertFalse(keep_cell(self.row('Plasma','uncertain','uncertain')))
        self.assertTrue(keep_cell(self.row('Plasma','broadly_compatible','disagreement')))
        self.assertTrue(keep_cell(self.row('Plasma','uncertain','compatible')))
    def test_other_groups_preserved(self):
        for label in ['B','CD4_T','NK','Monocyte_macrophage','Neutrophil','Erythroid','Platelet']:
            self.assertTrue(keep_cell(self.row(label,'uncertain','uncertain')))
if __name__=='__main__':unittest.main()

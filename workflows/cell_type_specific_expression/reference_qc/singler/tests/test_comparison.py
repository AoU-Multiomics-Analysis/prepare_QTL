import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from comparison import compare

class CoverageTests(unittest.TestCase):
    def test_missing_coverage(self):
        self.assertEqual(compare('Erythroid', 'Neutrophils', 'Neutrophils', 'Monaco_main'), 'not_covered')
        self.assertEqual(compare('Platelet', 'Naive B cells', '', 'Monaco_fine'), 'not_covered')
    def test_broad_T_not_subtype_evidence(self):
        self.assertEqual(compare('CD8_T', 'T_cells', 'T_cells', 'HPCA_main'), 'broadly_compatible')
        self.assertEqual(compare('CD8_T', 'MAIT cells', 'MAIT cells', 'Monaco_fine'), 'T_subtype_unresolved')
    def test_basophil_disagreement(self):
        self.assertEqual(compare('CD8_T', 'Low-density basophils', 'Low-density basophils', 'Monaco_fine'), 'disagreement')
    def test_uncertainty(self):
        self.assertEqual(compare('CD8_T', 'Effector memory CD8 T cells', '', 'Monaco_fine'), 'uncertain')
    def test_plasma_coverage(self):
        self.assertEqual(compare('Plasma', 'B cells', 'B cells', 'Monaco_main'), 'broadly_compatible')
        self.assertEqual(compare('Plasma', 'Plasmablasts', 'Plasmablasts', 'Monaco_fine'), 'compatible')
        self.assertEqual(compare('B', 'Plasmablasts', 'Plasmablasts', 'Monaco_fine'), 'disagreement')

if __name__ == '__main__':
    unittest.main()

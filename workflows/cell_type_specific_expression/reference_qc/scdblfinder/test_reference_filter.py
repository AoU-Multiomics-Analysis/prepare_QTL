import tempfile, unittest, gzip
from pathlib import Path
from scripts.reference_filter import filter_reference

class ReferenceFilterTest(unittest.TestCase):
    def test_preserves_counts_labels_and_order(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            with gzip.open(p/'in.gz','wt') as f:
                f.write('Gene\tCD8_T\tNK\tCD8_T\nA\t1\t7\t3\nB\t0\t9\t0\n')
            result=filter_reference(p/'in.gz',p/'out.gz',[True,False,True],['CD8_T','NK','CD8_T'])
            with gzip.open(p/'out.gz','rt') as f:
                self.assertEqual(f.read(),'Gene\tCD8_T\tCD8_T\nA\t1\t3\n')
            self.assertEqual(result['genes'],1)
            self.assertEqual(result['zero_genes_removed'],['B'])
    def test_rejects_misaligned_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            with gzip.open(p/'in.gz','wt') as f: f.write('Gene\tNK\nA\t1\n')
            with self.assertRaises(ValueError):
                filter_reference(p/'in.gz',p/'out.gz',[True],['CD8_T'])
if __name__=='__main__': unittest.main()

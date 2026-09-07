import csv
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'tools/phenotype_pc_qtl/prepare_pc_scan.py'

class PrepareTests(unittest.TestCase):
    def run_case(self, matrix, ids='exprPC1\nexprPC2\n', chromosome='chrY'):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'cov.tsv').write_text(matrix)
            (p/'ids.txt').write_text(ids)
            r=subprocess.run([sys.executable,str(SCRIPT),'--covariates',str(p/'cov.tsv'),'--phenotype-pc-ids',str(p/'ids.txt'),'--chromosome',chromosome],cwd=p,text=True,capture_output=True)
            files={f.name:f.read_text() for f in p.glob('pc_scan.*')}
            return r,files

    def test_split_preserves_samples_values_and_non_pc_covariates(self):
        r,f=self.run_case('ID\tS3\tS1\tS2\tS4\tS5\tS6\n'
            'genPC1\t1\t2\t3\t4\t5\t6\nexprPC2\t2.0\t-1\t0\t3\t4\t5\n'
            'sex\t0\t1\t0\t1\t0\t1\nexprPC1\t-2\t1e-3\t1\t2\t3\t4\n')
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(f['pc_scan.phenotypes.bed'],'#chr\tstart\tend\tphenotype_id\tS3\tS1\tS2\tS4\tS5\tS6\nchrY\t0\t1\texprPC1\t-2\t1e-3\t1\t2\t3\t4\nchrY\t0\t1\texprPC2\t2.0\t-1\t0\t3\t4\t5\n')
        self.assertEqual(f['pc_scan.covariates.tsv'],'ID\tS3\tS1\tS2\tS4\tS5\tS6\ngenPC1\t1\t2\t3\t4\t5\t6\nsex\t0\t1\t0\t1\t0\t1\n')

    def test_rejects_bad_inputs(self):
        good='ID\tA\tB\tC\tD\ngenPC1\t1\t2\t3\t4\nexprPC1\t4\t2\t1\t3\n'
        cases=[(good,'missing\n'),(good,'exprPC1\nexprPC1\n'),(good.replace('A\tB','A\tA'),'exprPC1\n'),(good.replace('4\t2','NaN\t2'),'exprPC1\n'),(good.replace('4\t2\t1\t3','1\t1\t1\t1'),'exprPC1\n'),(good,'exprPC1\ngenPC1\n')]
        for matrix,ids in cases:
            with self.subTest(matrix=matrix,ids=ids):
                r,_=self.run_case(matrix,ids)
                self.assertNotEqual(r.returncode,0)
                self.assertIn('ERROR:',r.stderr)

    def test_rejects_unlocalized_cloud_uri(self):
        r=subprocess.run([sys.executable,str(SCRIPT),'--covariates','gs://bucket/cov.tsv','--phenotype-pc-ids','missing.txt'],text=True,capture_output=True)
        self.assertNotEqual(r.returncode,0)
        self.assertIn('localization error',r.stderr)

if __name__=='__main__': unittest.main()

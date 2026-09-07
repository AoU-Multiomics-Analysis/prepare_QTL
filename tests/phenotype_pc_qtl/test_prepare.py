from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'tools/phenotype_pc_qtl/prepare_pc_scan.py'

class PrepareTests(unittest.TestCase):
    def run_case(self, matrix):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d); (p/'cov.tsv').write_text(matrix)
            r=subprocess.run([sys.executable,str(SCRIPT),'--covariates',str(p/'cov.tsv'),'--chromosome','chrY'],cwd=p,text=True,capture_output=True)
            return r,{f.name:f.read_text() for f in p.glob('pc_scan.*')}

    def test_auto_selects_only_exact_pc_digit_names(self):
        names=['GENETICPC1','PC2','sex','PC1','PC10','PC1_extra','pc3','PC','myPC4']
        header='ID\tS3\tS1\tS2\tS4\tS5\tS6\tS7\tS8\tS9\tS10\n'
        values='1\t2\t3\t4\t5\t6\t7\t8\t9\t1e-3\n'
        r,f=self.run_case(header+''.join(n+'\t'+values for n in names))
        self.assertEqual(r.returncode,0,r.stderr)
        self.assertEqual(f['pc_scan.phenotypes.bed'],'#chr\tstart\tend\tphenotype_id\tS3\tS1\tS2\tS4\tS5\tS6\tS7\tS8\tS9\tS10\n'+''.join('chrY\t0\t1\t'+n+'\t'+values for n in ['PC2','PC1','PC10']))
        self.assertEqual(f['pc_scan.covariates.tsv'],header+''.join(n+'\t'+values for n in names if n not in ['PC2','PC1','PC10']))

    def test_rejects_bad_inputs(self):
        good='ID\tA\tB\tC\tD\nGENETICPC1\t1\t2\t3\t4\nPC1\t4\t2\t1\t3\n'
        cases=[good.replace('\nPC1','\nexprPC1'),good+'PC1\t1\t2\t3\t4\n',good.replace('A\tB','A\tA'),good.replace('4\t2','NaN\t2'),good.replace('4\t2\t1\t3','1\t1\t1\t1'),good.replace('GENETICPC1','PC2')]
        for matrix in cases:
            with self.subTest(matrix=matrix):
                r,_=self.run_case(matrix)
                self.assertNotEqual(r.returncode,0)
                self.assertIn('ERROR:',r.stderr)

    def test_rejects_unlocalized_cloud_uri(self):
        r=subprocess.run([sys.executable,str(SCRIPT),'--covariates','gs://bucket/cov.tsv'],text=True,capture_output=True)
        self.assertNotEqual(r.returncode,0)
        self.assertIn('localization error',r.stderr)

if __name__=='__main__': unittest.main()

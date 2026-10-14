from pathlib import Path
import csv
import gzip
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / 'tools/regenie_export/export_regenie.py'
BED = '#chr\tstart\tend\tphenotype_id\t001\tS2\tS3\nchr1\t0\t1\tENSG1.2\t1.23000\tNA\t-2e-08\nchr2\t3\t4\tchr2:3:4:clu_1\t3\t4\t5\nchr2\t8\t9\tprotein,1\t6\t7\t8\n'


class ExportTests(unittest.TestCase):
    def run_case(self, bed=BED, cov=None, orientation='sample', size=2):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "donor's $(touch UNSAFE) bed.gz"
            with gzip.open(source, 'wt') as stream:
                stream.write(bed)
            args = [sys.executable, str(SCRIPT), '--bed', str(source),
                    '--output-prefix', 'qtl', '--phenotypes-per-file', str(size)]
            if cov is not None:
                path = root / 'covariates.tsv'
                path.write_text(cov)
                args += ['--covariates', str(path), '--covariate-format', orientation]
            result = subprocess.run(args, cwd=root, capture_output=True, text=True)
            files = {p.name: p.read_text() for p in (root / 'outputs').glob('*.tsv')}
            self.assertFalse((root / 'UNSAFE').exists())
            return result, files

    def test_transpose_exact_values_and_split_all_traits(self):
        result, files = self.run_case()
        self.assertEqual(result.returncode, 0, result.stderr)
        traits = list(csv.DictReader(files['qtl.traits.tsv'].splitlines(), delimiter='\t'))
        self.assertEqual([r['source_id'] for r in traits], ['ENSG1.2', 'chr2:3:4:clu_1', 'protein,1'])
        self.assertEqual([r['chromosome'] for r in traits], ['chr1', 'chr2', 'chr2'])
        names = [r['phenotype_id'] for r in traits]
        for name in names:
            self.assertRegex(name, r'^pheno_[a-f0-9]{24}$')
        self.assertEqual(files['qtl.phenotypes.000001.tsv'],
                         'FID\tIID\t' + '\t'.join(names[:2]) + '\n0\t001\t1.23000\t3\n0\tS2\tNA\t4\n0\tS3\t-2e-08\t5\n')
        self.assertEqual(files['qtl.phenotypes.000002.tsv'],
                         'FID\tIID\t' + names[2] + '\n0\t001\t6\n0\tS2\t7\n0\tS3\t8\n')
        batches = list(csv.DictReader(files['qtl.batches.tsv'].splitlines(), delimiter='\t'))
        self.assertEqual(batches[0]['phenotype_columns'], ','.join(names[:2]))
        self.assertNotIn('qtl.covariates.tsv', files)
        other, reordered = self.run_case(bed=BED.splitlines()[0] + '\n' + '\n'.join(reversed(BED.splitlines()[1:])) + '\n')
        self.assertEqual(other.returncode, 0, other.stderr)
        mapping = {r['source_id']: r['phenotype_id'] for r in csv.DictReader(reordered['qtl.traits.tsv'].splitlines(), delimiter='\t')}
        self.assertEqual(mapping, dict(zip([r['source_id'] for r in traits], names)))

    def test_covariate_overlap_follows_bed_order_and_reports_exclusions(self):
        cov = 'sample_id\tage\tbatch\nS3\t38.00\tB\nextra\t42\tC\n001\t40\tA\n'
        result, files = self.run_case(cov=cov)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(files['qtl.covariates.tsv'], 'FID\tIID\tage\tbatch\n0\t001\t40\tA\n0\tS3\t38.00\tB\n')
        self.assertEqual(files['qtl.excluded_samples.tsv'], 'IID\treason\nS2\tmissing_covariate_sample\n')
        rows = files['qtl.phenotypes.000001.tsv'].splitlines()
        self.assertEqual([r.split('\t')[1] for r in rows[1:]], ['001', 'S3'])

    def test_tensor_covariates_transpose_without_rounding(self):
        cov = 'ID\tS3\t001\tS2\nGENETICPC1\t1e-9\t2.000\t-3\nPC1\t4\t5\t6\n'
        result, files = self.run_case(cov=cov, orientation='tensor')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(files['qtl.covariates.tsv'], 'FID\tIID\tGENETICPC1\tPC1\n0\t001\t2.000\t5\n0\tS2\t-3\t6\n0\tS3\t1e-9\t4\n')

    def test_missing_phenotypes_are_na(self):
        result, files = self.run_case(bed=BED.replace('1.23000', 'NaN').replace('-2e-08', '.'))
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('0\t001\tNA\t3', files['qtl.phenotypes.000001.tsv'])
        self.assertIn('0\tS3\tNA\t5', files['qtl.phenotypes.000001.tsv'])

    def test_accepts_numeric_padding_from_existing_r_covariate_merge(self):
        # t(data.frame) formats numeric columns before mixing them with sample IDs.
        cov = 'ID\t001\tS2\tS3\nage\t20\t60\t 9\nPC1\t-1.21000\t 0.00013\t 3.91000\n'
        result, files = self.run_case(cov=cov, orientation='tensor')
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(files['qtl.covariates.tsv'], 'FID\tIID\tage\tPC1\n0\t001\t20\t-1.21000\n0\tS2\t60\t0.00013\n0\tS3\t9\t3.91000\n')

    def test_rejects_invalid_beds(self):
        for bed in [BED.replace('001\tS2', '001\t001'), BED.replace('protein,1', 'ENSG1.2'),
                    BED.replace('1.23000', 'Inf'), BED.replace('1.23000', 'oops'),
                    BED.replace('1.23000\tNA', '1.23000'), BED.splitlines()[0] + '\n',
                    BED.replace('001\tS2', 'bad id\tS2')]:
            with self.subTest(bed=bed):
                result, _ = self.run_case(bed=bed)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('ERROR:', result.stderr)

    def test_rejects_invalid_covariates(self):
        for cov in ['ID\tage\n001\tNA\n', 'ID\tage\n001\tInf\n',
                    'ID\tage\n001\t1\n001\t2\n', 'ID\tage\tage\n001\t1\t2\n',
                    'ID\tage\nother\t1\n', 'ID\tage\n001\t\n']:
            with self.subTest(cov=cov):
                result, _ = self.run_case(cov=cov)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('ERROR:', result.stderr)

    def test_rejects_cloud_uri_and_bad_batch_size(self):
        result = subprocess.run([sys.executable, str(SCRIPT), '--bed', 'gs://bucket/bed.gz'], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('localization error', result.stderr)
        result, _ = self.run_case(size=0)
        self.assertNotEqual(result.returncode, 0)


if __name__ == '__main__':
    unittest.main()

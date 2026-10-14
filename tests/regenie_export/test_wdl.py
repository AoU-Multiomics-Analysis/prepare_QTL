from pathlib import Path
import subprocess
import tempfile
import unittest

import WDL

ROOT = Path(__file__).resolve().parents[2]
DESCRIPTOR = ROOT / 'workflows/common/ExportRegenie.wdl'


class WdlTests(unittest.TestCase):
    def test_embedded_program_matches_source(self):
        source = DESCRIPTOR.read_text().split('# BEGIN REGENIE EXPORT\n', 1)[1].split('# END REGENIE EXPORT', 1)[0]
        self.assertEqual(source, (ROOT / 'tools/regenie_export/export_regenie.py').read_text())

    def test_file_localization_and_optional_covariates_in_real_command(self):
        task = WDL.load(str(DESCRIPTOR)).tasks[0]
        for present in (False, True):
            with self.subTest(present=present), tempfile.TemporaryDirectory() as directory:
                root = Path(directory)
                bed = root / "new donor's $literal $(touch UNSAFE).bed"
                bed.write_text('#chr\tstart\tend\tID\t001\tS2\nchr1\t1\t2\tx:y\t1\t2\n')
                cov = root / "cov donor's `touch UNSAFE`.tsv"
                cov.write_text('ID\t001\tS2\nage\t20\t30\n')
                stdlib = WDL.StdLib.Base('1.0')
                env = WDL.Env.Bindings()
                for declaration in task.inputs:
                    if declaration.name == 'Bed':
                        value = WDL.Value.File('gs://incoming/bed.bed')
                    elif declaration.name == 'Covariates':
                        value = WDL.Value.File('gs://generated/merge/cov.tsv') if present else WDL.Value.Null()
                    elif declaration.name == 'CovariateFormat':
                        value = WDL.Value.String('tensor')
                    else:
                        value = declaration.expr.eval(env, stdlib)
                    env = env.bind(declaration.name, value.coerce(declaration.type))
                # Evaluate declarations with cloud paths, then localize only File values.
                for declaration in task.postinputs:
                    env = env.bind(declaration.name, declaration.expr.eval(env, stdlib).coerce(declaration.type))
                env = WDL.Value.rewrite_env_paths(env, lambda value: str(bed if value.value.endswith('.bed') else cov))
                command = task.command.eval(env, stdlib).value
                self.assertNotIn('gs://', command)
                result = subprocess.run(['bash'], input=command, cwd=root, text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertFalse((root / 'UNSAFE').exists())
                self.assertEqual((root / 'outputs/qtl.covariates.tsv').exists(), present)
                self.assertIn('Complete traits=1 samples=2', (root / 'export.log').read_text())

    def test_public_preparation_workflows_expose_regenie_batches(self):
        for path in ['expression/prepare_eQTL.wdl', 'proteomics/prepare_pQTL.wdl',
                     'splicing/prepare_sQTL.wdl', 'methylation/prepare_mQTL.wdl',
                     'methylation/AggregateMethylationCohort.wdl',
                     'methylation/merge_methylation.wdl',
                     'cell_type_specific_expression/prepare_cell_type_eQTL.wdl',
                     'phenotype_pc_qtl/phenotype_pc_qtl.wdl']:
            with self.subTest(path=path):
                document = WDL.load(str(ROOT / 'workflows' / path))
                output_types = {d.name: str(d.type) for d in document.workflow.outputs}
                batch_names = [n for n in output_types if n.endswith('RegeniePhenotypes')]
                self.assertTrue(batch_names, path)
                self.assertTrue(all('File' in output_types[n] and output_types[n].startswith('Array[') for n in batch_names))


if __name__ == '__main__':
    unittest.main()

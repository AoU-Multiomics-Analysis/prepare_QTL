"""Check executable entry points after migration into stage directories."""
from pathlib import Path
import importlib.util
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]


class StageLayoutTest(unittest.TestCase):
    def test_release_smoke_selects_candidate_script_not_trusted_layout(self):
        spec = importlib.util.spec_from_file_location('rnaseqc_smoke', ROOT / 'tests/rnaseqc2_aggregation/smoke_container.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertTrue(callable(getattr(module, 'script_for_wdl', None)))
        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / 'candidate.wdl'
            for path in ('/opt/prepare_qtl/scripts/expression/merge_rnaseqc.py',
                         '/opt/prepare_qtl/scripts/expression/rnaseqc/merge_rnaseqc.py'):
                candidate.write_text('python3 ' + path + ' validate-manifest\n')
                self.assertEqual(module.script_for_wdl(candidate), path)
            candidate.write_text('python3 /unexpected/tool.py\n')
            with self.assertRaises(ValueError):
                module.script_for_wdl(candidate)

    def test_rnaseqc_cli_from_stage_directory(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / 'scripts/expression/rnaseqc/merge_rnaseqc.py'), '--help'],
            capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('validate-manifest', result.stdout)


if __name__ == '__main__':
    unittest.main()

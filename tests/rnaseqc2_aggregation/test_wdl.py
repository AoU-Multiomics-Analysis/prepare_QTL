import json
import shutil
import subprocess
import unittest
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[2]
WDL = PROJECT_DIR / "workflows" / "expression" / "rnaseqc2_aggregate_batched.wdl"


class WdlSourceTest(unittest.TestCase):
    def test_optional_insert_size_output_does_not_create_an_empty_gzip(self) -> None:
        source = WDL.read_text()
        self.assertNotIn("touch ", source)
        self.assertIn(
            "Array[File] insert_size_hists = " 'glob("*.insert_size_hists.txt.gz")',
            source,
        )
        self.assertIn('if [[ "~{include_insert_sizes}" == "true" ]]', source)

    def test_output_prefix_reaches_commands_through_a_file(self) -> None:
        source = WDL.read_text()
        self.assertIn("File prefix_file = write_lines([prefix])", source)
        self.assertNotIn('"~{prefix}', source)


@unittest.skipUnless(shutil.which("miniwdl"), "miniwdl is not installed")
class WdlTest(unittest.TestCase):
    def test_workflow_passes_static_type_checking(self) -> None:
        result = subprocess.run(
            ["miniwdl", "check", str(WDL)],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_workflow_uses_one_combined_sample_manifest(self) -> None:
        result = subprocess.run(
            ["miniwdl", "input_template", str(WDL)],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        inputs = json.loads(result.stdout)
        self.assertIn("rnaseqc2_aggregate_batched_workflow.sample_manifest", inputs)
        self.assertNotIn("rnaseqc2_aggregate_batched_workflow.tpm_gcts_list", inputs)
        self.assertNotIn("rnaseqc2_aggregate_batched_workflow.sample_ids_list", inputs)


if __name__ == "__main__":
    unittest.main()

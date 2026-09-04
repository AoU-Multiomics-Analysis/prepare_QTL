import json
import shutil
import subprocess
import unittest
from pathlib import Path

try:
    import WDL as miniwdl
except ImportError:
    miniwdl = None

PROJECT_DIR = Path(__file__).resolve().parents[2]
WDL = PROJECT_DIR / "workflows" / "expression" / "rnaseqc2_aggregate_batched.wdl"


def workflow_file_writes(node) -> list[str]:
    """Inspect workflow syntax only; a call's children exclude its task body."""
    writes = []
    if isinstance(node, miniwdl.Expr.Apply) and node.function_name.startswith("write_"):
        writes.append(f"{node.function_name} at line {node.pos.line}")
    for child in node.children:
        writes.extend(workflow_file_writes(child))
    return writes


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
        self.assertNotIn('"~{prefix}', source)


@unittest.skipUnless(miniwdl, "miniwdl is not installed")
class TerraFileScopeTest(unittest.TestCase):
    def test_workflow_does_not_need_engine_side_file_writes(self) -> None:
        document = miniwdl.load(str(WDL))
        self.assertEqual(
            workflow_file_writes(document.workflow), [],
            "Terra file creation must occur in task scope, not workflow scope",
        )

    def test_scope_check_catches_nested_workflow_file_writes(self) -> None:
        cases = [
            'input { File f = write_lines(["sample"]) }',
            'File f = write_lines(["sample"])',
            'call sink { input: f = write_tsv([["sample"]]) }',
            'output { File f = write_map({"sample": "path"}) }',
            'scatter (s in ["sample"]) { '
            'if (true) { call sink { input: f = write_json([s]) } } }',
            'scatter (s in read_lines(write_lines(["sample"]))) { String x = s }',
            'if (size(write_lines(["sample"])) > 0) { String x = "sample" }',
        ]
        for body in cases:
            with self.subTest(body=body):
                document = miniwdl.parse_document(
                    'version 1.0\n'
                    'task sink { input { File f } command { echo "~{f}" } }\n'
                    f'workflow example {{ {body} }}'
                )
                self.assertEqual(len(workflow_file_writes(document.workflow)), 1)

    def test_scope_check_allows_task_file_writes(self) -> None:
        document = miniwdl.parse_document(
            'version 1.0\n'
            'task sink { input { String s } File f = write_lines([s]) '
            'command { cat "~{f}" "~{write_json([s])}" } }\n'
            'workflow example { call sink { input: s = "sample" } }'
        )
        self.assertEqual(workflow_file_writes(document.workflow), [])


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

    def test_workflow_does_not_require_a_script_upload(self) -> None:
        result = subprocess.run(
            ["miniwdl", "input_template", str(WDL)],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertEqual(
            set(json.loads(result.stdout)),
            {
                "rnaseqc2_aggregate_batched_workflow.sample_manifest",
                "rnaseqc2_aggregate_batched_workflow.prefix",
                "rnaseqc2_aggregate_batched_workflow.merge_disk_space_gb",
            },
        )


if __name__ == "__main__":
    unittest.main()

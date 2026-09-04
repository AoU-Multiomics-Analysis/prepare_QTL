import gzip
import json
import os
import shlex
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

try:
    import WDL as miniwdl
except ImportError:
    miniwdl = None

PROJECT_DIR = Path(__file__).resolve().parents[2]
WDL = PROJECT_DIR / "workflows" / "expression" / "rnaseqc2_aggregate_batched.wdl"
SCRIPT = PROJECT_DIR / "scripts" / "expression" / "merge_rnaseqc.py"


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


@unittest.skipUnless(miniwdl, "miniwdl is not installed")
class BatchCommandTest(unittest.TestCase):
    def run_batch(self, temp: Path, merge_exons: bool, mismatch: bool):
        """Run the WDL command; replace only cloud transfers with fixture copies."""
        fixtures = temp / "fixtures"
        fixtures.mkdir()
        manifest_rows = []
        for sample in ("sample_1", "sample_2"):
            sample_dir = fixtures / sample
            sample_dir.mkdir()
            gct = "#1.2\n1\t1\nName\tDescription\toriginal\ngene_1\tGene one\t5\n"
            for kind in ("tpm", "count", "exon"):
                text = gct
                if kind == "exon" and sample == "sample_2" and mismatch:
                    text = gct.replace("1\t1\n", "2\t1\n", 1)
                    text += "gene_2\tGene two\t7\n"
                (sample_dir / f"{kind}.gct").write_text(text)
            (sample_dir / "metrics.tsv").write_text(
                "Sample\toriginal\nTotal Reads\t10\n"
            )
            manifest_rows.append("\t".join([sample] + [
                f"gs://fixtures/{sample}/{name}" for name in
                ("tpm.gct", "count.gct", "exon.gct", "metrics.tsv")
            ]))
        manifest = temp / "samples.tsv"
        manifest.write_text(
            "sample_id\ttpm_gct\tcount_gct\texon_count_gct\tmetrics_tsv\n"
            + "\n".join(manifest_rows) + "\n"
        )
        prefix_file = temp / "prefix.txt"
        prefix_file.write_text("cohort\n")
        bindings = miniwdl.Env.Bindings()
        values = {
            "sample_manifest": miniwdl.Value.File(str(manifest)),
            "prefix_file": miniwdl.Value.File(str(prefix_file)),
            "batch_index": miniwdl.Value.Int(0),
            "batch_size": miniwdl.Value.Int(100),
            "num_threads": miniwdl.Value.Int(2),
            "include_insert_sizes": miniwdl.Value.Boolean(False),
            "merge_exons": miniwdl.Value.Boolean(merge_exons),
        }
        for key, value in values.items():
            bindings = bindings.bind(key, value)
        task = next(task for task in miniwdl.load(str(WDL)).tasks
                    if task.name == "aggregate_rnaseqc_batch")
        command = task.command.eval(bindings, miniwdl.StdLib.Base("1.0")).value
        command = command.replace(
            "/opt/prepare_qtl/scripts/expression/merge_rnaseqc.py",
            shlex.quote(str(SCRIPT)),
        )
        transfer_stub = r'''
gsutil() {
    [[ "$#" == 3 && "$1" == cp && "$2" == gs://fixtures/* ]] || return 99
    printf '%s\n' "$2" >> "$FIXTURE_ROOT/transfers-observed.txt"
    cp "$FIXTURE_ROOT/${2#gs://fixtures/}" "$3"
}
export -f gsutil
'''
        return subprocess.run(
            ["bash"], input=transfer_stub + command, cwd=temp,
            env={**os.environ, "FIXTURE_ROOT": str(fixtures)},
            capture_output=True, text=True, check=False,
        )

    def test_batch_skips_incompatible_exons_and_keeps_qc(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            temp = Path(directory)
            result = self.run_batch(temp, merge_exons=False, mismatch=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertFalse(list(temp.glob("*.exon_reads.gct.gz")))
            self.assertFalse(list((temp / "individual_outputs").glob("*.exon_reads.*")))
            observed = (temp / "fixtures/transfers-observed.txt").read_text()
            self.assertNotIn("exon.gct", observed)
            for kind in ("gene_tpm", "gene_reads"):
                with gzip.open(temp / f"cohort.batch_0.{kind}.gct.gz", "rt") as handle:
                    self.assertEqual(handle.read(),
                                     "#1.2\n1\t2\nName\tDescription\tsample_1\tsample_2\n"
                                     "gene_1\tGene one\t5\t5\n")
            with gzip.open(temp / "cohort.batch_0.metrics.txt.gz", "rt") as handle:
                self.assertEqual(handle.read(),
                                 "sample_id\tTotal Reads\nsample_1\t10\nsample_2\t10\n")

    def test_batch_merges_compatible_exons_when_enabled(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            temp = Path(directory)
            result = self.run_batch(temp, merge_exons=True, mismatch=False)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            with gzip.open(temp / "cohort.batch_0.exon_reads.gct.gz", "rt") as handle:
                self.assertEqual(handle.read(),
                                 "#1.2\n1\t2\nName\tDescription\tsample_1\tsample_2\n"
                                 "gene_1\tGene one\t5\t5\n")

    def test_batch_rejects_incompatible_exons_when_enabled(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            temp = Path(directory)
            result = self.run_batch(temp, merge_exons=True, mismatch=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("GCT row count differs", result.stderr)
            self.assertFalse(list(temp.glob("*.exon_reads.gct.gz")))


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
                "rnaseqc2_aggregate_batched_workflow.merge_exons",
            },
        )


if __name__ == "__main__":
    unittest.main()

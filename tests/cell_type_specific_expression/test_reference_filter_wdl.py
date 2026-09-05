"""Typed contracts for reference preparation and post-export BED filtering."""
import gzip
import hashlib
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import WDL


ROOT = Path(__file__).resolve().parents[2]


class TaskStdLib(WDL.StdLib.Base):
    def _virtualize_filename(self, filename):
        return filename


class CloudGeneratedFileStdLib(TaskStdLib):
    """Model write_lines returning a cloud File, not an already-local filename."""
    def __init__(self, directory):
        super().__init__("1.0", write_dir=directory)
        self.generated_paths = {}

    def _virtualize_filename(self, filename):
        uri = "gs://test-bucket/call-FilterCellTypeBeds/" + Path(filename).name
        self.generated_paths[uri] = filename
        return uri


def render_after_localization(task, local_inputs, directory):
    """Evaluate declarations with cloud Files, then localize inputs for the command.

    Serialized file contents are deliberately NOT rewritten. This models the
    boundary which local-input-only tests miss; it is not a Terra execution.
    """
    paths = {}

    def cloud_path(value):
        uri = "gs://test-bucket/inputs" + value.value
        paths[uri] = value.value
        return uri

    env = WDL.Env.Bindings()
    for binding in local_inputs:
        cloud_value = WDL.Value.rewrite_paths(binding.value, cloud_path)
        env = env.bind(binding.name, cloud_value)
    stdlib = TaskStdLib("1.0", write_dir=directory)
    for decl in task.postinputs:
        env = env.bind(decl.name, decl.expr.eval(env, stdlib))
    # Cromwell's WomFileMapper recursively maps Files in objects, arrays and
    # optionals too. It cannot map strings embedded in a previously written file.
    localized = WDL.Env.Bindings()
    for binding in env:
        localized = localized.bind(binding.name, WDL.Value.rewrite_paths(
            binding.value, lambda value: paths.get(value.value, value.value)
        ))
    return task.command.eval(localized, stdlib).value


class ReferenceFilterWdlTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.deconvolution = WDL.load(str(
            ROOT / "workflows/cell_type_specific_expression/deconvolution.wdl"
        ))
        cls.qtl = WDL.load(str(
            ROOT / "workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl"
        ))

    def workflow_call(self, document, name):
        def calls(nodes):
            for node in nodes:
                if isinstance(node, WDL.Tree.Call):
                    yield node
                elif isinstance(node, (WDL.Tree.Conditional, WDL.Tree.Scatter)):
                    yield from calls(node.body)
        return next(node for node in calls(document.workflow.body) if node.name == name)

    def test_deconvolution_prepares_optional_reference_and_always_filters_beds(self):
        inputs = {decl.name: decl for decl in self.deconvolution.workflow.inputs}
        self.assertIsInstance(inputs["haemopedia_counts"].type, WDL.Type.File)
        self.assertTrue(inputs["haemopedia_counts"].type.optional)
        self.assertEqual(str(inputs["reference_min_mean_log2_cpm1"].type), "Float")
        self.assertIsInstance(inputs["reference_residual_cutoff"].type, WDL.Type.Float)
        self.assertTrue(inputs["reference_residual_cutoff"].type.optional)

        prepare = self.workflow_call(self.deconvolution, "PrepareHaemopedia")
        conditional_calls = {
            node.name for block in self.deconvolution.workflow.body
            if isinstance(block, WDL.Tree.Conditional)
            for node in block.body if isinstance(node, WDL.Tree.Call)
        }
        self.assertIn(prepare.name, conditional_calls)
        filter_call = self.workflow_call(self.deconvolution, "FilterCellTypeBeds")
        self.assertNotIn(filter_call.name, conditional_calls)
        self.assertEqual(str(filter_call.inputs["cell_type_beds"]),
                         "ExportTcaBeds.cell_type_beds")
        self.assertEqual(str(filter_call.inputs["cell_type_bed_inventory"]),
                         "ExportTcaBeds.cell_type_bed_inventory")

    def test_filtered_outputs_are_public_and_route_to_eqtl_scatter(self):
        outputs = {decl.name: str(decl.expr) for decl in self.deconvolution.workflow.outputs}
        required = {
            "filtered_cell_type_beds", "filtered_cell_type_bed_inventory",
            "negative_expression_summary", "reference_gene_comparison",
            "reference_filter_metrics", "reference_filter_plots", "reference_filter_log",
            "haemopedia_reference_summary", "haemopedia_reference_samples",
            "haemopedia_reference_metadata",
        }
        self.assertTrue(required.issubset(outputs))

        scatter = self.workflow_call(self.qtl, "PrepareScatterInputs")
        self.assertEqual(str(scatter.inputs["cell_type_beds"]),
                         "CellTypeDeconvolution.filtered_cell_type_beds")
        self.assertEqual(str(scatter.inputs["cell_type_bed_inventory"]),
                         "CellTypeDeconvolution.filtered_cell_type_bed_inventory")

    def test_filter_task_has_typed_inputs_and_logging(self):
        task_path = ROOT / "workflows/cell_type_specific_expression/tasks/reference_filter.wdl"
        self.assertTrue(task_path.exists())
        document = WDL.load(str(task_path))
        tasks = {task.name: task for task in document.tasks}
        self.assertEqual(set(tasks), {"PrepareHaemopedia", "FilterCellTypeBeds"})
        filter_task = tasks["FilterCellTypeBeds"]
        self.assertEqual(str({decl.name: decl.type for decl in filter_task.inputs}["reference_summary"]),
                         "File?")
        command = str(filter_task.command)
        self.assertIn("filter_cell_type_beds.R", command)
        self.assertIn("status=failed", command)
        self.assertIn("completion_time", command)

    def test_prepare_command_passes_localized_counts_as_literal_shell_data(self):
        document = WDL.load(str(
            ROOT / "workflows/cell_type_specific_expression/tasks/reference_filter.wdl"
        ))
        task = next(task for task in document.tasks if task.name == "PrepareHaemopedia")
        with tempfile.TemporaryDirectory() as directory:
            counts = Path(directory) / "counts ' $(touch unexpected_side_effect) $literal.tsv"
            counts.touch()
            env = WDL.Env.Bindings().bind("counts", WDL.Value.File(str(counts)))
            command = render_after_localization(task, env, directory)
            # Exercise the real command's path setup without needing R on the
            # WDL-only CI host. The integration tests below also run the R CLI.
            setup = command.split("Rscript ", 1)[0]
            result = subprocess.run(
                ["bash", "-c", setup + '\nprintf "%s\\n" "$counts_path"'],
                cwd=directory, text=True, capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertEqual(result.stdout.splitlines()[-1], str(counts))
            self.assertFalse((Path(directory) / "unexpected_side_effect").exists())

    def test_filter_command_passes_localized_arguments_and_optional_values(self):
        document = WDL.load(str(
            ROOT / "workflows/cell_type_specific_expression/tasks/reference_filter.wdl"
        ))
        task = next(task for task in document.tasks if task.name == "FilterCellTypeBeds")
        for use_reference in (False, True):
            with self.subTest(reference=use_reference), tempfile.TemporaryDirectory() as directory:
                inventory = directory + "/inventory's.tsv"
                beds = [directory + "/first/b_cells.bed.gz",
                        directory + '/second $literal "quoted"/cd4_t_cells.bed.gz']
                reference = directory + "/reference's.tsv.gz" if use_reference else None
                env = WDL.Env.Bindings().bind(
                    "cell_type_bed_inventory", WDL.Value.File(inventory)
                ).bind(
                    "cell_type_beds", WDL.Value.Array(WDL.Type.File(), [WDL.Value.File(p) for p in beds])
                ).bind(
                    "reference_summary", WDL.Value.File(reference) if reference else WDL.Value.Null()
                ).bind("min_mean_log2_cpm1", WDL.Value.Float(0.01)).bind(
                    "residual_cutoff", WDL.Value.Float(3.5) if use_reference else WDL.Value.Null()
                )
                command = render_after_localization(task, env, directory)
                stub = """Rscript() {
                  printf '%s\\0' "$@" > captured_args
                  printf 'cell_group\\n' > outputs/filtered_inventory.tsv
                }
                """
                result = subprocess.run(["bash", "-c", stub + command], cwd=directory,
                                        text=True, capture_output=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                tokens = (Path(directory) / "captured_args").read_bytes().decode().split(chr(0))[:-1]
                self.assertIn("--inventory", tokens)
                self.assertEqual(tokens[tokens.index("--inventory") + 1], inventory)
                self.assertEqual(float(tokens[tokens.index("--min-mean-log2-cpm1") + 1]), 0.01)
                bed_list = Path(directory) / tokens[tokens.index("--bed-list") + 1]
                self.assertEqual(bed_list.read_text().splitlines(), beds)
                if use_reference:
                    self.assertEqual(tokens[tokens.index("--reference-summary") + 1], reference)
                    self.assertEqual(float(tokens[tokens.index("--residual-cutoff") + 1]), 3.5)
                else:
                    self.assertNotIn("--reference-summary", tokens)
                    self.assertNotIn("--residual-cutoff", tokens)
                cloud_stdlib = CloudGeneratedFileStdLib(directory)
                argument = next(part for part in task.command.parts
                    if isinstance(part, WDL.Expr.Placeholder) and "write_lines" in str(part.expr))
                generated = argument.expr.eval(env, cloud_stdlib)
                self.assertIsInstance(generated, WDL.Value.File,
                    "Generated BED list must remain File-typed for Cromwell localization")
                localized = WDL.Value.rewrite_paths(generated,
                    lambda value: cloud_stdlib.generated_paths[value.value])
                self.assertEqual(Path(localized.value).read_text().splitlines(), beds)

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is required")
    def test_rendered_tasks_use_localized_paths_and_return_filter_outputs(self):
        self.run_rendered_tasks(use_reference=True)

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is required")
    def test_rendered_filter_without_reference_preserves_bed_values(self):
        self.run_rendered_tasks(use_reference=False)

    def run_rendered_tasks(self, use_reference):
        dependency_check = subprocess.run(
            ["Rscript", "-e", 'quit(status=as.integer(!requireNamespace("edgeR", quietly=TRUE)))'],
            capture_output=True,
        )
        if dependency_check.returncode != 0:
            self.skipTest("edgeR is unavailable")

        task_document = WDL.load(str(
            ROOT / "workflows/cell_type_specific_expression/tasks/reference_filter.wdl"
        ))
        tasks = {task.name: task for task in task_document.tasks}
        script_root = ROOT / "scripts/cell_type_specific_expression"
        with tempfile.TemporaryDirectory(prefix="reference-filter-task-") as directory:
            work = Path(directory)
            (work / "pipeline").symlink_to(script_root, target_is_directory=True)
            counts = work / "counts with apostrophe's $(touch unexpected_side_effect).tsv"
            populations = [
                "NveB", "MemB", "CD4T", "CD8T", "NK", "Mono",
                "MonoNonClassical", "Neut", "Eo", "myDC", "myDC123", "pDC",
            ]
            counts.write_text("gene_id\t" + "\t".join(f"{x}.1" for x in populations) + "\n" +
                              "ENSG000001\t" + "\t".join(["10"] * 12) + "\n" +
                              "ENSG000002\t" + "\t".join(["20"] * 12) + "\n" +
                              "ENSG000003\t" + "\t".join(["40"] * 12) + "\n")
            prepare_env = WDL.Env.Bindings().bind("counts", WDL.Value.File(str(counts)))
            prepare_command = render_after_localization(
                tasks["PrepareHaemopedia"], prepare_env, directory
            ).replace(
                "/opt/prepare_qtl/scripts/cell_type_specific_expression/prepare_haemopedia.R",
                "pipeline/prepare_haemopedia.R",
            )
            process_env = dict(os.environ, CELL_TYPE_SPECIFIC_EXPRESSION_ROOT=str(script_root))
            result = subprocess.run(["bash", "-c", prepare_command], cwd=work,
                                    env=process_env, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertFalse((work / "unexpected_side_effect").exists())
            reference_summary = work / "outputs/reference_summary.tsv.gz"
            self.assertTrue(reference_summary.exists())

            bed = work / "b_cells.bed.gz"
            with gzip.open(bed, "wt") as handle:
                handle.write("#chr\tstart\tend\tgene_id\ts1\ts2\n")
                handle.write("chr1\t0\t1\tENSG000001\t1\t2\n")
                handle.write("chr1\t1\t2\tENSG000002\t2\t3\n")
                handle.write("chr1\t2\t3\tENSG000003\t4\t5\n")
            second_bed = work / "cd4_t_cells.bed.gz"
            shutil.copyfile(bed, second_bed)
            inventory = work / "inventory's.tsv"
            bed_sha256 = hashlib.sha256(bed.read_bytes()).hexdigest()
            inventory.write_text(
                "logical_name\tpath\tsha256\tn_genes\tn_samples\tscale\tcell_group\tslug\n"
                f"b_cells\t{bed.name}\t{bed_sha256}\t3\t2\tcpm\tB cells\tb_cells\n"
                f"cd4_t_cells\t{second_bed.name}\t{bed_sha256}\t3\t2\tcpm\tCD4 T cells\tcd4_t_cells\n"
            )
            filter_env = WDL.Env.Bindings().bind(
                "cell_type_bed_inventory", WDL.Value.File(str(inventory))
            ).bind(
                "cell_type_beds", WDL.Value.Array(WDL.Type.File(), [
                    WDL.Value.File(str(bed)), WDL.Value.File(str(second_bed))
                ])
            ).bind(
                "reference_summary", WDL.Value.File(str(reference_summary))
                if use_reference else WDL.Value.Null()
            ).bind("min_mean_log2_cpm1", WDL.Value.Float(0.01)).bind(
                "residual_cutoff", WDL.Value.Null()
            )
            filter_command = render_after_localization(
                tasks["FilterCellTypeBeds"], filter_env, directory
            ).replace(
                "/opt/prepare_qtl/scripts/cell_type_specific_expression/filter_cell_type_beds.R",
                "pipeline/filter_cell_type_beds.R",
            )
            result = subprocess.run(["bash", "-c", filter_command], cwd=work,
                                    env=process_env, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            for relative_path in (
                "outputs/filtered_beds.txt", "outputs/filtered_inventory.tsv",
                "outputs/negative_summary.tsv.gz", "outputs/gene_comparison.tsv.gz",
                "outputs/filter_metrics.tsv",
            ):
                self.assertTrue((work / relative_path).exists(), relative_path)
            output_paths = (work / "outputs/filtered_beds.txt").read_text().splitlines()
            self.assertEqual(len(output_paths), 2)
            for output_path in output_paths:
                with gzip.open(work / output_path, "rt") as handle, gzip.open(bed, "rt") as original:
                    self.assertEqual(handle.read(), original.read())
            self.assertFalse((work / "unexpected_side_effect").exists())


if __name__ == "__main__":
    unittest.main()

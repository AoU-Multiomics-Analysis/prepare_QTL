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

    def test_filter_task_serializes_config_only_inside_task(self):
        task_path = ROOT / "workflows/cell_type_specific_expression/tasks/reference_filter.wdl"
        self.assertTrue(task_path.exists())
        document = WDL.load(str(task_path))
        tasks = {task.name: task for task in document.tasks}
        self.assertEqual(set(tasks), {"PrepareHaemopedia", "FilterCellTypeBeds"})
        filter_task = tasks["FilterCellTypeBeds"]
        self.assertEqual(str({decl.name: decl.type for decl in filter_task.inputs}["reference_summary"]),
                         "File?")
        source = task_path.read_text()
        command = str(filter_task.command)
        self.assertIn("File config_json = write_json(object", source)
        self.assertIn("filter_cell_type_beds.R", command)
        self.assertIn("status=failed", command)
        self.assertIn("completion_time", command)

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is required")
    def test_rendered_tasks_use_localized_paths_and_return_filter_outputs(self):
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
            counts = work / "counts with apostrophe's.tsv"
            populations = [
                "NveB", "MemB", "CD4T", "CD8T", "NK", "Mono",
                "MonoNonClassical", "Neut", "Eo", "myDC", "myDC123", "pDC",
            ]
            counts.write_text("gene_id\t" + "\t".join(f"{x}.1" for x in populations) + "\n" +
                              "ENSG000001\t" + "\t".join(["10"] * 12) + "\n" +
                              "ENSG000002\t" + "\t".join(["20"] * 12) + "\n" +
                              "ENSG000003\t" + "\t".join(["40"] * 12) + "\n")
            prepare_env = WDL.Env.Bindings().bind("counts", WDL.Value.File(str(counts)))
            stdlib = TaskStdLib("1.0", write_dir=directory)
            counts_path_decl = next(decl for decl in tasks["PrepareHaemopedia"].postinputs
                                    if decl.name == "counts_path_file")
            prepare_env = prepare_env.bind(
                "counts_path_file", counts_path_decl.expr.eval(prepare_env, stdlib)
            )
            prepare_command = tasks["PrepareHaemopedia"].command.eval(
                prepare_env, stdlib
            ).value.replace(
                "/opt/prepare_qtl/scripts/cell_type_specific_expression/prepare_haemopedia.R",
                "pipeline/prepare_haemopedia.R",
            )
            process_env = dict(os.environ, CELL_TYPE_SPECIFIC_EXPRESSION_ROOT=str(script_root))
            result = subprocess.run(["bash", "-c", prepare_command], cwd=work,
                                    env=process_env, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            reference_summary = work / "outputs/reference_summary.tsv.gz"
            self.assertTrue(reference_summary.exists())

            bed = work / "b_cells.bed.gz"
            with gzip.open(bed, "wt") as handle:
                handle.write("#chr\tstart\tend\tgene_id\ts1\ts2\n")
                handle.write("chr1\t0\t1\tENSG000001\t1\t2\n")
                handle.write("chr1\t1\t2\tENSG000002\t2\t3\n")
                handle.write("chr1\t2\t3\tENSG000003\t4\t5\n")
            inventory = work / "inventory.tsv"
            bed_sha256 = hashlib.sha256(bed.read_bytes()).hexdigest()
            inventory.write_text(
                "logical_name\tpath\tsha256\tn_genes\tn_samples\tscale\tcell_group\tslug\n"
                f"b_cells\t{bed.name}\t{bed_sha256}\t3\t2\tcpm\tB cells\tb_cells\n"
            )
            filter_env = WDL.Env.Bindings().bind(
                "cell_type_bed_inventory", WDL.Value.File(str(inventory))
            ).bind(
                "cell_type_beds", WDL.Value.Array(WDL.Type.File(), [WDL.Value.File(str(bed))])
            ).bind(
                "reference_summary", WDL.Value.File(str(reference_summary))
            ).bind("min_mean_log2_cpm1", WDL.Value.Float(0.01)).bind(
                "residual_cutoff", WDL.Value.Null()
            )
            stdlib = TaskStdLib("1.0", write_dir=directory)
            config_decl = next(decl for decl in tasks["FilterCellTypeBeds"].postinputs
                               if decl.name == "config_json")
            config_value = config_decl.expr.eval(filter_env, stdlib)
            filter_env = filter_env.bind("config_json", config_value)
            filter_command = tasks["FilterCellTypeBeds"].command.eval(filter_env, stdlib).value.replace(
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


if __name__ == "__main__":
    unittest.main()

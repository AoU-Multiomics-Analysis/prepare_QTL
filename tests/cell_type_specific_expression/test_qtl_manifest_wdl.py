"""Exercise the manifest task with upstream cloud File values, without cloud IO."""
import csv
import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import WDL

from test_reference_filter_wdl import CloudGeneratedFileStdLib, render_after_localization

ROOT = Path(__file__).resolve().parents[2]
SOURCES = {
    "int_beds": "PrepareCellTypeEqtl.IntBedFile",
    "scaled_beds": "PrepareCellTypeEqtl.ScaledBedFile",
    "int_phenotype_pcs": "PrepareCellTypeEqtl.IntPhenotypePCsOut",
    "int_phenotype_pcs_all": "PrepareCellTypeEqtl.IntPhenotypePCsAllOut",
    "scaled_phenotype_pcs": "PrepareCellTypeEqtl.ScaledPhenotypePCsOut",
    "scaled_phenotype_pcs_all": "PrepareCellTypeEqtl.ScaledPhenotypePCsAllOut",
    "int_merged_covariates": "int_merged_covariate",
    "scaled_merged_covariates": "scaled_merged_covariate",
    "int_connectivity_outliers": "PrepareCellTypeEqtl.IntConnectivityOutliers",
    "scaled_connectivity_outliers": "PrepareCellTypeEqtl.ScaledConnectivityOutliers",
}
COLUMNS = dict(zip(SOURCES, (
    "int_bed", "scaled_bed", "int_phenotype_pcs", "int_phenotype_pcs_all",
    "scaled_phenotype_pcs", "scaled_phenotype_pcs_all", "int_merged_covariates",
    "scaled_merged_covariates", "int_connectivity_outliers", "scaled_connectivity_outliers",
)))


class QtlManifestWdlTest(unittest.TestCase):
    def setUp(self):
        doc = WDL.load(str(ROOT / "workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl"))
        self.call = next(n for n in doc.workflow.body
                         if isinstance(n, WDL.Tree.Call) and n.name == "BuildQtlManifest")
        self.task = self.call.callee

    def task_inputs(self, directory=None):
        source_env = WDL.Env.Bindings()
        expected = {}
        for name, source in SOURCES.items():
            paths = [f"gs://test-bucket/submission/shard-{i}/{name}.tsv" for i in (1, 0)]
            expected[COLUMNS[name]] = paths
            source_env = source_env.bind(
                source, WDL.Value.Array(WDL.Type.File(), [WDL.Value.File(p) for p in paths])
            )
        source_env = source_env.bind(
            "PrepareScatterInputs.cell_types", WDL.Value.Array(
                WDL.Type.String(), [WDL.Value.String(s) for s in ("Monocytes", "CD4 T cells")])
        ).bind(
            "PrepareScatterInputs.cell_type_slugs", WDL.Value.Array(
                WDL.Type.String(), [WDL.Value.String(s) for s in ("monocytes", "cd4_t_cells")])
        )
        metadata = {
            "source_beds": ["gs://test-bucket/raw/monocytes.bed.gz",
                            "gs://test-bucket/raw/cd4_t_cells.bed.gz"],
            "filtered_beds": ["gs://test-bucket/filtered/monocytes.filtered.bed.gz",
                              "gs://test-bucket/filtered/cd4_t_cells.filtered.bed.gz"],
            "negative_summary": "gs://test-bucket/reports/negative_summary.tsv.gz",
            "gene_comparison": "gs://test-bucket/reports/gene_comparison.tsv.gz",
            "filter_metrics": "gs://test-bucket/reports/filter_metrics.tsv",
        }
        for name in ("source_beds", "filtered_beds"):
            expected["source_cpm_bed" if name == "source_beds" else "filtered_cpm_bed"] = metadata[name]
            source_env = source_env.bind(
                f"CellTypeDeconvolution.{'cell_type_beds' if name == 'source_beds' else 'filtered_cell_type_beds'}",
                WDL.Value.Array(WDL.Type.File(), [WDL.Value.File(path) for path in metadata[name]])
            )
        for name, output_name in (("negative_summary", "negative_expression_summary"),
                                  ("gene_comparison", "reference_gene_comparison"),
                                  ("filter_metrics", "reference_filter_metrics")):
            expected[output_name] = [metadata[name], metadata[name]]
            source_env = source_env.bind(
                f"CellTypeDeconvolution.{output_name}", WDL.Value.File(metadata[name]))
        inventory_dir = Path(directory or "/tmp")
        for filtered, expression in (
            (False, "CellTypeDeconvolution.cell_type_bed_inventory"),
            (True, "CellTypeDeconvolution.filtered_cell_type_bed_inventory"),
        ):
            path = inventory_dir / ("filtered_inventory.tsv" if filtered else "source_inventory.tsv")
            if directory:
                suffix = ".filtered" if filtered else ""
                path.write_text("path\tslug\nmonocytes%s.bed.gz\tmonocytes\ncd4_t_cells%s.bed.gz\tcd4_t_cells\n" %
                                (suffix, suffix))
            source_env = source_env.bind(expression, WDL.Value.File(str(path)))
        task_env = WDL.Env.Bindings()
        for decl in self.task.inputs:
            if decl.name in ("docker_image", "cpu", "memory", "disk_gb",
                             "preemptible_attempts", "max_retries"):
                continue
            value = self.call.inputs[decl.name].eval(source_env, WDL.StdLib.Base("1.0"))
            value = value.coerce(decl.type)
            # Path metadata must not be File-typed: localization would replace cloud URLs.
            if decl.name not in ("source_bed_inventory", "filtered_bed_inventory"):
                if isinstance(value.type, WDL.Type.Array):
                    self.assertIsInstance(value.type.item_type, WDL.Type.String)
                else:
                    self.assertIsInstance(value.type, WDL.Type.String)
            task_env = task_env.bind(decl.name, value)
        return task_env, expected

    def test_task_boundary_preserves_upstream_cloud_paths_as_metadata(self):
        env, expected = self.task_inputs()
        for name, column in COLUMNS.items():
            self.assertEqual(env[name].json, expected[column])

    def test_generated_metadata_files_remain_localizable_and_keep_cloud_urls(self):
        with tempfile.TemporaryDirectory() as directory:
            env, _ = self.task_inputs(directory)
            stdlib = CloudGeneratedFileStdLib(directory)
            for part in self.task.command.parts:
                if not isinstance(part, WDL.Expr.Placeholder) or "write_json" not in str(part.expr):
                    continue
                with self.subTest(expression=str(part.expr)):
                    generated = part.expr.eval(env, stdlib)
                    self.assertIsInstance(generated, WDL.Value.File)
                    localized = WDL.Value.rewrite_paths(generated,
                        lambda value: stdlib.generated_paths[value.value])
                    expected = part.expr.arguments[0].eval(env, stdlib).json
                    self.assertEqual(json.loads(Path(localized.value).read_text()), expected)

    def test_inventory_arguments_are_localized_before_shell_execution(self):
        # Serializing File paths in task declarations loses localization inside
        # the generated file. Exercise that boundary, not just local input paths.
        with tempfile.TemporaryDirectory(prefix="qtl-manifest-localization-") as directory:
            env, _ = self.task_inputs(directory)
            expected = {}
            for name in ("source_bed_inventory", "filtered_bed_inventory"):
                path = Path(directory) / f"{name} donor's $literal `echo unsafe` $(touch unexpected_side_effect).tsv"
                path.write_text(Path(env[name].value).read_text())
                env = env.bind(name, WDL.Value.File(str(path)))
                expected[name] = str(path)
            command = render_after_localization(self.task, env, directory)
            capture = '''Rscript() {
              printf '%s\\0' "$@" > captured_args
              printf 'entity:cell_type_id\\nmonocytes\\ncd4_t_cells\\n' > outputs/cell_type_qtl_manifest.tsv
            }
            '''
            result = subprocess.run(["bash", "-c", capture + command], cwd=directory,
                                    text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            arguments = (Path(directory) / "captured_args").read_bytes().decode().split(chr(0))[:-1]
            for name, flag in (("source_bed_inventory", "--source-inventory"),
                               ("filtered_bed_inventory", "--filtered-inventory")):
                self.assertEqual(arguments[arguments.index(flag) + 1], expected[name])
            self.assertFalse((Path(directory) / "unexpected_side_effect").exists())

    def test_unreadable_or_unlocalized_inventory_fails_before_r(self):
        for name in ("source_bed_inventory", "filtered_bed_inventory"):
            for kind in ("missing", "cloud"):
                with self.subTest(input=name, kind=kind), tempfile.TemporaryDirectory() as directory:
                    env, _ = self.task_inputs(directory)
                    path = ("gs://test-bucket/not-localized.tsv" if kind == "cloud"
                            else str(Path(directory) / "missing.tsv"))
                    env = env.bind(name, WDL.Value.File(path))
                    command = render_after_localization(self.task, env, directory)
                    result = subprocess.run(["bash", "-c",
                        'Rscript() { touch r_started; return 99; }\n' + command],
                        cwd=directory, text=True, capture_output=True)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertFalse((Path(directory) / "r_started").exists())
                    diagnostic = "localization_error" if kind == "cloud" else "inventory_not_readable"
                    self.assertIn(diagnostic, result.stdout + result.stderr)
                    self.assertIn(name, result.stdout + result.stderr)

    @unittest.skipUnless(shutil.which("Rscript"), "Rscript is required to execute the task command")
    def test_rendered_task_writes_full_paths_and_ids_without_opening_cloud_files(self):
        # GitHub's WDL-only validation host need not have the R image's packages.
        dependencies = subprocess.run(
            ["Rscript", "-e",
             'quit(status=as.integer(!all(vapply(c("optparse","jsonlite","purrr","tibble","readr"), '
             'requireNamespace, logical(1), quietly=TRUE))))'],
            text=True, capture_output=True,
        )
        if dependencies.returncode != 0:
            self.skipTest("R manifest packages are unavailable; the container smoke runs the task")
        with tempfile.TemporaryDirectory(prefix="qtl-manifest-task-") as directory:
            env, expected = self.task_inputs(directory)
            quoted_path = "gs://test-bucket/path with spaces/it's-$(touch unexpected_side_effect).tsv"
            expected["int_bed"][0] = quoted_path
            env = env.bind("int_beds", WDL.Value.Array(
                WDL.Type.String(), [WDL.Value.String(p) for p in expected["int_bed"]]))
            command = render_after_localization(self.task, env, directory)
            script = ROOT / "scripts/cell_type_specific_expression/downstream/build_qtl_manifest.R"
            # Rscript can encode spaces in --file before the script sees it.
            # A relative test-only link also matches the container's space-free path.
            (Path(directory) / "pipeline").symlink_to(script.parent, target_is_directory=True)
            command = command.replace("/cell_type_specific_expression/downstream/", "/cell_type_specific_expression/").replace(
                "/opt/prepare_qtl/scripts/cell_type_specific_expression/build_qtl_manifest.R",
                "pipeline/build_qtl_manifest.R",
            )
            process_env = dict(os.environ, CELL_TYPE_SPECIFIC_EXPRESSION_ROOT=str(script.parent.parent))
            result = subprocess.run(["bash", "-c", command], cwd=directory,
                                    env=process_env, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertFalse((Path(directory) / "unexpected_side_effect").exists())
            with open(Path(directory) / "outputs/cell_type_qtl_manifest.tsv") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                self.assertEqual(reader.fieldnames[0], "entity:cell_type_id")
                rows = list(reader)
            self.assertEqual([r["entity:cell_type_id"] for r in rows], ["monocytes", "cd4_t_cells"])
            for column, paths in expected.items():
                self.assertEqual([r[column] for r in rows], paths)
            self.assertIn("validated_cell_count:2", result.stdout)


if __name__ == "__main__":
    unittest.main()

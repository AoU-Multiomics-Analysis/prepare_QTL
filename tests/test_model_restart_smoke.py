"""Contract tests for the reusable cell-type model-restart smoke runner."""
import importlib.util
import gzip
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
RESTART_RUNNER = ROOT / "tests/cell_type_specific_expression/smoke/run_model_restart.py"
PINNED_RUNNER = ROOT / "tests/cell_type_specific_expression/smoke/run_pinned_images.py"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ModelRestartSmokeTest(unittest.TestCase):
    def test_focused_restart_checks_skips_and_beds_without_qtl_outputs(self):
        module = load_module('focused_restart', RESTART_RUNNER)
        prefix = module.PREFIX
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            beds = []
            for name in ('b_cells', 'cd4_t_cells'):
                bed = root / (name + '.bed.gz')
                with gzip.open(bed, 'wt') as stream:
                    stream.write('gene\tS1\ng1\t1\n')
                beds.append(str(bed))
            parameters = root / 'parameters.json'
            parameters.write_text(json.dumps({'proportion_mode': 'precomputed_model'}))
            images = {'fit': 'fit-digest', 'export': 'export-digest', 'qtl': 'qtl-digest'}
            baseline = {prefix + 'tca_model': 'model.rds', prefix + 'stage_images': images,
                        prefix + 'cell_type_beds': beds, prefix + 'filtered_cell_type_beds': beds}
            result = {'CellTypeDeconvolution.' + name: None for name in (
                'estimated_proportions', 'tca_model_unfiltered', 'fit_tca_log',
                'proportions_lm22', 'proportions_combined', 'gene_type_filter_log')}
            result.update({'CellTypeDeconvolution.cell_type_beds': beds,
                           'CellTypeDeconvolution.filtered_cell_type_beds': beds,
                           'CellTypeDeconvolution.stage_images': {k: v for k, v in images.items() if k != 'qtl'},
                           'CellTypeDeconvolution.effective_parameters_file': str(parameters)})
            baseline_path, inputs_path = root / 'baseline.json', root / 'inputs.json'
            baseline_path.write_text(json.dumps(baseline))
            inputs_path.write_text(json.dumps({prefix + 'lm22': 'reference.tsv', prefix + 'expression': 'bulk.bed'}))
            output = root / 'restart'
            output.mkdir()
            (output / 'outputs.json').write_text(json.dumps(result))
            commands = []
            with mock.patch.object(module.subprocess, 'run', side_effect=lambda command, **kw: commands.append(command)):
                argv = ['--baseline-outputs', str(baseline_path), '--baseline-inputs', str(inputs_path),
                        '--restart-inputs', str(root / 'restart.json'), '--output-directory', str(output),
                        '--deconvolution-only']
                module.main(argv)
                self.assertEqual(commands[0][2], 'workflows/cell_type_specific_expression/deconvolution.wdl')
                generated = json.loads((root / 'restart.json').read_text())
                self.assertEqual(generated['CellTypeDeconvolution.precomputed_tca_model'], 'model.rds')
                self.assertEqual(generated['CellTypeDeconvolution.covariates'], 'reference.tsv')
                result['CellTypeDeconvolution.fit_tca_log'] = 'unexpected-fit.log'
                (output / 'outputs.json').write_text(json.dumps(result))
                with self.assertRaisesRegex(AssertionError, 'fit_tca_log should be skipped'):
                    module.main(argv)

    def test_restart_runner_reads_flat_miniwdl_outputs(self):
        module = load_module("model_restart_smoke", RESTART_RUNNER)
        expected = {"PrepareCellTypeEqtlWorkflow.tca_model": "model.rds"}
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "outputs.json"
            output.write_text(json.dumps(expected))
            self.assertEqual(module.read_outputs(output), expected)

    def test_restart_runner_exposes_path_overrides(self):
        result = subprocess.run(
            [sys.executable, str(RESTART_RUNNER), "--help"],
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        for option in (
            "--baseline-outputs",
            "--baseline-inputs",
            "--restart-inputs",
            "--output-directory",
        ):
            self.assertIn(option, result.stdout)

    def test_pinned_runner_passes_its_precomputed_paths_to_restart(self):
        # WDL is used only by main(); this test exercises the pure command builder.
        with mock.patch.dict(sys.modules, {"WDL": types.ModuleType("WDL")}):
            module = load_module("pinned_restart_smoke", PINNED_RUNNER)
        self.assertTrue(hasattr(module, "restart_command"))
        command = module.restart_command(
            Path("ci-runs/pinned-precomputed/outputs.json"),
            Path("ci-runs/pinned-precomputed.inputs.json"),
            Path("ci-runs/pinned-model-restart.inputs.json"),
            Path("ci-runs/pinned-model-restart"),
        )
        self.assertEqual(
            command,
            [
                sys.executable,
                "tests/cell_type_specific_expression/smoke/run_model_restart.py",
                "--baseline-outputs",
                "ci-runs/pinned-precomputed/outputs.json",
                "--baseline-inputs",
                "ci-runs/pinned-precomputed.inputs.json",
                "--restart-inputs",
                "ci-runs/pinned-model-restart.inputs.json",
                "--output-directory",
                "ci-runs/pinned-model-restart",
            ],
        )


if __name__ == "__main__":
    unittest.main()

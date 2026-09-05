"""Contract tests for the reusable cell-type model-restart smoke runner."""
import importlib.util
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

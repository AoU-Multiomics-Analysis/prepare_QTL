"""Ensure smoke inputs cannot silently replace immutable workflow defaults."""
import importlib.util
from pathlib import Path
import unittest
import tempfile
import json

ROOT = Path(__file__).resolve().parents[1]


class PinnedSmokeTest(unittest.TestCase):
    def test_reads_flat_miniwdl_output_file(self):
        path = ROOT / 'tests/cell_type_specific_expression/smoke/run_pinned_images.py'
        spec = importlib.util.spec_from_file_location('pinned_smoke', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertTrue(hasattr(module, 'read_outputs'))
        expected = {'PrepareCellTypeEqtlWorkflow.stage_images': {'fit': 'pinned-image'}}
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'outputs.json'
            output.write_text(json.dumps(expected))
            self.assertEqual(module.read_outputs(output), expected)

    def test_fixture_image_overrides_are_replaced_but_analysis_inputs_are_kept(self):
        path = ROOT / 'tests/cell_type_specific_expression/smoke/run_pinned_images.py'
        self.assertTrue(path.exists(), 'Pinned smoke runner is not implemented')
        spec = importlib.util.spec_from_file_location('pinned_smoke', path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        prefix = 'PrepareCellTypeEqtlWorkflow.'
        fixture = {prefix + 'expression': 'expression.bed', prefix + 'fit_docker_image': 'local:test'}
        pins = {'fit_docker_image': 'ghcr.io/example/image@sha256:' + 'a' * 64}
        result = module.pinned_inputs(fixture, pins)
        self.assertEqual(result, {prefix + 'expression': 'expression.bed', prefix + 'fit_docker_image': pins['fit_docker_image']})
        self.assertEqual(fixture[prefix + 'fit_docker_image'], 'local:test')


if __name__ == '__main__':
    unittest.main()

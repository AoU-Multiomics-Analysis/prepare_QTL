"""Check task image inputs and metadata across the nested cell workflows."""
import os
from pathlib import Path
import sys
import unittest
import WDL

POLICY_ROOT = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('RELEASE_SOURCE_ROOT', POLICY_ROOT)).resolve()
sys.path.insert(0, str(POLICY_ROOT / 'ci'))
from wdl_stage_routing import validate_routing

class TaskImagesTest(unittest.TestCase):
    def test_nested_defaults_and_metadata(self):
        base = ROOT / 'workflows/cell_type_specific_expression'
        parent = WDL.load(str(base / 'prepare_cell_type_eQTL.wdl')).workflow
        child = WDL.load(str(base / 'deconvolution.wdl')).workflow
        def pins(w):
            return {d.name: d.expr.literal.value for d in w.inputs if d.name.endswith('_image')}
        for name, image in pins(child).items():
            self.assertEqual(image, pins(parent)[name])
            self.assertRegex(image, r'^ghcr\.io/[^\s@]+@sha256:[0-9a-f]{64}$')
        for workflow in (parent, child):
            env = WDL.Env.Bindings()
            expected = {}
            for name in pins(workflow):
                value = name + '-independent-test-image'
                env = env.bind(name, WDL.Value.String(value))
                expected[name.removesuffix('_image')] = value
            output = next(d for d in workflow.outputs if d.name == 'task_images')
            self.assertEqual(output.expr.eval(env, WDL.StdLib.Base('1.0')).json, expected)
        self.assertEqual(validate_routing(ROOT, POLICY_ROOT), [])

if __name__ == '__main__':
    unittest.main()

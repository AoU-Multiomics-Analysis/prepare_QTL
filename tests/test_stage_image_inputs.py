"""Catch image cross-wiring and mutable defaults at actual WDL call boundaries."""
from pathlib import Path
import os
import sys
import unittest
import WDL

POLICY_ROOT = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get('RELEASE_SOURCE_ROOT', POLICY_ROOT)).resolve()
sys.path.insert(0, str(POLICY_ROOT / 'ci'))
from wdl_stage_routing import validate_routing
BASE = ROOT / 'workflows/cell_type_specific_expression'
STAGES = ('estimation', 'fit', 'export', 'downstream')


class StageImagesTest(unittest.TestCase):
    def setUp(self):
        self.child = WDL.load(str(BASE / 'deconvolution.wdl')).workflow
        self.parent = WDL.load(str(BASE / 'prepare_cell_type_eQTL.wdl')).workflow

    def environment(self, downstream='downstream-image'):
        env = WDL.Env.Bindings()
        for stage in (*STAGES, 'qtl'):
            env = env.bind(stage + '_docker_image', WDL.Value.String(
                downstream if stage == 'downstream' else stage + '-image'))
        return env

    def test_defaults_are_immutable_and_entry_points_agree(self):
        child = {d.name: d for d in self.child.inputs}
        parent = {d.name: d for d in self.parent.inputs}
        for stage in (*STAGES, 'qtl'):
            name = stage + '_docker_image'
            self.assertIn(name, parent)
            value = parent[name].expr.eval(WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value
            self.assertRegex(value, r'^ghcr\.io/[^\s@]+@sha256:[0-9a-f]{64}$')
            if stage != 'qtl':
                self.assertIn(name, child)
                self.assertEqual(value, child[name].expr.eval(WDL.Env.Bindings(), WDL.StdLib.Base('1.0')).value)

    def test_each_call_receives_only_its_stage_image(self):
        # The validator follows every call, including future tasks, and derives
        # stage identity from trusted pin inputs even when image digests match.
        self.assertEqual(validate_routing(ROOT, POLICY_ROOT), [])

    def test_output_records_selected_stage_images(self):
        for workflow, groups in [(self.child, STAGES), (self.parent, (*STAGES, 'qtl'))]:
            declaration = next(d for d in workflow.outputs if d.name == 'stage_images')
            value = declaration.expr.eval(self.environment('updated-downstream'), WDL.StdLib.Base('1.0')).json
            self.assertEqual(value, {stage: 'updated-downstream' if stage == 'downstream'
                                    else stage + '-image' for stage in groups})


if __name__ == '__main__':
    unittest.main()

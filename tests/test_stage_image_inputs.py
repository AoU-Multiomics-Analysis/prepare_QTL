"""Catch image cross-wiring and mutable defaults at actual WDL call boundaries."""
from pathlib import Path
import unittest
import WDL

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / 'workflows/cell_type_specific_expression'
STAGES = ('estimation', 'fit', 'export', 'downstream')


def calls(nodes):
    for node in nodes:
        if isinstance(node, WDL.Tree.Call):
            yield node
        elif isinstance(node, (WDL.Tree.Scatter, WDL.Tree.Conditional)):
            yield from calls(node.body)


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
        groups = {
            'estimation': ['FilterExpressionGenes', 'ValidateProportionMode', 'PrepareHspeBatches',
                           'RunHspeBatch', 'MergeHspeBatches', 'ProcessProportions'],
            'fit': ['FitTca', 'CleanTcaModel'],
            'export': ['ExportTcaBeds'],
            'downstream': ['SummarizeCellTypeBeds', 'PrepareHaemopedia', 'FilterCellTypeBeds', 'BuildManifest'],
        }
        expected = {name: group for group, names in groups.items() for name in names}
        actual_calls = list(calls(self.child.body))
        self.assertEqual(set(expected), {c.name for c in actual_calls})
        for downstream in ('downstream-image', 'replacement-image'):
            for call in actual_calls:
                group = expected[call.name]
                want = downstream if group == 'downstream' else group + '-image'
                self.assertEqual(call.inputs['docker_image'].eval(
                    self.environment(downstream), WDL.StdLib.Base('1.0')).value, want, call.name)

    def test_parent_routes_images_to_child_and_qtl(self):
        for call in calls(self.parent.body):
            if call.name == 'CellTypeDeconvolution':
                for stage in STAGES:
                    name = stage + '_docker_image'
                    self.assertIn(name, call.inputs)
                    self.assertEqual(call.inputs[name].eval(self.environment(), WDL.StdLib.Base('1.0')).value,
                                     stage + '-image')
            else:
                key = 'DockerImage' if call.name == 'PrepareCellTypeEqtl' else 'docker_image'
                want = 'qtl-image' if key == 'DockerImage' else 'downstream-image'
                self.assertEqual(call.inputs[key].eval(self.environment(), WDL.StdLib.Base('1.0')).value, want)

    def test_output_records_selected_stage_images(self):
        for workflow, groups in [(self.child, STAGES), (self.parent, (*STAGES, 'qtl'))]:
            declaration = next(d for d in workflow.outputs if d.name == 'stage_images')
            value = declaration.expr.eval(self.environment('updated-downstream'), WDL.StdLib.Base('1.0')).json
            self.assertEqual(value, {stage: 'updated-downstream' if stage == 'downstream'
                                    else stage + '-image' for stage in groups})


if __name__ == '__main__':
    unittest.main()

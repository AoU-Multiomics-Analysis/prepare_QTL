"""Evaluate restart branch guards with model/proportion input combinations."""
import unittest
import shlex
from pathlib import Path
import WDL

ROOT = Path(__file__).resolve().parents[2]


class ModelRestartTest(unittest.TestCase):
    def test_task_commands_enable_restart_only_when_requested(self):
        doc = WDL.load(str(ROOT / 'workflows/cell_type_specific_expression/tasks/tca.wdl'))
        for task in doc.tasks:
            if task.name not in ('CleanTcaModel', 'ExportTcaBeds'):
                continue
            for enabled in (False, True):
                env = WDL.Env.Bindings().bind('parallel_argument', WDL.Value.String(''))
                for decl in task.inputs:
                    if decl.name == 'reuse_model':
                        value = WDL.Value.Boolean(enabled)
                    elif decl.expr is not None:
                        value = decl.expr.eval(env, WDL.StdLib.Base('1.0'))
                    elif decl.type.optional:
                        value = WDL.Value.Null()
                    elif isinstance(decl.type, WDL.Type.File):
                        value = WDL.Value.File('/localized/input')
                    else:
                        value = WDL.Value.String('image')
                    env = env.bind(decl.name, value)
                command = task.command.eval(env, WDL.StdLib.Base('1.0')).value
                with self.subTest(task=task.name, enabled=enabled):
                    self.assertEqual('--reuse-model' in shlex.split(command), enabled)

    def test_restart_skips_estimation_and_fit(self):
        doc = WDL.load(str(ROOT / 'workflows/cell_type_specific_expression/deconvolution.wdl'))
        inputs = {n.name: n for n in doc.workflow.inputs}
        self.assertIn('precomputed_tca_model', inputs)
        self.assertEqual(str(inputs['precomputed_tca_model'].type), 'File?')
        for model in (False, True):
            for proportions in (False, True):
                env = WDL.Env.Bindings().bind('haemopedia_counts', WDL.Value.Null())
                for name, present in [('precomputed_tca_model', model),
                                      ('precomputed_proportions', proportions)]:
                    env = env.bind(name, WDL.Value.File('/input.rds') if present else WDL.Value.Null())
                calls = set()
                for node in doc.workflow.body:
                    if isinstance(node, WDL.Tree.Call):
                        calls.add(node.name)
                    elif isinstance(node, WDL.Tree.Conditional):
                        if node.expr.eval(env, WDL.StdLib.Base('1.0')).value:
                            calls.update(n.name for n in node.body if isinstance(n, WDL.Tree.Call))
                with self.subTest(model=model, proportions=proportions):
                    self.assertEqual('FitTca' in calls, not model)
                    self.assertEqual('ProcessProportions' in calls, not model)
                    self.assertEqual('PrepareHspeBatches' in calls, not model and not proportions)
                    self.assertIn('CleanTcaModel', calls)
                    self.assertIn('ExportTcaBeds', calls)

    def test_parent_forwards_model(self):
        doc = WDL.load(str(ROOT / 'workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl'))
        call = next(n for n in doc.workflow.body if isinstance(n, WDL.Tree.Call))
        self.assertIn('precomputed_tca_model', call.inputs)


if __name__ == '__main__':
    unittest.main()

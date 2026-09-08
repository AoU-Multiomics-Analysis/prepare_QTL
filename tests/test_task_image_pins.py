"""Task-level release impact and pin-isolation regression tests."""
from pathlib import Path
import sys
import copy
import tempfile
import shutil
import unittest
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'ci'))
from plan_image_updates import plan_changes, validate_registry
from propose_image_pins import propose
from wdl_stage_routing import validate_routing

class TaskPins(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = yaml.safe_load((ROOT / 'ci/image-stages.yml').read_text())
        cls.targets = yaml.safe_load((ROOT / 'ci/release-pins.yml').read_text())

    def test_hspe_edit_preserves_filter_pin(self):
        plan = plan_changes(self.config, ['scripts/cell_type_specific_expression/estimation/run_hspe_batch.R'])
        self.assertEqual(plan['stages'], ['hspe__run_hspe_batch'])
        files = {p['path']: (ROOT / p['path']).read_text()
                 for targets in self.targets['stages'].values() for p in targets}
        new_image = self.config['images']['cell_type']['repository'] + '@sha256:' + 'a' * 64
        updated = propose(self.config, self.targets, plan, files, {'cell_type': new_image})
        filter_path = 'workflows/cell_type_specific_expression/tasks/expression.wdl'
        self.assertNotIn(filter_path, updated)
        self.assertIn('workflows/cell_type_specific_expression/tasks/hspe.wdl', updated)

    def test_environment_updates_all_consumers_only(self):
        plan = plan_changes(self.config, ['envs/CellTypeSpecificExpression/Dockerfile'])
        self.assertEqual(set(plan['stages']), {k for k, v in self.config['stages'].items() if v['image'] == 'cell_type'})

    def test_shared_helpers_update_all_users(self):
        plan = plan_changes(self.config, ['scripts/methylation/MethylationUtils.R'])
        expected = {k for k, v in self.config['stages'].items() if 'scripts/methylation/MethylationUtils.R' in v['sources']}
        self.assertGreater(len(expected), 1)
        self.assertEqual(set(plan['stages']), expected)

    def test_unknown_runtime_script_requires_mapping(self):
        path = 'scripts/cell_type_specific_expression/estimation/new_method.R'
        plan = plan_changes(self.config, [path])
        self.assertEqual(plan['unmapped'], [path])

    def test_hspe_helper_does_not_update_expression_filter(self):
        plan = plan_changes(self.config, ['scripts/cell_type_specific_expression/R/hspe_stage.R'])
        self.assertNotIn('expression__filter_expression_genes', plan['stages'])
        self.assertIn('hspe__run_hspe_batch', plan['stages'])

    def test_trans_ld_script_updates_only_its_task(self):
        plan = plan_changes(self.config, ['tools/trans_ld_regions/scripts/trans_ld_regions/ld.py'])
        self.assertEqual(plan['stages'], ['trans_ld_ancestry'])
        self.assertEqual(plan['builds'], ['trans_ld'])

    def test_missing_transitive_module_is_rejected(self):
        from task_dependencies import validate_task_dependencies
        config = copy.deepcopy(self.config)
        config['stages']['expression__filter_expression_genes']['sources'].remove(
            'scripts/cell_type_specific_expression/R/expression.R')
        errors = validate_task_dependencies(config, ROOT)
        self.assertTrue(any('missing transitive R modules' in e for e in errors), errors)

    def test_new_helper_call_requires_dependency_update(self):
        from task_dependencies import validate_task_dependencies
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for folder in ('scripts', 'tools', 'rust'):
                shutil.copytree(ROOT / folder, root / folder)
            helper = root / 'scripts/cell_type_specific_expression/R/expression.R'
            helper.write_text(helper.read_text() + '\nnew_alias <- function() read_cell_type_mapping(NULL)\n')
            errors = validate_task_dependencies(self.config, root)
            self.assertTrue(any('expression__filter_expression_genes' in e and
                                'reference_mapping.R' in e for e in errors), errors)

    def test_complete_registry_and_routing(self):
        self.assertEqual(validate_registry(self.config, ROOT), [])
        self.assertEqual(validate_routing(ROOT, ROOT), [])

if __name__ == '__main__':
    unittest.main()

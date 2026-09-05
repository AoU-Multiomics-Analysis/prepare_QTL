import importlib.util
import subprocess
import sys
from pathlib import Path
import unittest

import yaml

ROOT = Path(__file__).resolve().parents[1]


class ImagePlanTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location('image_plan', ROOT / 'ci/plan_image_updates.py')
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)
        cls.config = yaml.safe_load((ROOT / 'ci/image-stages.yml').read_text())

    def test_filter_edit_changes_only_downstream_pin(self):
        plan = self.module.plan_changes(self.config, ['scripts/cell_type_specific_expression/filter_cell_type_beds.R'])
        self.assertEqual(plan['stages'], ['cell_downstream'])
        self.assertEqual(plan['builds'], ['cell_type', 'standard'])
        self.assertEqual(plan['unmapped'], [])

    def test_shared_cell_io_changes_all_cell_stages(self):
        plan = self.module.plan_changes(self.config, ['scripts/cell_type_specific_expression/R/io.R'])
        self.assertEqual(plan['stages'], ['cell_downstream', 'cell_estimation', 'cell_export', 'cell_fit'])

    def test_standard_expression_and_rust_changes_are_separate(self):
        plan = self.module.plan_changes(self.config, ['scripts/expression/PrepareExpression.R', 'rust/methylation_filter/src/main.rs'])
        self.assertEqual(plan['stages'], ['expression', 'methylation_rust'])
        self.assertEqual(plan['builds'], ['methylation_rust', 'standard'])

    def test_environment_changes_update_all_consumers(self):
        plan = self.module.plan_changes(self.config, ['envs/PhenotypePCs/Dockerfile'])
        self.assertEqual(plan['stages'], ['expression', 'methylation', 'proteomics', 'splicing'])

    def test_wdl_and_tests_do_not_request_builds_or_pins(self):
        paths = ['workflows/genotype/prepare_VCF.wdl', 'tests/example_test.py']
        plan = self.module.plan_changes(self.config, paths)
        self.assertEqual(plan['builds'], [])
        self.assertEqual(plan['stages'], [])
        self.assertEqual(plan['wdl_checks'], [paths[0]])
        self.assertEqual(plan['test_changes'], [paths[1]])

    def test_new_source_requires_explicit_mapping(self):
        plan = self.module.plan_changes(self.config, ['scripts/cell_type_specific_expression/new_method.R'])
        self.assertEqual(plan['unmapped'], ['scripts/cell_type_specific_expression/new_method.R'])

    def test_all_tracked_sources_and_workflows_are_classified(self):
        errors = self.module.validate_registry(self.config, ROOT)
        self.assertEqual(errors, [])

    def test_unknown_workflow_requires_review(self):
        plan = self.module.plan_changes(self.config, ['workflows/new_assay/run.wdl'])
        self.assertEqual(plan['unmapped'], ['workflows/new_assay/run.wdl'])

    def test_cli_reports_selected_group(self):
        result = subprocess.run([sys.executable, str(ROOT / 'ci/plan_image_updates.py'),
                                 '--changed', 'scripts/cell_type_specific_expression/filter_cell_type_beds.R'],
                                capture_output=True, text=True, cwd=ROOT)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn('`cell_downstream`', result.stdout)
        self.assertNotIn('`cell_fit`', result.stdout)

    def test_cli_rejects_unmapped_source(self):
        result = subprocess.run([sys.executable, str(ROOT / 'ci/plan_image_updates.py'),
                                 '--changed', 'scripts/new_assay/tool.R'],
                                capture_output=True, text=True, cwd=ROOT)
        self.assertEqual(result.returncode, 1)
        self.assertIn('scripts/new_assay/tool.R', result.stdout)


if __name__ == '__main__':
    unittest.main()

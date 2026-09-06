import copy
import json
import shutil
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'ci'))
import test_release_images as runner
from plan_image_updates import plan_changes


class StageTestSelection(unittest.TestCase):
    def setUp(self):
        self.config = yaml.safe_load((ROOT / 'ci/image-stages.yml').read_text())

    def select(self, stages, paths):
        self.assertTrue(hasattr(runner, 'runtime_test_plan'), 'Stage-specific test planning is missing')
        return runner.runtime_test_plan(self.config, set(stages), paths)

    def test_fit_only_does_not_run_estimation_or_qtl(self):
        plan = self.select(['cell_fit'], ['scripts/cell_type_specific_expression/fit/fit_tca.R'])
        self.assertEqual(plan['stages'], ['cell_fit'])
        self.assertFalse(plan['integration'])
        self.assertIn('test-tca-stage.R', plan['cell_tests']['cell_fit'])
        self.assertIn('test-fit-cli.R', plan['cell_tests']['cell_fit'])

    def test_plot_change_does_not_fit_tca(self):
        plan = self.select(['cell_export', 'cell_downstream'], ['scripts/cell_type_specific_expression/R/qc.R'])
        self.assertFalse(plan['integration'])
        self.assertNotIn('cell_fit', plan['cell_tests'])

    def test_real_source_mapping_selects_only_hspe_tests(self):
        paths = ['scripts/cell_type_specific_expression/R/hspe_stage.R']
        changes = plan_changes(self.config, paths)
        plan = self.select(changes['stages'], paths)
        self.assertEqual(plan['stages'], ['cell_estimation'])
        self.assertIn('test-hspe-batches.R', plan['cell_tests']['cell_estimation'])
        self.assertFalse(plan['integration'])

    def test_manual_override_and_missing_evidence_keep_integration(self):
        self.assertTrue(hasattr(runner, 'runtime_test_plan'))
        self.assertTrue(runner.runtime_test_plan(self.config, {'cell_fit'},
            ['scripts/cell_type_specific_expression/fit/fit_tca.R'], True)['integration'])
        self.assertTrue(self.select(['cell_fit'], [])['integration'])

    def test_all_registered_stage_tests_exist_and_exemptions_are_exact_existing_files(self):
        self.assertTrue(hasattr(runner, 'runtime_test_plan'))
        plan = runner.runtime_test_plan(self.config, set(self.config['stages']), [], True)
        for files in plan['cell_tests'].values():
            for filename in files:
                self.assertTrue((ROOT / 'tests/cell_type_specific_expression/testthat' / filename).is_file())
        for path in self.config['stage_only_test_paths']:
            self.assertTrue((ROOT / path).is_file(), path)

    def test_bed_interface_shared_code_and_unknown_source_keep_integration(self):
        for path in ['scripts/cell_type_specific_expression/R/bed_outputs.R',
                     'scripts/cell_type_specific_expression/shared/new.R',
                     'scripts/cell_type_specific_expression/fit/new.R',
                     'envs/CellTypeSpecificExpression/environment.yml',
                     'workflows/cell_type_specific_expression/tasks/tca.wdl']:
            with self.subTest(path=path):
                self.assertTrue(self.select(['cell_export'], [path])['integration'])

    def test_non_cell_images_do_not_start_cell_pipeline(self):
        self.assertFalse(self.select(['rnaseqc'], ['envs/RNASeQCAggregation/Dockerfile'])['integration'])

    def test_missing_stage_gate_fails_closed(self):
        self.assertTrue(hasattr(runner, 'runtime_test_plan'))
        config = copy.deepcopy(self.config)
        config['stages']['cell_fit'].pop('runtime_tests', None)
        with self.assertRaisesRegex(ValueError, 'runtime tests'):
            runner.runtime_test_plan(config, {'cell_fit'}, [])

    def test_release_execution_uses_selected_image_without_launching_workflow(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'source'
            source.mkdir()
            for name in ('workflows', 'scripts', 'ci'):
                shutil.copytree(ROOT / name, source / name)
            runner.subprocess.run(['git', 'init', str(source)], check=True, capture_output=True)
            runner.subprocess.run(['git', '-C', str(source), 'add', '.'], check=True, capture_output=True)
            record = Path(directory) / 'release.json'
            record.write_text(json.dumps({'plan': {'stages': ['cell_fit'], 'wdl_checks': []},
                'pr': {'base': {'sha': 'a' * 40}, 'head': {'sha': 'b' * 40}}}))
            # Docker and runner processes are external boundaries. Keep plan
            # construction, WDL routing, pin parsing and artifact writing real.
            original_run = runner.subprocess.run
            original_output = runner.subprocess.check_output

            def output(command, **kwargs):
                if 'diff' in command:
                    return b'scripts/cell_type_specific_expression/fit/fit_tca.R\0'
                return original_output(command, **kwargs)

            def run(command, **kwargs):
                if command[0] == 'git':
                    return original_run(command, **kwargs)
                return runner.subprocess.CompletedProcess(command, 0)

            with mock.patch.object(sys, 'argv', ['runner', '--source', str(source), '--record', str(record)]), \
                 mock.patch.object(runner.subprocess, 'check_output', side_effect=output), \
                 mock.patch.object(runner.subprocess, 'run', side_effect=run) as execute:
                runner.main()
            plan = json.loads((source / 'ci-runs/runtime-test-plan.json').read_text())
            self.assertFalse(plan['integration'])
            commands = [call.args[0] for call in execute.call_args_list]
            stage_commands = [cmd for cmd in commands if 'tests/release/test_cell_stage.R' in cmd]
            self.assertEqual(len(stage_commands), 1)
            self.assertIn('test-fit-cli.R', stage_commands[0])
            self.assertNotIn('test-hspe-batches.R', stage_commands[0])
            self.assertTrue(any('@sha256:' in arg for arg in stage_commands[0]))
            self.assertFalse(any(any('run_pinned_images.py' in arg for arg in cmd) for cmd in commands))


if __name__ == '__main__':
    unittest.main()

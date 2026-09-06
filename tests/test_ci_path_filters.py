"""Check which expensive CI jobs a changed file selects (positive globs only)."""
from fnmatch import fnmatchcase
from pathlib import Path
import unittest
import os
import subprocess
import tempfile

import yaml

ROOT = Path(__file__).resolve().parents[1]


class CiPathTests(unittest.TestCase):
    def setUp(self):
        self.workflows = [yaml.load(p.read_text(), Loader=yaml.BaseLoader)
                          for p in (ROOT / '.github/workflows').glob('*.yml')]

    def heavy_jobs(self, path, event):
        jobs = set()
        for workflow in self.workflows:
            if event not in workflow['on']:
                continue
            trigger = workflow['on'][event] or {}
            patterns = trigger.get('paths')
            if patterns is not None and not any(fnmatchcase(path, p) for p in patterns):
                continue
            for name, job in workflow['jobs'].items():
                if any('docker/build-push-action@' in step.get('uses', '') or
                       'docker build ' in step.get('run', '')
                       for step in job.get('steps', [])):
                    jobs.add(name)
        return jobs

    def test_changed_files_select_only_relevant_builds(self):
        cases = {
            'workflows/cell_type_specific_expression/deconvolution.wdl': set(),
            'workflows/expression/rnaseqc2_aggregate_batched.wdl': set(),
            'docs/terra-file-paths.md': set(),
            'scripts/cell_type_specific_expression/run_hspe.R': set(),
            'scripts/common/MergeCovariates.R': set(),
            'scripts/expression/merge_rnaseqc.py': set(),
            'rust/methylation_merge/src/main.rs': set(),
            'envs/CellTypeSpecificExpression/environment.yml': set(),
            'envs/PhenotypePCs/Dockerfile': set(),
            'envs/MethylationRust/Dockerfile': set(),
            'envs/RNASeQCAggregation/environment.yml': set(),
            'tests/cell_type_specific_expression/fixtures/hspe-e2e.inputs.json': {'smoke'},
            'tests/test_prepare_expression_sample_list.R': {'smoke'},
            'tests/test_prepare_methylation.R': {'smoke'},
            'tests/rnaseqc2_aggregation/smoke_container.py': {'container'},
            'tests/cell_type_specific_expression/test_reference_filter_wdl.py': set(),
            '.dockerignore': set(),
        }
        for event in ('push', 'pull_request'):
            for path, expected in cases.items():
                with self.subTest(event=event, path=path):
                    self.assertEqual(self.heavy_jobs(path, event), expected if event == 'pull_request' else set())

    def test_manual_dispatch_runs_all_builds(self):
        self.assertEqual(self.heavy_jobs('', 'workflow_dispatch'),
                         {'build_and_push', 'build_cell_type_specific_expression',
                          'build_methylation_rust', 'smoke', 'container'})


class IntegrationSelectionTests(unittest.TestCase):
    def setUp(self):
        workflow = yaml.load(
            (ROOT / '.github/workflows/cell-type-specific-expression-ci.yml').read_text(),
            Loader=yaml.BaseLoader)
        self.steps = workflow['jobs']['smoke']['steps']

    def test_full_workflow_steps_share_gate_but_builds_and_r_tests_do_not(self):
        gated = {step['name'] for step in self.steps if 'if' in step
                 and step['if'] == "steps.integration.outputs.run == 'true'"}
        self.assertEqual(gated, {
            'Install MiniWDL 1.15.0',
            'Run the hspe end-to-end smoke workflow',
            'Assert the hspe deconvolution outputs', 'Assert the hspe QTL outputs',
            'Run the precomputed end-to-end smoke workflow',
            'Assert the precomputed deconvolution outputs', 'Assert the precomputed QTL outputs',
            'Restart from the cleaned model and compare outputs',
        })
        for step in self.steps:
            if 'docker build ' in step.get('run', '') or step.get('name', '').startswith(
                    ('Run the cell-type-specific R suite', 'Test pre-normalized', 'Test methylation')):
                self.assertNotIn('if', step)

    def test_real_selector(self):
        command = next(step['run'] for step in self.steps if step.get('id') == 'integration')
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)

            def git(*args):
                return subprocess.check_output(['git', '-C', directory, *args], text=True).strip()

            git('init', '-q')
            git('config', 'user.name', 'Test')
            git('config', 'user.email', 'test@example.invalid')

            def commit(path):
                target = repo / path
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text('fixture')
                git('add', path)
                git('commit', '-qm', 'fixture')
                return git('rev-parse', 'HEAD')

            base = commit('README')
            cases = [
                ('tests/test_prepare_expression_sample_list.R', 'false'),
                ('tests/cell_type_specific_expression/testthat/test-fit.R', 'false'),
                ('.github/workflows/cell-type-specific-expression-ci.yml', 'false'),
                ('tests/cell_type_specific_expression/smoke/assert_qtl_outputs.R', 'true'),
                ('tests/cell_type_specific_expression/fixtures/input.json', 'true'),
                ('tests/cell_type_specific_expression/generate_reference_fixture.R', 'true'),
            ]

            def select(before, head, event='pull_request'):
                output = repo / 'job-output'
                output.write_text('')
                subprocess.run(['bash', '-c', command], cwd=repo, check=True,
                               env={**os.environ, 'BASE': before, 'HEAD': head,
                                    'EVENT': event, 'GITHUB_OUTPUT': str(output)},
                               capture_output=True, text=True)
                return output.read_text().strip()

            for path, expected in cases:
                with self.subTest(path=path):
                    git('checkout', '-q', '--detach', base)
                    head = commit(path)
                    self.assertEqual(select(base, head), f'run={expected}')
            # A later unrelated commit must not hide fixture changes in the PR.
            later = commit('docs/example.md')
            self.assertEqual(select(base, later), 'run=true')
            self.assertEqual(select('', '', 'workflow_dispatch'), 'run=true')
            self.assertEqual(select('f' * 40, later), 'run=true')


if __name__ == '__main__':
    unittest.main()

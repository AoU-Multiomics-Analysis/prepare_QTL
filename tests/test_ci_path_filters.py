"""Check which expensive CI jobs a changed file selects (positive globs only)."""
from fnmatch import fnmatchcase
from pathlib import Path
import unittest

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
            'scripts/cell_type_specific_expression/run_hspe.R':
                {'build_and_push', 'build_cell_type_specific_expression', 'smoke'},
            'scripts/common/MergeCovariates.R': {'build_and_push', 'smoke'},
            'scripts/expression/merge_rnaseqc.py': {'build_and_push', 'container'},
            'rust/methylation_merge/src/main.rs': {'build_methylation_rust'},
            'envs/CellTypeSpecificExpression/environment.yml':
                {'build_cell_type_specific_expression', 'smoke'},
            'envs/PhenotypePCs/Dockerfile': {'build_and_push', 'smoke'},
            'envs/MethylationRust/Dockerfile': {'build_methylation_rust'},
            'envs/RNASeQCAggregation/environment.yml': {'container'},
            'tests/cell_type_specific_expression/fixtures/hspe-e2e.inputs.json': {'smoke'},
            'tests/test_prepare_expression_sample_list.R': {'build_and_push', 'smoke'},
            'tests/rnaseqc2_aggregation/smoke_container.py': {'container'},
            'tests/cell_type_specific_expression/test_reference_filter_wdl.py': set(),
            '.dockerignore': {'build_and_push', 'build_cell_type_specific_expression',
                              'build_methylation_rust', 'smoke', 'container'},
        }
        for event in ('push', 'pull_request'):
            for path, expected in cases.items():
                with self.subTest(event=event, path=path):
                    self.assertEqual(self.heavy_jobs(path, event), expected)

    def test_manual_dispatch_runs_all_builds(self):
        self.assertEqual(self.heavy_jobs('', 'workflow_dispatch'),
                         {'build_and_push', 'build_cell_type_specific_expression',
                          'build_methylation_rust', 'smoke', 'container'})


if __name__ == '__main__':
    unittest.main()

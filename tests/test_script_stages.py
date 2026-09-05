import copy
import importlib.util
from pathlib import Path
import subprocess
import tempfile
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]


class ScriptStagesTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        spec = importlib.util.spec_from_file_location(
            'script_stages', ROOT / 'ci/script_stages.py')
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)
        cls.config = yaml.safe_load((ROOT / 'ci/image-stages.yml').read_text())

    def test_source_patterns_adds_owned_directory_globs(self):
        stage = {
            'sources': ['scripts/legacy.R'],
            'script_roots': ['scripts/assay/prepare', 'rust/filter'],
        }
        self.assertEqual(
            self.module.source_patterns(stage),
            ['scripts/legacy.R', 'scripts/assay/prepare/**', 'rust/filter/**'],
        )

    def test_new_scripts_in_each_owned_root_select_the_owner(self):
        cases = {
            'scripts/cell_type_specific_expression/estimation/new_tool.R': {'cell_estimation'},
            'scripts/cell_type_specific_expression/fit/new_tool.R': {'cell_fit'},
            'scripts/cell_type_specific_expression/export/new_tool.R': {'cell_export'},
            'scripts/cell_type_specific_expression/downstream/new_tool.R': {'cell_downstream'},
            'scripts/expression/prepare/new_tool.R': {'expression'},
            'scripts/expression/rnaseqc/new_tool.py': {'rnaseqc'},
            'scripts/proteomics/new_tool.R': {'proteomics'},
            'scripts/splicing/new_tool.R': {'splicing'},
            'scripts/methylation/new_tool.R': {'methylation'},
            'rust/methylation_filter/src/new_module.rs': {'methylation_rust'},
            'rust/methylation_merge/src/new_module.rs': {'methylation_rust'},
        }
        for source, expected in cases.items():
            with self.subTest(source=source):
                self.assertEqual(
                    self.module.stages_for_script(self.config, source), expected)

    def test_new_shared_cell_script_selects_all_declared_consumers(self):
        self.assertEqual(
            self.module.stages_for_script(
                self.config,
                'scripts/cell_type_specific_expression/shared/new_helper.R',
            ),
            {'cell_estimation', 'cell_fit', 'cell_export', 'cell_downstream'},
        )

    def test_unknown_script_root_has_no_owner(self):
        self.assertEqual(
            self.module.stages_for_script(self.config, 'scripts/new_assay/tool.R'),
            set(),
        )

    def test_source_paths_must_be_safe_normalized_posix_paths(self):
        invalid_paths = [
            '/scripts/proteomics/tool.R',
            'scripts/proteomics/../splicing/tool.R',
            'scripts//proteomics/tool.R',
            'scripts/proteomics/./tool.R',
            'scripts\\proteomics\\tool.R',
        ]
        for source in invalid_paths:
            with self.subTest(source=source):
                with self.assertRaises(ValueError):
                    self.module.stages_for_script(self.config, source)

    def test_validation_rejects_overlapping_owned_roots(self):
        config = copy.deepcopy(self.config)
        config['stages']['splicing']['script_roots'] = ['scripts/proteomics/nested']
        errors = self.module.validate_script_roots(config)
        self.assertTrue(any('overlap' in error.lower() for error in errors), errors)

    def test_validation_rejects_unknown_shared_consumers(self):
        config = copy.deepcopy(self.config)
        config['shared'][0]['stages'].append('missing_stage')
        errors = self.module.validate_script_roots(config)
        self.assertTrue(any('missing_stage' in error for error in errors), errors)

    def test_validation_rejects_root_outside_its_image_build_paths(self):
        config = copy.deepcopy(self.config)
        config['stages']['rnaseqc']['script_roots'] = ['scripts/outside_rnaseqc']
        errors = self.module.validate_script_roots(config)
        self.assertTrue(any('build' in error.lower() for error in errors), errors)

    def test_validation_rejects_unsafe_roots_and_alias_targets(self):
        config = copy.deepcopy(self.config)
        config['stages']['proteomics']['script_roots'].append('scripts//unsafe')
        config['legacy_runtime_roots'].append('../outside')
        config['runtime_aliases'] = {'/tmp/Ambiguous.R': '../outside/Ambiguous.R'}
        errors = self.module.validate_script_roots(config)
        self.assertGreaterEqual(len(errors), 3, errors)

    def test_canonical_runtime_prefix_resolves_repository_source(self):
        self.assertEqual(
            self.module.source_for_runtime(
                self.config,
                '/opt/prepare_qtl/scripts/expression/merge_rnaseqc.py',
                ROOT,
            ),
            'scripts/expression/merge_rnaseqc.py',
        )

    def test_canonical_runtime_path_can_contain_spaces_and_quotes(self):
        self.assertEqual(
            self.module.source_for_runtime(
                self.config,
                '"/opt/prepare_qtl/scripts/proteomics/new tool.R"',
                ROOT,
            ),
            'scripts/proteomics/new tool.R',
        )

    def test_legacy_tmp_runtime_resolves_unique_tracked_standard_source(self):
        self.assertEqual(
            self.module.source_for_runtime(
                self.config, '/tmp/PrepareExpression.R', ROOT),
            'scripts/expression/PrepareExpression.R',
        )

    def test_legacy_tmp_runtime_rejects_unknown_source(self):
        with self.assertRaisesRegex(ValueError, 'Unknown runtime script'):
            self.module.source_for_runtime(self.config, '/tmp/NotTracked.R', ROOT)

    def test_legacy_tmp_runtime_rejects_ambiguous_basename(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp)
            subprocess.run(['git', 'init', '-q'], cwd=repo, check=True)
            for relative in ['scripts/one/Duplicate.R', 'scripts/two/Duplicate.R']:
                path = repo / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('# fixture\n')
            subprocess.run(
                ['git', 'add', 'scripts/one/Duplicate.R', 'scripts/two/Duplicate.R'],
                cwd=repo,
                check=True,
            )
            config = {
                'images': {
                    'standard': {'build_paths': ['scripts/**']},
                },
                'stages': {
                    'one': {
                        'image': 'standard',
                        'sources': [],
                        'script_roots': ['scripts/one'],
                    },
                    'two': {
                        'image': 'standard',
                        'sources': [],
                        'script_roots': ['scripts/two'],
                    },
                },
                'shared': [],
                'legacy_runtime_roots': ['scripts/one', 'scripts/two'],
            }
            with self.assertRaisesRegex(ValueError, 'Ambiguous runtime script'):
                self.module.source_for_runtime(config, '/tmp/Duplicate.R', repo)

    def test_explicit_runtime_alias_resolves_ambiguous_basename(self):
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp)
            subprocess.run(['git', 'init', '-q'], cwd=repo, check=True)
            for relative in ['scripts/one/Duplicate.R', 'scripts/two/Duplicate.R']:
                path = repo / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text('# fixture\n')
            subprocess.run(['git', 'add', 'scripts'], cwd=repo, check=True)
            config = {
                'images': {
                    'standard': {'build_paths': ['scripts/**']},
                },
                'stages': {
                    'one': {
                        'image': 'standard',
                        'sources': [],
                        'script_roots': ['scripts/one'],
                    },
                    'two': {
                        'image': 'standard',
                        'sources': [],
                        'script_roots': ['scripts/two'],
                    },
                },
                'shared': [],
                'legacy_runtime_roots': ['scripts/one', 'scripts/two'],
                'runtime_aliases': {
                    '/tmp/Duplicate.R': 'scripts/two/Duplicate.R',
                },
            }
            self.assertEqual(
                self.module.source_for_runtime(config, '/tmp/Duplicate.R', repo),
                'scripts/two/Duplicate.R',
            )

    def test_runtime_path_rejects_traversal_and_unknown_prefixes(self):
        paths = [
            '/opt/prepare_qtl/scripts/proteomics/../splicing/tool.R',
            '/opt/prepare_qtl/scripts//proteomics/tool.R',
            '/usr/local/bin/tool.R',
        ]
        for runtime_path in paths:
            with self.subTest(runtime_path=runtime_path):
                with self.assertRaises(ValueError):
                    self.module.source_for_runtime(
                        self.config, runtime_path, ROOT)


if __name__ == '__main__':
    unittest.main()

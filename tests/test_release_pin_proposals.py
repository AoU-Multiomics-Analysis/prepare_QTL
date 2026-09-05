"""Test stage-limited candidate patches without registry or GitHub writes."""
import importlib.util
from pathlib import Path
import sys
import unittest
import json
import subprocess
import tempfile
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'ci'))
OLD = 'ghcr.io/example/cell@sha256:' + 'a' * 64
NEW = 'ghcr.io/example/cell@sha256:' + 'b' * 64
TEXT = '''version 1.0
workflow Demo {
  input {
    String fit_image = "%s"
    String export_image = "%s"
  }
  output { String selected = export_image }
}
''' % (OLD, OLD)


class ProposalTests(unittest.TestCase):
    def setUp(self):
        path = ROOT / 'ci/propose_image_pins.py'
        self.assertTrue(path.exists(), 'Proposal tool is not implemented')
        spec = importlib.util.spec_from_file_location('proposals', path)
        self.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.module)
        self.config = {'images': {'cell': {'repository': 'ghcr.io/example/cell'}},
                       'stages': {'fit': {'image': 'cell'}, 'export': {'image': 'cell'}}}
        self.targets = {'version': 1, 'stages': {
            'fit': [{'path': 'workflows/a.wdl', 'input': 'fit_image'}],
            'export': [{'path': 'workflows/a.wdl', 'input': 'export_image'},
                       {'path': 'workflows/b.wdl', 'input': 'export_image'}]}}
        self.files = {'workflows/a.wdl': TEXT, 'workflows/b.wdl': TEXT}
        self.plan = {'stages': ['export'], 'unmapped': []}

    def propose(self, candidate=NEW):
        return self.module.propose(self.config, self.targets, self.plan,
                                   self.files, {'cell': candidate})

    def test_export_edit_preserves_fit_and_all_other_text(self):
        result = self.propose()
        expected = TEXT.replace('export_image = "' + OLD, 'export_image = "' + NEW)
        self.assertEqual(result, {'workflows/a.wdl': expected, 'workflows/b.wdl': expected})
        self.assertEqual(self.files['workflows/a.wdl'], TEXT)

    def test_same_digest_has_no_patch(self):
        self.assertEqual(self.propose(OLD), {})

    def test_unaffected_stages_have_no_patch(self):
        self.plan['stages'] = []
        self.assertEqual(self.propose(), {})

    def test_invalid_or_wrong_repository_reference_is_rejected(self):
        for value in ('ghcr.io/example/cell:main', NEW.replace('/cell@', '/other@'),
                      NEW + '\n', 'sha256:' + 'b' * 64):
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.propose(value)

    def test_duplicate_target_is_rejected(self):
        self.targets['stages']['fit'] += self.targets['stages']['export'][:1]
        with self.assertRaises(ValueError):
            self.propose()

    def test_missing_input_is_rejected(self):
        self.targets['stages']['export'][0]['input'] = 'absent'
        with self.assertRaises(ValueError):
            self.propose()

    def test_expression_instead_of_literal_default_is_rejected(self):
        self.files['workflows/a.wdl'] = TEXT.replace('"' + OLD + '"', '"prefix" + "suffix"')
        with self.assertRaises(ValueError):
            self.propose()

    def test_inconsistent_entry_point_defaults_are_rejected(self):
        self.files['workflows/b.wdl'] = TEXT.replace(OLD, NEW)
        with self.assertRaises(ValueError):
            self.propose()

    def test_unmapped_or_unsupported_stage_is_rejected(self):
        self.plan['unmapped'] = ['scripts/new.R']
        with self.assertRaises(ValueError):
            self.propose()
        self.plan['unmapped'] = []
        self.plan['stages'] = ['not_configured']
        with self.assertRaises(ValueError):
            self.propose()

    def test_missing_candidate_is_rejected(self):
        with self.assertRaises(ValueError):
            self.module.propose(self.config, self.targets, self.plan, self.files, {})

    def test_task_default_cannot_be_selected_as_workflow_input(self):
        self.files['workflows/a.wdl'] = TEXT.replace('String export_image =', 'String renamed =').replace(
            'selected = export_image', 'selected = renamed')
        with self.assertRaises(ValueError):
            self.propose()

    def test_real_repository_targets_update_only_downstream(self):
        config = yaml.safe_load((ROOT / 'ci/image-stages.yml').read_text())
        targets = yaml.safe_load((ROOT / 'ci/release-pins.yml').read_text())
        files = {item['path']: (ROOT / item['path']).read_text()
                 for locations in targets['stages'].values() for item in locations}
        candidate = config['images']['cell_type']['repository'] + '@sha256:' + 'c' * 64
        result = self.module.propose(config, targets, {'stages': ['cell_downstream'], 'unmapped': []},
                                     files, {'cell_type': candidate})
        self.assertEqual(len(result), 2)
        for path, updated in result.items():
            for name in ('estimation_docker_image', 'fit_docker_image', 'export_docker_image'):
                self.assertEqual(self.module.literal_span(updated, name)[2],
                                 self.module.literal_span(files[path], name)[2])
            self.assertEqual(self.module.literal_span(updated, 'downstream_docker_image')[2], candidate)


class ProposalCliTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.repo = Path(self.temp.name) / 'repo'
        self.repo.mkdir()
        self.git('init', '-q')
        self.git('config', 'user.email', 'test@example.invalid')
        self.git('config', 'user.name', 'Proposal Test')
        config = {'version': 1, 'images': {'cell': {'repository': 'ghcr.io/example/cell',
                    'build_paths': ['scripts/**'], 'environment_paths': []}},
                  'stages': {'export': {'image': 'cell', 'sources': ['scripts/export.R']}},
                  'shared': [], 'ignored_sources': [],
                  'workflow_groups': [{'paths': ['workflows/**'], 'stages': ['export']}]}
        targets = {'version': 1, 'stages': {'export': [
            {'path': 'workflows/demo.wdl', 'input': 'export_image'}]}}
        for name, value in [('ci/image-stages.yml', yaml.safe_dump(config)),
                            ('ci/release-pins.yml', yaml.safe_dump(targets)),
                            ('workflows/demo.wdl', TEXT), ('scripts/export.R', '# original\n')]:
            path = self.repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(value)
        self.git('add', '.')
        self.git('commit', '-qm', 'baseline')
        self.base = self.git('rev-parse', 'HEAD').strip()
        self.git('branch', 'release-base', self.base)
        (self.repo / 'scripts/export.R').write_text('# changed\n')
        self.git('add', '.')
        self.git('commit', '-qm', 'export change')
        self.head = self.git('rev-parse', 'HEAD').strip()
        self.output = Path(self.temp.name) / 'candidate'

    def git(self, *args):
        return subprocess.check_output(['git', '-C', str(self.repo), *args], text=True)

    def run_cli(self):
        return subprocess.run([sys.executable, str(ROOT / 'ci/propose_image_pins.py'),
            '--repo', str(self.repo), '--base', 'release-base', '--head', 'HEAD',
            '--expected-base', self.base, '--expected-head', self.head,
            '--candidate-image', 'cell=' + NEW, '--output-dir', str(self.output)],
            capture_output=True, text=True)

    def test_patch_applies_but_tracked_files_are_not_modified(self):
        result = self.run_cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        record = json.loads((self.output / 'release-candidate.json').read_text())
        self.assertEqual(record['plan']['stages'], ['export'])
        self.assertEqual(record['status'], 'candidate_not_validated')
        self.assertEqual(record['head'], self.head)
        self.git('apply', '--check', str(self.output / 'pins.patch'))
        self.assertEqual(self.git('status', '--porcelain'), '')
        self.assertEqual((self.repo / 'workflows/demo.wdl').read_text(), TEXT)
        self.assertNotEqual(self.run_cli().returncode, 0, 'Existing output must not be overwritten')

    def test_new_head_rejects_stale_candidate_before_output_creation(self):
        self.git('commit', '--allow-empty', '-qm', 'new source')
        result = self.run_cli()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('stale', result.stderr)
        self.assertFalse(self.output.exists())

    def test_uncommitted_config_is_not_used(self):
        (self.repo / 'ci/image-stages.yml').write_text('invalid: local change\n')
        result = self.run_cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual((self.repo / 'ci/image-stages.yml').read_text(), 'invalid: local change\n')

    def test_moved_base_rejects_stale_candidate(self):
        self.git('update-ref', 'refs/heads/release-base', self.head)
        result = self.run_cli()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(self.output.exists())

    def test_patch_handles_wdl_without_final_newline(self):
        (self.repo / 'workflows/demo.wdl').write_text(TEXT.rstrip('\n'))
        self.git('add', '.')
        self.git('commit', '-qm', 'no final newline')
        self.head = self.git('rev-parse', 'HEAD').strip()
        result = self.run_cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.git('apply', '--check', str(self.output / 'pins.patch'))


if __name__ == '__main__':
    unittest.main()

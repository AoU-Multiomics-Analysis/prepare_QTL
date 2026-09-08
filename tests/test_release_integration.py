"""Release policy tests; no credentials, registry writes, or container builds."""
import sys
from pathlib import Path
import unittest
from unittest import mock
import tempfile
import subprocess
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'ci'))


class ReleasePolicyTests(unittest.TestCase):
    def setUp(self):
        import release_integration as release
        self.release = release
        self.pr = {'number': 7, 'state': 'open',
                   'base': {'ref': 'main', 'sha': 'a' * 40, 'repo': {'full_name': 'owner/repo'}},
                   'head': {'ref': 'feature', 'sha': 'b' * 40, 'repo': {'full_name': 'owner/repo'}}}

    def test_accepts_only_open_same_repo_pr_against_main(self):
        self.release.validate_pr(self.pr, 'owner/repo')
        for field, replacement in [('state', 'closed'), ('head', {**self.pr['head'], 'repo': {'full_name': 'fork/repo'}}),
                                   ('base', {**self.pr['base'], 'ref': 'other'})]:
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.release.validate_pr({**self.pr, field: replacement}, 'owner/repo')

    def test_queued_release_skips_when_trigger_head_has_changed(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'release.json'
            with mock.patch.object(self.release, 'api', return_value=self.pr), \
                 mock.patch.object(self.release.subprocess, 'run') as execute:
                self.release.snapshot(Path(directory), 'owner/repo', 7, output,
                                      expected_head='c' * 40)
            self.assertFalse(output.exists())
            execute.assert_not_called()

    def test_rejects_stale_head_or_base(self):
        self.release.validate_pr(self.pr, 'owner/repo', self.pr)
        for field in ('base', 'head'):
            changed = {**self.pr, field: {**self.pr[field], 'sha': 'c' * 40}}
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.release.validate_pr(changed, 'owner/repo', self.pr)

    def test_source_changes_allowed_but_ci_policy_changes_rejected(self):
        self.release.check_policy_changes(['scripts/fit.R', 'workflows/fit.wdl', 'envs/Dockerfile', 'tests/test.R'])
        for path in ['ci/release-pins.yml', '.github/workflows/image-release.yml']:
            with self.subTest(path=path), self.assertRaises(ValueError):
                self.release.check_policy_changes([path])

    def test_publication_requires_known_source_fingerprint(self):
        result = {'fingerprint': 'a' * 64, 'reference': 'ghcr.io/owner/image@sha256:' + 'b' * 64}
        spec = {'fingerprint': 'a' * 64, 'repository': 'ghcr.io/owner/image'}
        self.release.validate_build_result(spec, result)
        with self.assertRaises(ValueError):
            self.release.validate_build_result(spec, {**result, 'fingerprint': 'c' * 64})
        with self.assertRaises(ValueError):
            self.release.validate_build_result(spec, {**result, 'reference': 'ghcr.io/owner/image:main'})

    def test_commit_payload_has_expected_head_and_only_supplied_files(self):
        payload = self.release.commit_payload('owner/repo', 'feature', 'a' * 40, {'workflows/a.wdl': 'version 1.0\n'})
        self.assertEqual(payload['expectedHeadOid'], 'a' * 40)
        self.assertEqual(payload['branch'], {'repositoryNameWithOwner': 'owner/repo', 'branchName': 'feature'})
        self.assertEqual([x['path'] for x in payload['fileChanges']['additions']], ['workflows/a.wdl'])
        self.assertNotIn('deletions', payload['fileChanges'])

    def test_release_workflow_separates_test_and_write_credentials(self):
        root = Path(__file__).resolve().parents[1]
        flow = yaml.load((root / '.github/workflows/image-release.yml').read_text(), Loader=yaml.BaseLoader)
        self.assertEqual(set(flow['on']), {'workflow_dispatch'})
        self.assertIn("github.ref == 'refs/heads/main'", flow['jobs']['snapshot']['if'])
        self.assertEqual(flow['jobs']['build']['environment'], 'release-publish')
        self.assertEqual(flow['jobs']['commit']['environment'], 'release-commit')
        self.assertEqual(flow['permissions']['contents'], 'read')
        self.assertNotIn('packages', flow['permissions'])
        self.assertNotIn('test', flow['jobs'])
        self.assertEqual(flow['jobs']['commit']['needs'], ['snapshot', 'build'])
        self.assertIn("needs.build.result == 'success'", flow['jobs']['commit']['if'])
        for name in ('source-unit-checks.yml', 'pinned-image-smoke.yml'):
            manual = yaml.load((root / '.github/workflows' / name).read_text(), Loader=yaml.BaseLoader)
            self.assertEqual(set(manual['on']), {'workflow_dispatch'})

    def test_automatic_dispatch_never_checks_out_pr_code(self):
        root = Path(__file__).resolve().parents[1]
        flow = yaml.load((root / '.github/workflows/image-release-dispatch.yml').read_text(), Loader=yaml.BaseLoader)
        job = flow['jobs']['dispatch']
        self.assertIn("vars.RELEASE_AUTO_TRIGGER == 'true'", job['if'])
        self.assertEqual(len(job['steps']), 1)
        self.assertTrue(job['steps'][0]['uses'].startswith('actions/github-script@'))
        self.assertIn("ref: 'main'", job['steps'][0]['with']['script'])


class FingerprintTests(unittest.TestCase):
    def test_only_image_build_inputs_change_fingerprint(self):
        import release_integration as release
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)
            def git(*args):
                return subprocess.check_output(['git', '-C', directory, *args], text=True).strip()
            git('init', '-q')
            git('config', 'user.name', 'Test')
            git('config', 'user.email', 'test@example.invalid')
            (repo / 'scripts').mkdir()
            (repo / 'scripts/fit.R').write_text('# fit\n')
            (repo / 'a.wdl').write_text('version 1.0\n')
            git('add', '.')
            git('commit', '-qm', 'first')
            first = release.source_fingerprint(repo, 'HEAD', ['scripts/**'])
            (repo / 'a.wdl').write_text('version 1.0\n# new image pin\n')
            git('add', '.')
            git('commit', '-qm', 'pin')
            self.assertEqual(first, release.source_fingerprint(repo, 'HEAD', ['scripts/**']))
            (repo / 'scripts/fit.R').write_text('# new code\n')
            git('add', '.')
            git('commit', '-qm', 'code')
            self.assertNotEqual(first, release.source_fingerprint(repo, 'HEAD', ['scripts/**']))


class ReleaseTestSelection(unittest.TestCase):
    def test_wdl_changes_do_not_select_builds_or_runtime_tests(self):
        from test_release_images import selected_stages
        from plan_image_updates import plan_changes
        root = Path(__file__).resolve().parents[1]
        config = yaml.safe_load((root / 'ci/image-stages.yml').read_text())
        plan = plan_changes(config, ['workflows/methylation/cohort_aggregation.wdl'])
        self.assertEqual(plan['builds'], [])
        self.assertEqual(selected_stages({'plan': plan}, config), set())

    def test_each_registered_stage_has_a_runtime_gate(self):
        from test_release_images import SUPPORTED_STAGES
        root = Path(__file__).resolve().parents[1]
        config = yaml.safe_load((root / 'ci/image-stages.yml').read_text())
        self.assertEqual({v.get('test_group', k) for k, v in config['stages'].items()}, SUPPORTED_STAGES)

    def test_release_runner_rejects_candidate_policy_override(self):
        import test_release_images as release_tests
        from wdl_stage_routing import validate_routing

        with tempfile.TemporaryDirectory() as directory:
            candidate = Path(directory) / 'candidate'
            trusted = Path(directory) / 'trusted'

            def write(root, relative, text):
                path = root / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)

            digest = 'ghcr.io/example/cell@sha256:' + 'a' * 64
            trusted_config = {
                'images': {'cell': {'repository': 'ghcr.io/example/cell'}},
                'stages': {
                    'cell_export': {
                        'image': 'cell',
                        'script_roots': ['scripts/cell_type_specific_expression/export'],
                    },
                    'cell_fit': {
                        'image': 'cell',
                        'script_roots': ['scripts/cell_type_specific_expression/fit'],
                    },
                },
                'shared': [],
            }
            pins = {'stages': {
                'cell_export': [{'path': 'workflows/main.wdl', 'input': 'export_docker_image'}],
                'cell_fit': [{'path': 'workflows/main.wdl', 'input': 'fit_docker_image'}],
            }}
            wdl = f'''version 1.0

task Export {{
  input {{ String docker_image }}
  command <<<
    Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/export/new.R
  >>>
  runtime {{ docker: docker_image }}
}}

workflow Main {{
  input {{
    String export_docker_image = "{digest}"
    String fit_docker_image = "{digest}"
  }}
  call Export {{ input: docker_image = fit_docker_image }}
}}
'''
            write(candidate, 'workflows/main.wdl', wdl)
            write(candidate, 'scripts/cell_type_specific_expression/export/new.R', '# fixture\n')
            write(trusted, 'ci/image-stages.yml', yaml.safe_dump(trusted_config))
            write(trusted, 'ci/release-pins.yml', yaml.safe_dump(pins))

            candidate_config = yaml.safe_load(yaml.safe_dump(trusted_config))
            candidate_config['stages']['cell_export']['script_roots'], \
                candidate_config['stages']['cell_fit']['script_roots'] = (
                    candidate_config['stages']['cell_fit']['script_roots'],
                    candidate_config['stages']['cell_export']['script_roots'],
                )
            write(candidate, 'ci/image-stages.yml', yaml.safe_dump(candidate_config))
            write(candidate, 'ci/release-pins.yml', yaml.safe_dump(pins))

            self.assertEqual(validate_routing(candidate, candidate), [])
            self.assertTrue(validate_routing(candidate, trusted))
            with mock.patch.object(release_tests, '__file__', str(trusted / 'ci/test_release_images.py')), \
                    mock.patch.object(sys, 'argv', ['test_release_images.py', '--source',
                                                   str(candidate), '--all-stages']), \
                    mock.patch.object(release_tests.subprocess, 'run') as run:
                with self.assertRaisesRegex(ValueError, 'stage routing'):
                    release_tests.main()
                run.assert_not_called()


if __name__ == '__main__':
    unittest.main()

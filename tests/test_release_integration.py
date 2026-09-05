"""Release policy tests; no credentials, registry writes, or container builds."""
import sys
from pathlib import Path
import unittest
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

    def test_rejects_stale_head_or_base(self):
        self.release.validate_pr(self.pr, 'owner/repo', self.pr)
        for field in ('base', 'head'):
            changed = {**self.pr, field: {**self.pr[field], 'sha': 'c' * 40}}
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.release.validate_pr(changed, 'owner/repo', self.pr)

    def test_source_changes_allowed_but_ci_policy_changes_rejected(self):
        self.release.check_policy_changes(['scripts/fit.R', 'workflows/fit.wdl', 'envs/Dockerfile'])
        for path in ['ci/release-pins.yml', '.github/workflows/image-release.yml', 'tests/test.R']:
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
        for step in flow['jobs']['test']['steps']:
            self.assertNotIn('create-github-app-token', step.get('uses', ''))
            self.assertNotIn('secrets.', str(step))
        self.assertIn("needs.test.result == 'success'", flow['jobs']['commit']['if'])

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
    def test_wdl_changes_select_consumers_without_builds(self):
        from test_release_images import selected_stages
        from plan_image_updates import plan_changes
        root = Path(__file__).resolve().parents[1]
        config = yaml.safe_load((root / 'ci/image-stages.yml').read_text())
        plan = plan_changes(config, ['workflows/methylation/cohort_aggregation.wdl'])
        self.assertEqual(plan['builds'], [])
        self.assertEqual(selected_stages({'plan': plan}, config), {'methylation', 'methylation_rust'})

    def test_each_registered_stage_has_a_runtime_gate(self):
        from test_release_images import SUPPORTED_STAGES
        root = Path(__file__).resolve().parents[1]
        config = yaml.safe_load((root / 'ci/image-stages.yml').read_text())
        self.assertEqual(set(config['stages']), SUPPORTED_STAGES)


if __name__ == '__main__':
    unittest.main()

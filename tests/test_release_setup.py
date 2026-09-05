"""Tests for the read-only release settings planner and GitHub transport."""
import base64
import copy
import contextlib
import io
import json
import os
import sys
from pathlib import Path
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'ci'))

from release_setup import (GitHubCLI, SetupError, apply_settings,
                           inspect_settings, plan_settings)
from setup_release import main


def empty_settings():
    return {
        'repository': {'full_name': 'owner/repo', 'id': 42,
                       'owner_type': 'Organization', 'admin': True,
                       'main_sha': 'a' * 40},
        'workflows_present': True,
        'reviewers': [{'type': 'User', 'id': 7, 'login': 'reviewer'}],
        'environments': {'release-publish': None, 'release-commit': None},
        'variables': {}, 'secret_present': False, 'label_present': False,
        'readiness_notices': [],
    }


def protected_environment(*, reviewers=None, branches=None, wait_timer=15, **extra):
    environment = {
        'wait_timer': wait_timer,
        'prevent_self_review': True,
        'reviewers': reviewers or [{'type': 'User', 'reviewer': {'id': 7}}],
        'deployment_branch_policy': {
            'protected_branches': False, 'custom_branch_policies': True,
        },
        **extra,
    }
    return {
        'environment': environment,
        'branches': branches or [{'name': 'main', 'type': 'branch'}],
    }


class SettingsPlannerTests(unittest.TestCase):
    def test_new_installation_never_enables_releases(self):
        plan = plan_settings(empty_settings(), None)
        values = {a['target']: a['value'] for a in plan['actions']
                  if a['kind'] == 'variable'}
        self.assertEqual(values, {'RELEASE_ENABLED': 'false',
                                  'RELEASE_AUTO_TRIGGER': 'false',
                                  'RELEASE_AUTO_COMMIT': 'false'})
        self.assertEqual(plan['conflicts'], [])

    def test_existing_true_variable_is_preserved(self):
        snapshot = empty_settings()
        snapshot['variables']['RELEASE_ENABLED'] = 'true'
        plan = plan_settings(snapshot, None)
        self.assertNotIn('RELEASE_ENABLED', {a['target'] for a in plan['actions']})
        self.assertTrue(any('RELEASE_ENABLED' in notice and 'true' in notice
                            for notice in plan['notices']))

    def test_different_app_id_is_a_conflict(self):
        snapshot = empty_settings()
        snapshot['variables']['RELEASE_APP_ID'] = '123'
        plan = plan_settings(snapshot, 456)
        self.assertTrue(any('RELEASE_APP_ID' in conflict for conflict in plan['conflicts']))
        self.assertFalse(any(a['target'] == 'RELEASE_APP_ID' for a in plan['actions']))

    def test_secret_without_app_id_is_a_conflict(self):
        snapshot = empty_settings()
        snapshot['secret_present'] = True
        plan = plan_settings(snapshot, None)
        self.assertTrue(any('secret' in conflict.lower() and 'RELEASE_APP_ID' in conflict
                            for conflict in plan['conflicts']))

    def test_repository_admin_and_workflows_are_required(self):
        for field in ('admin', 'workflows_present'):
            with self.subTest(field=field):
                snapshot = empty_settings()
                if field == 'admin':
                    snapshot['repository']['admin'] = False
                else:
                    snapshot[field] = False
                plan = plan_settings(snapshot, None)
                self.assertNotEqual(plan['conflicts'], [])
                self.assertEqual(plan['actions'], [])

    def test_one_to_six_explicit_user_reviewers_are_required(self):
        for reviewers in ([], [{'type': 'Team', 'id': 7, 'login': 'reviewer'}],
                          [{'type': 'User', 'id': i, 'login': f'user-{i}'}
                           for i in range(7)]):
            with self.subTest(reviewers=reviewers):
                snapshot = empty_settings()
                snapshot['reviewers'] = reviewers
                plan = plan_settings(snapshot, None)
                self.assertTrue(any('reviewer' in conflict.lower()
                                    for conflict in plan['conflicts']))
                self.assertEqual(plan['actions'], [])

    def test_new_environments_use_exact_protection_policy(self):
        plan = plan_settings(empty_settings(), None)
        environments = [a for a in plan['actions'] if a['kind'] == 'environment']
        branches = [a for a in plan['actions'] if a['kind'] == 'branch_rule']
        self.assertEqual([a['target'] for a in environments],
                         ['release-publish', 'release-commit'])
        self.assertEqual(environments[0]['expected'], None)
        self.assertEqual(environments[0]['value'], {
            'wait_timer': 0,
            'prevent_self_review': False,
            'reviewers': [{'type': 'User', 'id': 7}],
            'deployment_branch_policy': {
                'protected_branches': False, 'custom_branch_policies': True,
            },
        })
        self.assertEqual(branches[0]['value'], {'name': 'main', 'type': 'branch'})

    def test_existing_environment_preserves_wait_timer_and_unknown_rules(self):
        snapshot = empty_settings()
        snapshot['environments'] = {
            name: protected_environment(wait_timer=30,
                                        can_admins_bypass=False,
                                        custom_protection_rules=[{'id': 99}])
            for name in ('release-publish', 'release-commit')
        }
        plan = plan_settings(snapshot, None)
        self.assertEqual([a for a in plan['actions']
                          if a['kind'] in ('environment', 'branch_rule')], [])
        self.assertEqual(plan['conflicts'], [])

    def test_existing_environment_requires_custom_branch_policy_controls(self):
        snapshot = empty_settings()
        current = protected_environment()
        current['environment']['deployment_branch_policy'] = {
            'protected_branches': True, 'custom_branch_policies': False,
        }
        snapshot['environments']['release-publish'] = current
        plan = plan_settings(snapshot, None)
        self.assertTrue(any('release-publish' in conflict and 'deployment' in conflict
                            for conflict in plan['conflicts']))

    def test_direct_reviewer_ids_from_created_environment_are_reusable(self):
        snapshot = empty_settings()
        snapshot['environments'] = {
            name: protected_environment(reviewers=[{'type': 'User', 'id': 7}])
            for name in ('release-publish', 'release-commit')
        }
        snapshot['variables'] = {name: 'false' for name in (
            'RELEASE_ENABLED', 'RELEASE_AUTO_TRIGGER', 'RELEASE_AUTO_COMMIT')}
        snapshot['label_present'] = True
        plan = plan_settings(snapshot, None)
        self.assertEqual(plan['conflicts'], [])
        self.assertEqual(plan['actions'], [])

    def test_readiness_notices_pass_through_to_plan(self):
        snapshot = empty_settings()
        snapshot['readiness_notices'] = ['Package ghcr.io/owner/image is not public.']
        plan = plan_settings(snapshot, None)
        self.assertIn(snapshot['readiness_notices'][0], plan['notices'])

    def test_additional_branch_or_tag_policy_is_a_conflict(self):
        variants = [
            [{'name': 'main', 'type': 'branch'}, {'name': 'release/*', 'type': 'branch'}],
            [{'name': 'main', 'type': 'branch'}, {'name': 'v*', 'type': 'tag'}],
        ]
        for branches in variants:
            with self.subTest(branches=branches):
                snapshot = empty_settings()
                snapshot['environments']['release-publish'] = protected_environment(
                    branches=branches)
                plan = plan_settings(snapshot, None)
                self.assertTrue(any('release-publish' in conflict and 'main' in conflict
                                    for conflict in plan['conflicts']))

    def test_absent_main_branch_policy_is_a_conflict(self):
        snapshot = empty_settings()
        snapshot['environments']['release-publish'] = protected_environment(
            branches=[{'name': 'other', 'type': 'branch'}])
        plan = plan_settings(snapshot, None)
        self.assertTrue(any('release-publish' in conflict and 'main' in conflict
                            for conflict in plan['conflicts']))

    def test_mismatched_existing_reviewer_set_is_a_conflict(self):
        snapshot = empty_settings()
        snapshot['environments']['release-publish'] = protected_environment(
            reviewers=[{'type': 'User', 'reviewer': {'id': 8}}])
        plan = plan_settings(snapshot, None)
        self.assertTrue(any('release-publish' in conflict and 'reviewer' in conflict.lower()
                            for conflict in plan['conflicts']))

    def test_extra_or_malformed_environment_reviewer_is_a_conflict(self):
        invalid_reviewer_sets = [
            [
                {'type': 'User', 'reviewer': {'id': 7}},
                {'type': 'Team', 'reviewer': {'id': 88}},
            ],
            [
                {'type': 'User', 'reviewer': {'id': 7}},
                {'type': 'User', 'reviewer': {}},
            ],
        ]
        for reviewers in invalid_reviewer_sets:
            with self.subTest(reviewers=reviewers):
                snapshot = empty_settings()
                snapshot['environments']['release-publish'] = protected_environment(
                    reviewers=reviewers)
                plan = plan_settings(snapshot, None)
                self.assertTrue(any(
                    'release-publish' in conflict and 'reviewer' in conflict.lower()
                    for conflict in plan['conflicts']))

    def test_already_complete_snapshot_has_no_settings_writes(self):
        snapshot = empty_settings()
        snapshot['environments'] = {
            name: protected_environment()
            for name in ('release-publish', 'release-commit')
        }
        snapshot['variables'] = {
            'RELEASE_APP_ID': '123', 'RELEASE_ENABLED': 'false',
            'RELEASE_AUTO_TRIGGER': 'false', 'RELEASE_AUTO_COMMIT': 'false',
        }
        snapshot['secret_present'] = True
        snapshot['label_present'] = True
        plan = plan_settings(snapshot, 123)
        self.assertEqual(plan['actions'], [])
        self.assertEqual(plan['conflicts'], [])
        self.assertTrue(any('cannot be inspected' in notice for notice in plan['notices']))


def api_response(status, body):
    reason = 'OK' if status < 400 else 'Error'
    return (f'HTTP/2.0 {status} {reason}\r\ncontent-type: application/json\r\n\r\n'
            + json.dumps(body)).encode()


class GitHubCLITransportTests(unittest.TestCase):
    @mock.patch('release_setup.subprocess.run')
    def test_get_is_explicit_and_removes_inherited_debug(self, run):
        run.return_value = mock.Mock(returncode=0, stdout=api_response(200, {'ok': True}),
                                     stderr=b'')
        with mock.patch.dict(os.environ, {'GH_DEBUG': 'api'}, clear=False):
            result = GitHubCLI().get('/repos/owner/repo')
        command = run.call_args.args[0]
        options = run.call_args.kwargs
        self.assertEqual(command, ['gh', 'api', '--hostname', 'github.com',
                                   '--method', 'GET', '--include', '/repos/owner/repo'])
        self.assertNotIn('GH_DEBUG', options['env'])
        self.assertEqual(options['env']['GH_HOST'], 'github.com')
        self.assertEqual(result, {'ok': True})

    @mock.patch('release_setup.subprocess.run')
    def test_optional_404_requires_prior_repository_authorization(self, run):
        run.side_effect = [
            mock.Mock(returncode=0, stdout=api_response(200, {'id': 42}), stderr=b''),
            mock.Mock(returncode=1, stdout=api_response(404, {'message': 'Not Found'}),
                      stderr=b''),
        ]
        gh = GitHubCLI()
        gh.get('/repos/owner/repo')
        self.assertIsNone(gh.get('/repos/owner/repo/environments/missing', optional=True))

        run.reset_mock()
        run.side_effect = [mock.Mock(returncode=1,
                                     stdout=api_response(404, {'message': 'Not Found'}),
                                     stderr=b'')]
        with self.assertRaises(SetupError):
            GitHubCLI().get('/repos/owner/repo/environments/missing', optional=True)

    @mock.patch('release_setup.subprocess.run')
    def test_auth_rate_limit_and_malformed_responses_are_safe_errors(self, run):
        cases = [
            mock.Mock(returncode=1, stdout=api_response(401, {'token': 'secret-token'}),
                      stderr=b'secret-token'),
            mock.Mock(returncode=1, stdout=api_response(403, {'message': 'rate limited'}),
                      stderr=b'rate limited secret-token'),
            mock.Mock(returncode=0,
                      stdout=b'HTTP/2.0 200 OK\r\ncontent-type: application/json\r\n\r\n{bad',
                      stderr=b'secret-token'),
        ]
        for response in cases:
            with self.subTest(response=response):
                run.return_value = response
                with self.assertRaises(SetupError) as caught:
                    GitHubCLI().get('/repos/owner/repo')
                self.assertNotIn('secret-token', str(caught.exception))

    @mock.patch('release_setup.subprocess.run')
    def test_process_start_errors_are_redacted(self, run):
        run.side_effect = OSError('secret-token')
        with self.assertRaises(SetupError) as caught:
            GitHubCLI().get('/repos/owner/repo')
        self.assertNotIn('secret-token', str(caught.exception))

    def test_pages_uses_one_hundred_items_and_stops_after_short_page(self):
        gh = GitHubCLI()
        first = [{'id': i} for i in range(100)]
        with mock.patch.object(gh, 'get', side_effect=[{'items': first},
                                                       {'items': [{'id': 100}]}]) as get:
            self.assertEqual(len(gh.pages('/things?scope=x', 'items')), 101)
        self.assertEqual(get.call_args_list, [
            mock.call('/things?scope=x&per_page=100&page=1'),
            mock.call('/things?scope=x&per_page=100&page=2'),
        ])

    def test_pages_rejects_repeated_full_page(self):
        gh = GitHubCLI()
        page = {'items': [{'id': i} for i in range(100)]}
        with mock.patch.object(gh, 'get', side_effect=[page, copy.deepcopy(page)]):
            with self.assertRaises(SetupError):
                gh.pages('/things', 'items')

    @mock.patch('release_setup.subprocess.run')
    def test_write_sends_json_on_stdin(self, run):
        run.return_value = mock.Mock(returncode=0, stdout=api_response(201, {'id': 9}),
                                     stderr=b'')
        result = GitHubCLI().write('POST', '/repos/owner/repo/labels',
                                   {'name': 'release-ready'})
        command = run.call_args.args[0]
        self.assertEqual(command, ['gh', 'api', '--hostname', 'github.com',
                                   '--method', 'POST', '--include', '--input', '-',
                                   '/repos/owner/repo/labels'])
        self.assertEqual(json.loads(run.call_args.kwargs['input']),
                         {'name': 'release-ready'})
        self.assertEqual(result, {'id': 9})

    @mock.patch('release_setup.subprocess.run')
    def test_secret_upload_uses_stdin_and_redacts_failures(self, run):
        pem = b'private-key-material'
        run.return_value = mock.Mock(returncode=1, stdout=b'', stderr=pem)
        with self.assertRaises(SetupError) as caught:
            GitHubCLI().upload_key('owner/repo', pem)
        command = run.call_args.args[0]
        self.assertEqual(command, [
            'gh', 'secret', 'set', 'RELEASE_APP_PRIVATE_KEY',
            '--repo', 'owner/repo', '--env', 'release-commit',
        ])
        self.assertNotIn(pem.decode(), ' '.join(command))
        self.assertEqual(run.call_args.kwargs['input'], pem)
        self.assertEqual(run.call_args.kwargs['env']['GH_HOST'], 'github.com')
        self.assertNotIn(pem.decode(), str(caught.exception))

    @mock.patch('release_setup.urllib.request.urlopen')
    def test_public_manifest_reads_only_fixed_ghcr_host(self, urlopen):
        token = mock.MagicMock()
        token.__enter__.return_value = token
        token.read.return_value = b'{"token":"anonymous-token"}'
        tags = mock.MagicMock()
        tags.__enter__.return_value = tags
        tags.read.return_value = b'{"tags":["sha-123"]}'
        manifest = mock.MagicMock()
        manifest.__enter__.return_value = manifest
        manifest.read.return_value = b'{"schemaVersion":2}'
        urlopen.side_effect = [token, tags, manifest]
        self.assertTrue(GitHubCLI().public_manifest(
            'ghcr.io/owner/image'))
        requests = [call.args[0] for call in urlopen.call_args_list]
        self.assertEqual([request.full_url for request in requests], [
            ('https://ghcr.io/token?service=ghcr.io&'
             'scope=repository%3Aowner%2Fimage%3Apull'),
            'https://ghcr.io/v2/owner/image/tags/list?n=100',
            'https://ghcr.io/v2/owner/image/manifests/sha-123',
        ])
        self.assertTrue(all(request.host == 'ghcr.io' for request in requests))
        self.assertNotIn('Authorization', requests[0].headers)
        self.assertEqual(requests[1].headers['Authorization'],
                         'Bearer anonymous-token')


def workflow_content(text):
    return {'encoding': 'base64',
            'content': base64.b64encode(text.encode()).decode()}


PUBLISH_WORKFLOW = """
jobs:
  snapshot:
    if: vars.RELEASE_ENABLED == 'true'
  build:
    environment: release-publish
  commit:
    if: vars.RELEASE_AUTO_COMMIT == 'true'
    environment: release-commit
    steps:
      - with:
          app-id: ${{ vars.RELEASE_APP_ID }}
          private-key: ${{ secrets.RELEASE_APP_PRIVATE_KEY }}
"""

DISPATCH_WORKFLOW = """
jobs:
  dispatch:
    if: vars.RELEASE_ENABLED == 'true' && vars.RELEASE_AUTO_TRIGGER == 'true'
    steps:
      - env:
          COMMIT_PINS: ${{ vars.RELEASE_AUTO_COMMIT }}
        with:
          label: release-ready
          workflow_id: image-release.yml
"""


class InspectGitHub:
    def __init__(self, workflow_override=None):
        self.get_calls = []
        self.page_calls = []
        self.manifest_calls = []
        self.workflow_override = workflow_override or {}

    def get(self, route, *, optional=False):
        self.get_calls.append((route, optional))
        values = {
            '/repos/owner/repo': {
                'full_name': 'owner/repo', 'id': 42,
                'owner': {'type': 'Organization'}, 'permissions': {'admin': True},
            },
            '/repos/owner/repo/branches/main': {'commit': {'sha': 'a' * 40}},
            '/repos/owner/repo/contents/.github/workflows/image-release.yml?ref=' + 'a' * 40:
                workflow_content(PUBLISH_WORKFLOW),
            '/repos/owner/repo/contents/.github/workflows/image-release-dispatch.yml?ref=' + 'a' * 40:
                workflow_content(DISPATCH_WORKFLOW),
            '/repos/owner/repo/contents/ci/image-stages.yml?ref=' + 'a' * 40:
                workflow_content('images:\n  standard:\n    repository: ghcr.io/owner/image\n'),
            '/repos/owner/repo/collaborators/reviewer/permission': {
                'permission': 'read',
                'user': {'type': 'User', 'id': 7, 'login': 'reviewer'},
            },
            '/repos/owner/repo/environments/release-publish': None,
            '/repos/owner/repo/environments/release-commit': None,
            '/repos/owner/repo/labels/release-ready': None,
        }
        values.update(self.workflow_override)
        if route not in values:
            raise AssertionError(f'unexpected GET {route}')
        return copy.deepcopy(values[route])

    def pages(self, route, key):
        self.page_calls.append((route, key))
        values = {
            ('/repos/owner/repo/actions/variables', 'variables'): [],
            ('/repositories/42/environments/release-commit/secrets', 'secrets'): [],
        }
        if (route, key) not in values:
            raise AssertionError(f'unexpected pages call {(route, key)}')
        return copy.deepcopy(values[(route, key)])

    def public_manifest(self, repository):
        self.manifest_calls.append(repository)
        return True


class SettingsInspectionTests(unittest.TestCase):
    def test_inspection_resolves_main_once_and_returns_non_secret_snapshot(self):
        gh = InspectGitHub()
        snapshot = inspect_settings(gh, 'owner/repo', ['reviewer'])
        self.assertEqual(snapshot, empty_settings())
        branch_reads = [route for route, _ in gh.get_calls if route.endswith('/branches/main')]
        self.assertEqual(branch_reads, ['/repos/owner/repo/branches/main'])
        workflow_reads = [route for route, _ in gh.get_calls if '/contents/' in route]
        self.assertEqual(workflow_reads, [
            '/repos/owner/repo/contents/.github/workflows/image-release.yml?ref=' + 'a' * 40,
            '/repos/owner/repo/contents/.github/workflows/image-release-dispatch.yml?ref=' + 'a' * 40,
            '/repos/owner/repo/contents/ci/image-stages.yml?ref=' + 'a' * 40,
        ])
        self.assertEqual(gh.page_calls, [
            ('/repos/owner/repo/actions/variables', 'variables'),
        ])
        self.assertEqual(gh.manifest_calls, ['ghcr.io/owner/image'])

    def test_inspection_rejects_identity_change_and_non_reader(self):
        cases = [
            {'/repos/owner/repo': {
                'full_name': 'other/repo', 'id': 42,
                'owner': {'type': 'Organization'}, 'permissions': {'admin': True},
            }},
            {'/repos/owner/repo/collaborators/reviewer/permission': {
                'permission': 'none',
                'user': {'type': 'User', 'id': 7, 'login': 'reviewer'},
            }},
        ]
        for overrides in cases:
            with self.subTest(overrides=overrides), self.assertRaises(SetupError):
                inspect_settings(InspectGitHub(overrides), 'owner/repo', ['reviewer'])

    def test_inspection_marks_missing_or_incompatible_workflow(self):
        route = ('/repos/owner/repo/contents/.github/workflows/'
                 'image-release-dispatch.yml?ref=' + 'a' * 40)
        for value in (None, workflow_content('jobs: {}\n')):
            with self.subTest(value=value):
                snapshot = inspect_settings(InspectGitHub({route: value}),
                                            'owner/repo', ['reviewer'])
                self.assertFalse(snapshot['workflows_present'])

    def test_workflow_contract_tokens_in_comments_do_not_count(self):
        route = ('/repos/owner/repo/contents/.github/workflows/'
                 'image-release-dispatch.yml?ref=' + 'a' * 40)
        commented = workflow_content(
            'jobs: {}\n# RELEASE_ENABLED RELEASE_AUTO_TRIGGER RELEASE_AUTO_COMMIT\n'
            '# release-ready image-release.yml\n')
        snapshot = inspect_settings(InspectGitHub({route: commented}),
                                    'owner/repo', ['reviewer'])
        self.assertFalse(snapshot['workflows_present'])

    def test_inspection_accepts_wrapped_base64_content(self):
        route = ('/repos/owner/repo/contents/.github/workflows/'
                 'image-release.yml?ref=' + 'a' * 40)
        response = workflow_content(PUBLISH_WORKFLOW)
        response['content'] = '\n'.join(
            response['content'][index:index + 20]
            for index in range(0, len(response['content']), 20))
        snapshot = inspect_settings(InspectGitHub({route: response}),
                                    'owner/repo', ['reviewer'])
        self.assertTrue(snapshot['workflows_present'])

    def test_inaccessible_package_is_a_readiness_notice(self):
        gh = InspectGitHub()
        gh.public_manifest = mock.Mock(return_value=False)
        snapshot = inspect_settings(gh, 'owner/repo', ['reviewer'])
        self.assertTrue(any('ghcr.io/owner/image' in notice
                            for notice in snapshot['readiness_notices']))
        self.assertTrue(any('write' in notice.lower()
                            for notice in snapshot['readiness_notices']))

    def test_missing_image_registry_is_a_readiness_notice(self):
        route = '/repos/owner/repo/contents/ci/image-stages.yml?ref=' + 'a' * 40
        snapshot = inspect_settings(InspectGitHub({route: None}),
                                    'owner/repo', ['reviewer'])
        self.assertTrue(any('image-stages.yml' in notice
                            for notice in snapshot['readiness_notices']))

    def test_absent_commit_environment_does_not_list_its_secrets(self):
        class MissingEnvironmentGitHub(InspectGitHub):
            def pages(self, route, key):
                if route == '/repositories/42/environments/release-commit/secrets':
                    raise SetupError('The environment does not exist.')
                return super().pages(route, key)

        snapshot = inspect_settings(MissingEnvironmentGitHub(),
                                    'owner/repo', ['reviewer'])
        self.assertFalse(snapshot['secret_present'])


class StatefulGitHub(InspectGitHub):
    def __init__(self, *, corrupt_environment=False, created_reviewers=None):
        super().__init__()
        self.environments = {'release-publish': None, 'release-commit': None}
        self.variables = {}
        self.label = None
        self.writes = []
        self.corrupt_environment = corrupt_environment
        self.created_reviewers = created_reviewers

    def get(self, route, *, optional=False):
        for name in self.environments:
            if route == f'/repos/owner/repo/environments/{name}':
                value = self.environments[name]
                return copy.deepcopy(value['environment'] if value else None)
        if route == '/repos/owner/repo/labels/release-ready':
            return copy.deepcopy(self.label)
        return super().get(route, optional=optional)

    def pages(self, route, key):
        for name in self.environments:
            if route == f'/repos/owner/repo/environments/{name}/deployment-branch-policies':
                value = self.environments[name]
                return copy.deepcopy(value['branches'] if value else [])
        if route == '/repos/owner/repo/actions/variables':
            return [{'name': name, 'value': value}
                    for name, value in self.variables.items()]
        return super().pages(route, key)

    def write(self, method, route, body):
        self.writes.append((method, route, copy.deepcopy(body)))
        if route.endswith('/deployment-branch-policies'):
            name = route.split('/environments/', 1)[1].split('/', 1)[0]
            self.environments[name]['branches'].append(copy.deepcopy(body))
        elif '/environments/' in route:
            name = route.rsplit('/', 1)[1]
            stored = copy.deepcopy(body)
            if self.corrupt_environment:
                stored['reviewers'] = []
            if self.created_reviewers is not None:
                stored['reviewers'] = copy.deepcopy(self.created_reviewers)
            self.environments[name] = {'environment': stored, 'branches': []}
        elif route.endswith('/actions/variables'):
            self.variables[body['name']] = body['value']
        elif route.endswith('/labels'):
            self.label = copy.deepcopy(body)
        else:
            raise AssertionError(f'unexpected write {route}')
        return {}


class SettingsApplyTests(unittest.TestCase):
    def test_apply_creates_only_planned_settings_and_returns_fresh_snapshot(self):
        gh = StatefulGitHub()
        snapshot = empty_settings()
        plan = plan_settings(snapshot, None)
        completed = []
        result = apply_settings(gh, snapshot, plan, completed)
        self.assertEqual(result['variables'], {
            'RELEASE_ENABLED': 'false', 'RELEASE_AUTO_TRIGGER': 'false',
            'RELEASE_AUTO_COMMIT': 'false',
        })
        self.assertTrue(result['label_present'])
        self.assertEqual(result['environments']['release-publish']['branches'],
                         [{'name': 'main', 'type': 'branch'}])
        self.assertEqual([method for method, _, _ in gh.writes],
                         ['PUT', 'POST', 'PUT', 'POST', 'POST', 'POST', 'POST', 'POST'])
        self.assertEqual(completed, [
            'environment:release-publish', 'branch-rule:release-publish',
            'environment:release-commit', 'branch-rule:release-commit',
            'variable:RELEASE_ENABLED', 'variable:RELEASE_AUTO_TRIGGER',
            'variable:RELEASE_AUTO_COMMIT', 'label:release-ready',
        ])

    def test_apply_stops_if_a_target_changed_after_planning(self):
        gh = StatefulGitHub()
        snapshot = empty_settings()
        plan = plan_settings(snapshot, None)
        gh.variables['RELEASE_ENABLED'] = 'true'
        with self.assertRaises(SetupError):
            apply_settings(gh, snapshot, plan, [])
        self.assertFalse(any(route.endswith('/actions/variables')
                             for _, route, _ in gh.writes))

    def test_partial_environment_requires_manual_repair(self):
        gh = StatefulGitHub(corrupt_environment=True)
        completed = []
        with self.assertRaises(SetupError) as caught:
            apply_settings(gh, empty_settings(), plan_settings(empty_settings(), None),
                           completed)
        self.assertEqual(len(gh.writes), 1)
        self.assertEqual(completed, ['environment:release-publish'])
        self.assertIn('manual repair', str(caught.exception).lower())

    def test_apply_rejects_app_id_or_unplanned_targets_before_writes(self):
        gh = StatefulGitHub()
        plan = {
            'actions': [{'kind': 'variable', 'target': 'RELEASE_APP_ID',
                         'expected': None, 'value': '123'}],
            'conflicts': [], 'notices': [],
        }
        with self.assertRaises(SetupError):
            apply_settings(gh, empty_settings(), plan, [])
        self.assertEqual(gh.writes, [])

    def test_creation_rejects_extra_or_malformed_environment_reviewer(self):
        invalid_reviewer_sets = [
            [
                {'type': 'User', 'id': 7},
                {'type': 'Team', 'id': 88},
            ],
            [
                {'type': 'User', 'id': 7},
                {'type': 'User'},
            ],
        ]
        for reviewers in invalid_reviewer_sets:
            with self.subTest(reviewers=reviewers):
                gh = StatefulGitHub(created_reviewers=reviewers)
                completed = []
                with self.assertRaises(SetupError) as caught:
                    apply_settings(
                        gh, empty_settings(), plan_settings(empty_settings(), None),
                        completed)
                self.assertEqual(len(gh.writes), 1)
                self.assertEqual(completed, ['environment:release-publish'])
                self.assertIn('manual repair', str(caught.exception).lower())


class CLIGitHub(StatefulGitHub):
    """Persist fake API changes so repeated CLI runs use actual prior results."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.secret_present = False
        self.events = []
        self.upload_error = False
        self.app_id_error = False
        self.after_upload = None

    def pages(self, route, key):
        if route == '/repositories/42/environments/release-commit/secrets':
            return [{'name': 'RELEASE_APP_PRIVATE_KEY',
                     'created_at': '2026-09-05T00:00:00Z',
                     'updated_at': '2026-09-05T00:00:00Z'}] if self.secret_present else []
        return super().pages(route, key)

    def write(self, method, route, body):
        if body.get('name') == 'RELEASE_APP_ID' and self.app_id_error:
            raise SetupError('GitHub API write failed; no response details were retained.')
        result = super().write(method, route, body)
        self.events.append(('write', route, body.get('name')))
        return result

    def upload_key(self, repository, pem):
        if repository != 'owner/repo' or pem != b'FAKE-SECRET-MUST-NOT-ESCAPE':
            raise AssertionError('Wrong credential upload arguments')
        self.events.append(('upload',))
        if self.upload_error:
            raise SetupError('GitHub secret upload failed; the private key was not retained.')
        self.secret_present = True
        if self.after_upload:
            self.after_upload()


class CLIApp:
    def __init__(self, gh):
        self.gh = gh
        self.on_verify = None
        self.notices = ['Target repository access is verified; other selected repositories were not checked.']

    def register_app(self, repository, owner_type):
        if (repository, owner_type) != ('owner/repo', 'Organization'):
            raise AssertionError('Wrong registration owner')
        self.gh.events.append(('register',))
        return {'id': 123, 'slug': 'release-helper', 'pem': b'FAKE-SECRET-MUST-NOT-ESCAPE'}

    def verify_installation(self, repository, app_id, pem, *, newly_created=False):
        if (repository, app_id, pem) != ('owner/repo', 123, b'FAKE-SECRET-MUST-NOT-ESCAPE'):
            raise AssertionError('Wrong App verification arguments')
        self.gh.events.append(('verify', newly_created))
        if self.on_verify:
            self.on_verify()
        return {'id': 123, 'slug': 'release-helper', 'repository': repository,
                'permissions': {'contents': 'write', 'pull_requests': 'read', 'metadata': 'read'},
                'notices': self.notices}


class ReleaseSetupCLITests(unittest.TestCase):
    args = ['--repo', 'owner/repo', '--reviewer', 'reviewer']

    def invoke(self, args, gh, app_service=None):
        output = io.StringIO()
        with contextlib.redirect_stdout(output), contextlib.redirect_stderr(output):
            code = main(args, gh=gh, app_service=app_service)
        return code, output.getvalue()

    def test_dry_run_cannot_read_key_register_app_or_write(self):
        gh = StatefulGitHub()
        gh.write = mock.Mock(side_effect=AssertionError('GitHub write'))
        gh.upload_key = mock.Mock(side_effect=AssertionError('secret write'))
        forbidden = mock.Mock()
        forbidden.register_app.side_effect = AssertionError('registration')
        forbidden.verify_installation.side_effect = AssertionError('verification')
        with mock.patch('setup_release._read_private_key',
                        side_effect=AssertionError('key read')):
            for extra in ([], ['--create-app'],
                          ['--app-id', '123', '--private-key-file', '/not/read/in/dry-run']):
                with self.subTest(extra=extra):
                    code, output = self.invoke(self.args + extra, gh, forbidden)
                    self.assertEqual(code, 0)
                    self.assertIn('stage=plan status=complete', output)
        self.assertEqual(gh.writes, [])

    def test_invalid_arguments_fail_before_any_external_operation(self):
        invalid = [
            ['--repo', 'owner/repo'],
            ['--repo', 'owner/repo/extra', '--reviewer', 'reviewer'],
            ['--repo', 'owner/..', '--reviewer', 'reviewer'],
            ['--repo', 'owner/repo?x=1', '--reviewer', 'reviewer'],
            self.args + ['--reviewer', 'REVIEWER'],
            self.args + sum((['--reviewer', f'user-{i}'] for i in range(6)), []),
            self.args + ['--app-id', '123'],
            self.args + ['--private-key-file', '/unused'],
            self.args + ['--app-id', '0', '--private-key-file', '/unused'],
            self.args + ['--create-app', '--app-id', '123', '--private-key-file', '/unused'],
            self.args + ['--private-key', 'FAKE-SECRET-MUST-NOT-ESCAPE'],
        ]
        forbidden = mock.Mock()
        forbidden.get.side_effect = AssertionError('inspection')
        with mock.patch('setup_release._read_private_key', side_effect=AssertionError('key read')):
            for args in invalid:
                with self.subTest(args=args):
                    output = io.StringIO()
                    with contextlib.redirect_stderr(output), self.assertRaises(SystemExit) as caught:
                        main(args, gh=forbidden, app_service=forbidden)
                    self.assertEqual(caught.exception.code, 2)
                    self.assertNotIn('FAKE-SECRET-MUST-NOT-ESCAPE', output.getvalue())

    def test_help_does_not_construct_transports(self):
        with mock.patch('setup_release.GitHubCLI', side_effect=AssertionError('transport')):
            with contextlib.redirect_stdout(io.StringIO()), self.assertRaises(SystemExit) as caught:
                main(['--help'])
        self.assertEqual(caught.exception.code, 0)

    def apply(self, gh, app=None, *, create=False):
        credentials = (['--create-app'] if create else
                       ['--app-id', '123', '--private-key-file', '/fake/test-only.pem'])
        with mock.patch('setup_release._read_private_key', return_value=b'FAKE-SECRET-MUST-NOT-ESCAPE'):
            return self.invoke(self.args + credentials + ['--apply'], gh, app or CLIApp(gh))

    def test_apply_verifies_protections_and_app_before_secret_then_app_id(self):
        gh = CLIGitHub()
        code, output = self.apply(gh)
        self.assertEqual(code, 0)
        self.assertTrue(gh.secret_present)
        self.assertEqual(gh.variables, {'RELEASE_ENABLED': 'false',
                                      'RELEASE_AUTO_TRIGGER': 'false',
                                      'RELEASE_AUTO_COMMIT': 'false', 'RELEASE_APP_ID': '123'})
        self.assertEqual(gh.events[-3:], [('verify', False), ('upload',),
                         ('write', '/repos/owner/repo/actions/variables', 'RELEASE_APP_ID')])
        self.assertIn('stage=report status=complete', output)
        self.assertNotIn('FAKE-SECRET-MUST-NOT-ESCAPE', output)
        self.assertIn('protected workflow trial', output)

    def test_repeat_run_preserves_stored_secret_and_enabled_variables(self):
        gh = CLIGitHub()
        self.assertEqual(self.apply(gh)[0], 0)
        gh.variables['RELEASE_ENABLED'] = 'true'
        gh.variables['RELEASE_AUTO_TRIGGER'] = 'true'
        gh.variables['RELEASE_AUTO_COMMIT'] = 'true'
        before = copy.deepcopy(gh.writes)
        gh.events.clear()
        code, output = self.apply(gh)
        self.assertEqual(code, 0)
        self.assertEqual(gh.writes, before)
        self.assertEqual(gh.events, [('verify', False)])
        self.assertIn('stored secret was not verified', output)
        for name in ('RELEASE_ENABLED', 'RELEASE_AUTO_TRIGGER', 'RELEASE_AUTO_COMMIT'):
            self.assertEqual(gh.variables[name], 'true')
            self.assertIn(f'{name} is already true', output)

    def test_create_waits_for_explicit_installation_continuation(self):
        gh = CLIGitHub()
        def continue_installation(prompt):
            gh.events.append(('continue',))
            return 'continue'
        with mock.patch('builtins.input', side_effect=continue_installation):
            code, output = self.apply(gh, create=True)
        self.assertEqual(code, 0)
        self.assertEqual(gh.events[-5:], [('register',), ('continue',), ('verify', True),
                                         ('upload',), ('write', '/repos/owner/repo/actions/variables', 'RELEASE_APP_ID')])
        self.assertIn('App ID: 123', output)
        self.assertIn('https://github.com/apps/release-helper/installations/new', output)
        self.assertIn('owner/repo only', output)

    def test_create_without_continuation_preserves_settings_without_credentials(self):
        for response in ('', 'no', EOFError()):
            gh = CLIGitHub()
            options = {'side_effect': response} if isinstance(response, Exception) else {'return_value': response}
            with self.subTest(response=response), mock.patch('builtins.input', **options):
                code, output = self.apply(gh, create=True)
            self.assertEqual(code, 1)
            self.assertFalse(any(event[0] in ('verify', 'upload') for event in gh.events))
            self.assertNotIn('RELEASE_APP_ID', gh.variables)
            self.assertIn('App ID: 123', output)
            self.assertNotIn('stage=report status=complete', output)

    def test_apply_without_credentials_or_create_with_credentials_stops_before_writes(self):
        for app_id, secret, create in ((None, False, False), ('123', False, True),
                                      ('123', True, True), (None, True, True)):
            gh = CLIGitHub()
            if app_id:
                gh.variables['RELEASE_APP_ID'] = app_id
            gh.secret_present = secret
            if secret:
                gh.environments['release-commit'] = protected_environment()
            args = self.args + ['--apply'] + (['--create-app'] if create else [])
            with self.subTest(app_id=app_id, secret=secret, create=create):
                code, output = self.invoke(args, gh, CLIApp(gh))
                self.assertEqual(code, 1)
                self.assertEqual(gh.writes, [])
                self.assertEqual(gh.events, [])
                self.assertNotIn('stage=report status=complete', output)

    def test_mismatched_app_or_environment_stops_before_key_read(self):
        for target in ('app', 'environment'):
            gh = CLIGitHub()
            if target == 'app':
                gh.variables['RELEASE_APP_ID'] = '456'
            else:
                gh.environments['release-commit'] = protected_environment(
                    reviewers=[{'type': 'User', 'id': 88}])
            with self.subTest(target=target), mock.patch('setup_release._read_private_key',
                                                        side_effect=AssertionError('key read')):
                code, _ = self.invoke(self.args + ['--app-id', '123', '--private-key-file', '/never/read', '--apply'], gh, CLIApp(gh))
            self.assertEqual(code, 1)
            self.assertEqual(gh.writes, [])

    def test_apply_rechecks_all_settings_before_first_write(self):
        gh = CLIGitHub()
        get = gh.get
        reads = 0
        def changed(route, **kwargs):
            nonlocal reads
            if route == '/repos/owner/repo':
                reads += 1
                if reads == 2:
                    gh.variables['RELEASE_ENABLED'] = 'true'
            return get(route, **kwargs)
        gh.get = changed
        code, _ = self.apply(gh)
        self.assertEqual(code, 1)
        self.assertEqual(gh.writes, [])

    def test_stale_environment_app_id_or_new_secret_stops_before_upload(self):
        for target in ('environment', 'app', 'secret'):
            gh = CLIGitHub()
            app = CLIApp(gh)
            def change():
                if target == 'environment':
                    gh.environments['release-commit']['environment']['reviewers'] = []
                elif target == 'app':
                    gh.variables['RELEASE_APP_ID'] = '456'
                else:
                    gh.variables['RELEASE_APP_ID'] = '123'
                    gh.secret_present = True
            app.on_verify = change
            with self.subTest(target=target):
                code, _ = self.apply(gh, app)
            self.assertEqual(code, 1)
            self.assertNotIn(('upload',), gh.events)
            self.assertFalse(any(body.get('name') == 'RELEASE_APP_ID' for _, _, body in gh.writes))

    def test_upload_failure_reports_replacement_key_recovery_and_public_id(self):
        gh = CLIGitHub()
        gh.upload_error = True
        code, output = self.apply(gh)
        self.assertEqual(code, 1)
        self.assertNotIn('RELEASE_APP_ID', gh.variables)
        self.assertIn('App ID: 123', output)
        self.assertIn('Generate a replacement private key', output)
        self.assertIn('--app-id', output)
        self.assertIn('environment:release-publish', output)
        self.assertNotIn('FAKE-SECRET-MUST-NOT-ESCAPE', output)
        gh.upload_error = False
        writes = len(gh.writes)
        self.assertEqual(self.apply(gh)[0], 0)
        self.assertEqual(len(gh.writes), writes + 1)

    def test_app_id_write_failure_requires_manual_metadata_repair(self):
        gh = CLIGitHub()
        gh.app_id_error = True
        code, output = self.apply(gh)
        self.assertEqual(code, 1)
        self.assertTrue(gh.secret_present)
        self.assertNotIn('RELEASE_APP_ID', gh.variables)
        self.assertIn('App ID: 123', output)
        self.assertIn('manual metadata repair', output)
        self.assertNotIn('FAKE-SECRET-MUST-NOT-ESCAPE', output)
        before = list(gh.events)
        self.assertEqual(self.apply(gh)[0], 1)
        self.assertEqual(gh.events, before)

    def test_metadata_readback_failure_is_not_reported_as_complete(self):
        gh = CLIGitHub()
        gh.after_upload = lambda: setattr(gh, 'secret_present', False)
        code, output = self.apply(gh)
        self.assertEqual(code, 1)
        self.assertNotIn('stage=report status=complete', output)

    def test_app_id_changed_after_upload_is_never_replaced(self):
        gh = CLIGitHub()
        gh.after_upload = lambda: gh.variables.update(RELEASE_APP_ID='456')
        code, _ = self.apply(gh)
        self.assertEqual(code, 1)
        self.assertEqual(gh.variables['RELEASE_APP_ID'], '456')
        self.assertFalse(any(body.get('name') == 'RELEASE_APP_ID' for _, _, body in gh.writes))

    def test_partial_environment_and_missing_installation_do_not_upload(self):
        gh = CLIGitHub(corrupt_environment=True)
        code, output = self.apply(gh)
        self.assertEqual(code, 1)
        self.assertEqual(len(gh.writes), 1)
        self.assertIn('manual repair', output)
        self.assertIn('environment:release-publish', output)
        gh = CLIGitHub()
        app = CLIApp(gh)
        def missing():
            raise SetupError('GitHub App API access failed; no response details were retained.')
        app.on_verify = missing
        code, output = self.apply(gh, app)
        self.assertEqual(code, 1)
        self.assertNotIn(('upload',), gh.events)
        self.assertNotIn('stage=report status=complete', output)

    def test_package_and_broad_app_access_notices_are_reported(self):
        gh = CLIGitHub()
        gh.public_manifest = lambda repository: False
        app = CLIApp(gh)
        app.notices = ['The existing App installation can access all owner repositories.',
                       'The existing App has additional write permissions; access was not changed.']
        code, output = self.apply(gh, app)
        self.assertEqual(code, 0)
        self.assertIn('ghcr.io/owner/image', output)
        self.assertIn('all owner repositories', output)
        self.assertIn('additional write permissions', output)

    def test_raw_os_error_is_redacted_and_interrupt_is_not_success(self):
        for error in (OSError('FAKE-SECRET-MUST-NOT-ESCAPE'), KeyboardInterrupt()):
            gh = CLIGitHub()
            app = CLIApp(gh)
            def fail():
                raise error
            app.on_verify = fail
            with self.subTest(error=type(error).__name__):
                if isinstance(error, KeyboardInterrupt):
                    with self.assertRaises(KeyboardInterrupt):
                        self.apply(gh, app)
                else:
                    code, output = self.apply(gh, app)
                    self.assertEqual(code, 1)
                    self.assertNotIn('FAKE-SECRET-MUST-NOT-ESCAPE', output)


if __name__ == '__main__':
    unittest.main()

"""Plan and apply the repository settings used by image releases."""
import base64
import json
import os
import re
import subprocess
import urllib.error
import urllib.parse
import urllib.request

import yaml


class SetupError(Exception):
    """A fixed, non-sensitive message suitable for the CLI."""


ENVIRONMENTS = ('release-publish', 'release-commit')
DISABLED_VARIABLES = (
    'RELEASE_ENABLED',
    'RELEASE_AUTO_TRIGGER',
    'RELEASE_AUTO_COMMIT',
)


def _safe_environment():
    environment = os.environ.copy()
    environment.pop('GH_DEBUG', None)
    environment['GH_HOST'] = 'github.com'
    return environment


def _parse_api_response(process):
    output = process.stdout or b''
    if isinstance(output, str):
        output = output.encode()
    match = re.match(rb'HTTP/\S+\s+(\d{3})(?:\s[^\r\n]*)?[\r\n]', output)
    if not match:
        raise SetupError('GitHub returned an invalid API response.')
    status = int(match.group(1))
    separator = b'\r\n\r\n' if b'\r\n\r\n' in output else b'\n\n'
    parts = output.split(separator, 1)
    if len(parts) != 2:
        raise SetupError('GitHub returned an invalid API response.')
    body = parts[1].strip()
    try:
        value = json.loads(body) if body else {}
    except (UnicodeDecodeError, json.JSONDecodeError):
        raise SetupError('GitHub returned an invalid API response.') from None
    return status, value


class GitHubCLI:
    """Small, GitHub.com-only transport around the authenticated GitHub CLI."""

    def __init__(self):
        self._authorized_repositories = set()

    def _api(self, command, *, body=None):
        try:
            process = subprocess.run(
                command,
                input=body,
                capture_output=True,
                check=False,
                env=_safe_environment(),
            )
        except OSError:
            raise SetupError('GitHub CLI could not start.') from None
        status, value = _parse_api_response(process)
        return process.returncode, status, value

    def get(self, route: str, *, optional: bool = False):
        command = ['gh', 'api', '--hostname', 'github.com', '--method', 'GET',
                   '--include', route]
        returncode, status, value = self._api(command)
        repository_match = re.fullmatch(r'/repos/([^/]+/[^/?]+)', route)
        resource_match = re.match(r'/repos/([^/]+/[^/?]+)(?:/|\?)', route)
        if returncode == 0 and 200 <= status < 300:
            if repository_match:
                self._authorized_repositories.add(repository_match.group(1).lower())
            return value
        if (optional and status == 404 and resource_match
                and resource_match.group(1).lower() in self._authorized_repositories):
            return None
        raise SetupError('GitHub API access failed; check authentication and permissions.')

    def pages(self, route: str, key: str) -> list:
        values = []
        signatures = set()
        separator = '&' if '?' in route else '?'
        page_number = 1
        while True:
            response = self.get(
                f'{route}{separator}per_page=100&page={page_number}')
            if not isinstance(response, dict) or not isinstance(response.get(key), list):
                raise SetupError('GitHub returned an invalid paginated response.')
            page = response[key]
            signature = json.dumps(page, sort_keys=True, separators=(',', ':'))
            if signature in signatures:
                raise SetupError('GitHub returned a repeated pagination page.')
            signatures.add(signature)
            values.extend(page)
            if len(page) < 100:
                return values
            page_number += 1

    def write(self, method: str, route: str, body: dict):
        command = ['gh', 'api', '--hostname', 'github.com', '--method', method,
                   '--include', '--input', '-', route]
        encoded = json.dumps(body, separators=(',', ':')).encode()
        returncode, status, value = self._api(command, body=encoded)
        if returncode != 0 or not 200 <= status < 300:
            raise SetupError('GitHub API write failed; no response details were retained.')
        return value

    def upload_key(self, repository: str, pem: bytes) -> None:
        try:
            process = subprocess.run(
                ['gh', 'secret', 'set', 'RELEASE_APP_PRIVATE_KEY',
                 '--repo', repository, '--env', 'release-commit'],
                input=pem,
                capture_output=True,
                check=False,
                env=_safe_environment(),
            )
        except OSError:
            raise SetupError('GitHub secret upload could not start.') from None
        if process.returncode != 0:
            raise SetupError('GitHub secret upload failed; the private key was not retained.')

    def public_manifest(self, repository: str) -> bool:
        """Return whether one anonymous manifest is readable from fixed-host GHCR."""
        match = re.fullmatch(r'ghcr\.io/([a-z0-9][a-z0-9._/-]*)', repository)
        if not match or '..' in match.group(1).split('/'):
            raise SetupError('The image registry contains an invalid GHCR repository.')
        package = match.group(1)

        def read_json(url, accept, token=None):
            headers = {'Accept': accept, 'User-Agent': 'release-setup'}
            if token is not None:
                headers['Authorization'] = f'Bearer {token}'
            request = urllib.request.Request(
                url, headers=headers)
            try:
                with urllib.request.urlopen(request, timeout=10) as response:
                    return json.loads(response.read())
            except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError,
                    UnicodeDecodeError, json.JSONDecodeError, ValueError):
                return None

        token_query = urllib.parse.urlencode({
            'service': 'ghcr.io',
            'scope': f'repository:{package}:pull',
        })
        token_response = read_json(
            f'https://ghcr.io/token?{token_query}', 'application/json')
        if not isinstance(token_response, dict):
            return False
        token = token_response.get('token') or token_response.get('access_token')
        if not isinstance(token, str) or not token:
            return False
        tags = read_json(
            f'https://ghcr.io/v2/{package}/tags/list?n=100',
            'application/json', token)
        if not isinstance(tags, dict) or not isinstance(tags.get('tags'), list) or not tags['tags']:
            return False
        tag = tags['tags'][0]
        if not isinstance(tag, str) or not tag:
            return False
        manifest = read_json(
            f'https://ghcr.io/v2/{package}/manifests/{urllib.parse.quote(tag, safe="")}',
            ', '.join((
                'application/vnd.oci.image.index.v1+json',
                'application/vnd.oci.image.manifest.v1+json',
                'application/vnd.docker.distribution.manifest.list.v2+json',
                'application/vnd.docker.distribution.manifest.v2+json',
            )), token)
        return isinstance(manifest, dict) and manifest.get('schemaVersion') == 2


def _decode_yaml(response):
    if not isinstance(response, dict) or response.get('encoding') != 'base64':
        return None
    try:
        encoded = ''.join(response['content'].split())
        text = base64.b64decode(encoded, validate=True).decode()
        return yaml.safe_load(text)
    except (KeyError, ValueError, UnicodeDecodeError, yaml.YAMLError):
        return None


def _decode_workflow(response):
    document = _decode_yaml(response)
    if not isinstance(document, dict):
        return None
    return document


def _workflows_match(publish, dispatch):
    if not isinstance(publish, dict) or not isinstance(dispatch, dict):
        return False
    publish_jobs = publish.get('jobs') or {}
    snapshot_job = publish_jobs.get('snapshot') or {}
    build_job = publish_jobs.get('build') or {}
    commit_job = publish_jobs.get('commit') or {}
    dispatch_job = (dispatch.get('jobs') or {}).get('dispatch') or {}
    commit_text = json.dumps(commit_job, sort_keys=True, separators=(',', ':'))
    dispatch_text = json.dumps(dispatch_job, sort_keys=True, separators=(',', ':'))
    return (
        'RELEASE_ENABLED' in str(snapshot_job.get('if', ''))
        and build_job.get('environment') == 'release-publish'
        and commit_job.get('environment') == 'release-commit'
        and 'RELEASE_AUTO_COMMIT' in str(commit_job.get('if', ''))
        and 'RELEASE_APP_ID' in commit_text
        and 'RELEASE_APP_PRIVATE_KEY' in commit_text
        and 'RELEASE_ENABLED' in str(dispatch_job.get('if', ''))
        and 'RELEASE_AUTO_TRIGGER' in str(dispatch_job.get('if', ''))
        and 'release-ready' in dispatch_text
        and 'RELEASE_AUTO_COMMIT' in dispatch_text
        and 'image-release.yml' in dispatch_text
    )


def _normalize_environment(environment):
    normalized = dict(environment)
    if 'reviewers' not in normalized:
        normalized['reviewers'] = []
    if 'wait_timer' not in normalized:
        normalized['wait_timer'] = 0
    if 'prevent_self_review' not in normalized:
        normalized['prevent_self_review'] = False
    for rule in environment.get('protection_rules') or []:
        if rule.get('type') == 'wait_timer':
            normalized['wait_timer'] = rule.get('wait_timer', 0)
        elif rule.get('type') == 'required_reviewers':
            normalized['reviewers'] = rule.get('reviewers') or []
            normalized['prevent_self_review'] = rule.get('prevent_self_review', False)
    return normalized


def _branch_specs(branches):
    return [{'name': branch.get('name'), 'type': branch.get('type')}
            for branch in branches]


def inspect_settings(gh, repository: str, reviewers: list[str]) -> dict:
    """Read and normalize the non-sensitive release settings."""
    if not re.fullmatch(r'[^/\s]+/[^/\s]+', repository):
        raise SetupError('Repository must use the OWNER/REPO form.')
    if not 1 <= len(reviewers) <= 6 or len({name.lower() for name in reviewers}) != len(reviewers):
        raise SetupError('Specify one to six unique GitHub reviewer logins.')

    repo = gh.get(f'/repos/{repository}')
    full_name = repo.get('full_name')
    if not isinstance(full_name, str) or full_name.lower() != repository.lower():
        raise SetupError('GitHub returned a different repository identity.')
    repository_id = repo.get('id')
    owner_type = (repo.get('owner') or {}).get('type')
    if not isinstance(repository_id, int) or owner_type not in ('User', 'Organization'):
        raise SetupError('GitHub returned incomplete repository identity data.')

    branch = gh.get(f'/repos/{full_name}/branches/main')
    main_sha = ((branch.get('commit') or {}).get('sha')
                if isinstance(branch, dict) else None)
    if not isinstance(main_sha, str) or not re.fullmatch(r'[0-9a-fA-F]{40}', main_sha):
        raise SetupError('GitHub returned an invalid main branch revision.')

    workflow_texts = []
    for name in ('image-release.yml', 'image-release-dispatch.yml'):
        response = gh.get(
            f'/repos/{full_name}/contents/.github/workflows/{name}?ref={main_sha}',
            optional=True)
        workflow_texts.append(_decode_workflow(response) if response else None)

    readiness_notices = []
    registry_response = gh.get(
        f'/repos/{full_name}/contents/ci/image-stages.yml?ref={main_sha}',
        optional=True)
    registry = _decode_yaml(registry_response) if registry_response else None
    if not isinstance(registry, dict) or not isinstance(registry.get('images'), dict):
        readiness_notices.append(
            'ci/image-stages.yml is missing or invalid; GHCR readiness was not checked.')
    else:
        image_repositories = set()
        for image in registry['images'].values():
            if isinstance(image, dict) and isinstance(image.get('repository'), str):
                image_repositories.add(image['repository'])
            else:
                readiness_notices.append(
                    'ci/image-stages.yml has an image without a valid repository.')
        for image_repository in sorted(image_repositories):
            if not gh.public_manifest(image_repository):
                readiness_notices.append(
                    f'A public manifest for {image_repository} could not be read. '
                    'Resolve package access; this read-only check does not prove package-write access.')

    resolved_reviewers = []
    readable_permissions = {'read', 'triage', 'write', 'maintain', 'admin'}
    for login in reviewers:
        permission = gh.get(
            f'/repos/{full_name}/collaborators/{login}/permission')
        user = permission.get('user') or {}
        if (permission.get('permission') not in readable_permissions
                or user.get('type') != 'User'
                or not isinstance(user.get('id'), int)
                or str(user.get('login', '')).lower() != login.lower()):
            raise SetupError('A requested reviewer is not a readable repository user.')
        resolved_reviewers.append({
            'type': 'User', 'id': user['id'], 'login': user['login'],
        })

    environments = {}
    for name in ENVIRONMENTS:
        environment = gh.get(f'/repos/{full_name}/environments/{name}', optional=True)
        if environment is None:
            environments[name] = None
        else:
            branches = gh.pages(
                f'/repos/{full_name}/environments/{name}/deployment-branch-policies',
                'branch_policies')
            environments[name] = {
                'environment': _normalize_environment(environment),
                'branches': _branch_specs(branches),
            }

    variable_records = gh.pages(f'/repos/{full_name}/actions/variables', 'variables')
    variables = {}
    for record in variable_records:
        name = record.get('name')
        value = record.get('value')
        if isinstance(name, str) and isinstance(value, str):
            variables[name] = value
        else:
            raise SetupError('GitHub returned invalid repository variable data.')
    secret_records = gh.pages(
        f'/repositories/{repository_id}/environments/release-commit/secrets', 'secrets')
    secret_present = any(
        record.get('name') == 'RELEASE_APP_PRIVATE_KEY' for record in secret_records)
    label_present = gh.get(
        f'/repos/{full_name}/labels/release-ready', optional=True) is not None

    return {
        'repository': {
            'full_name': full_name,
            'id': repository_id,
            'owner_type': owner_type,
            'admin': bool((repo.get('permissions') or {}).get('admin')),
            'main_sha': main_sha,
        },
        'workflows_present': _workflows_match(*workflow_texts),
        'reviewers': resolved_reviewers,
        'environments': environments,
        'variables': variables,
        'secret_present': secret_present,
        'label_present': label_present,
        'readiness_notices': readiness_notices,
    }


def _reviewer_ids(reviewers):
    return {reviewer['id'] for reviewer in reviewers}


def _environment_reviewer_ids(environment):
    ids = set()
    for entry in environment.get('reviewers') or []:
        reviewer = entry.get('reviewer') or {}
        reviewer_id = entry.get('id', reviewer.get('id'))
        if entry.get('type') == 'User' and isinstance(reviewer_id, int):
            ids.add(reviewer_id)
    return ids


def plan_settings(snapshot: dict, requested_app_id: int | None) -> dict:
    """Return non-secret setting actions without changing existing settings."""
    actions = []
    conflicts = []
    notices = list(snapshot.get('readiness_notices') or [])
    repository = snapshot.get('repository') or {}
    reviewers = snapshot.get('reviewers') or []

    if not repository.get('admin'):
        conflicts.append('Repository administrator access is required.')
    if not snapshot.get('workflows_present'):
        conflicts.append('Required release workflows are missing or incompatible at main.')
    if not 1 <= len(reviewers) <= 6 or any(
            reviewer.get('type') != 'User'
            or not isinstance(reviewer.get('id'), int)
            or not reviewer.get('login') for reviewer in reviewers):
        conflicts.append('Specify one to six repository reviewers who are GitHub users.')
    elif len(_reviewer_ids(reviewers)) != len(reviewers):
        conflicts.append('Each repository reviewer must be unique.')

    variables = snapshot.get('variables') or {}
    current_app_id = variables.get('RELEASE_APP_ID')
    if current_app_id is not None and requested_app_id is not None:
        if str(current_app_id) != str(requested_app_id):
            conflicts.append('Existing RELEASE_APP_ID does not match the requested App ID.')
    if snapshot.get('secret_present') and current_app_id is None:
        conflicts.append(
            'The release App secret exists without RELEASE_APP_ID; manual repair is required.')
    if snapshot.get('secret_present'):
        notices.append(
            'The existing release App secret cannot be inspected; its value and validity are unknown.')

    if not conflicts:
        requested_reviewer_ids = _reviewer_ids(reviewers)
        protection = {
            'wait_timer': 0,
            'prevent_self_review': False,
            'reviewers': [{'type': 'User', 'id': reviewer['id']}
                          for reviewer in reviewers],
            'deployment_branch_policy': {
                'protected_branches': False,
                'custom_branch_policies': True,
            },
        }
        for name in ENVIRONMENTS:
            current = snapshot.get('environments', {}).get(name)
            if current is None:
                actions.append({
                    'kind': 'environment', 'target': name,
                    'expected': None, 'value': protection,
                })
                actions.append({
                    'kind': 'branch_rule', 'target': name,
                    'expected': [], 'value': {'name': 'main', 'type': 'branch'},
                })
                continue
            environment = current.get('environment') or {}
            if _environment_reviewer_ids(environment) != requested_reviewer_ids:
                conflicts.append(
                    f'Environment {name} does not have the requested reviewer set.')
            if environment.get('deployment_branch_policy') != {
                    'protected_branches': False, 'custom_branch_policies': True}:
                conflicts.append(
                    f'Environment {name} has an incompatible deployment branch policy.')
            branches = current.get('branches')
            if branches != [{'name': 'main', 'type': 'branch'}]:
                conflicts.append(
                    f'Environment {name} must have only the main branch deployment rule.')

        for name in DISABLED_VARIABLES:
            if name not in variables:
                actions.append({
                    'kind': 'variable', 'target': name,
                    'expected': None, 'value': 'false',
                })
            elif variables[name] == 'true':
                notices.append(f'{name} is already true and will be preserved.')

        if not snapshot.get('label_present'):
            actions.append({
                'kind': 'label', 'target': 'release-ready',
                'expected': False,
                'value': {'name': 'release-ready', 'color': '0e8a16'},
            })

    if conflicts:
        actions = []
    return {'actions': actions, 'conflicts': conflicts, 'notices': notices}


def _live_variable(gh, repository, name):
    records = gh.pages(f'/repos/{repository}/actions/variables', 'variables')
    for record in records:
        if record.get('name') == name:
            return record.get('value')
    return None


def _environment_payload_matches(environment, expected):
    normalized = _normalize_environment(environment)
    actual_reviewers = set()
    for entry in normalized.get('reviewers') or []:
        reviewer_id = entry.get('id')
        if reviewer_id is None:
            reviewer_id = (entry.get('reviewer') or {}).get('id')
        if entry.get('type') == 'User' and isinstance(reviewer_id, int):
            actual_reviewers.add(reviewer_id)
    expected_reviewers = {entry['id'] for entry in expected['reviewers']}
    return (
        normalized.get('wait_timer') == expected['wait_timer']
        and normalized.get('prevent_self_review') == expected['prevent_self_review']
        and actual_reviewers == expected_reviewers
        and normalized.get('deployment_branch_policy')
        == expected['deployment_branch_policy']
    )


def _action_is_allowed(action, snapshot):
    kind = action.get('kind')
    target = action.get('target')
    expected = action.get('expected')
    value = action.get('value')
    if kind == 'variable':
        return target in DISABLED_VARIABLES and expected is None and value == 'false'
    if kind == 'label':
        return (target == 'release-ready' and expected is False
                and value == {'name': 'release-ready', 'color': '0e8a16'})
    if kind == 'branch_rule':
        return (target in ENVIRONMENTS and expected == []
                and value == {'name': 'main', 'type': 'branch'})
    if kind == 'environment' and target in ENVIRONMENTS and expected is None:
        reviewers = snapshot.get('reviewers') or []
        required = {
            'wait_timer': 0,
            'prevent_self_review': False,
            'reviewers': [{'type': 'User', 'id': reviewer.get('id')}
                          for reviewer in reviewers],
            'deployment_branch_policy': {
                'protected_branches': False,
                'custom_branch_policies': True,
            },
        }
        return value == required
    return False


def apply_settings(gh, snapshot: dict, plan: dict, completed: list[str]) -> dict:
    """Apply only planned non-credential settings with target-level rechecks."""
    if plan.get('conflicts'):
        raise SetupError('Settings conflicts must be resolved before apply.')
    repository = (snapshot.get('repository') or {}).get('full_name')
    if not isinstance(repository, str):
        raise SetupError('The settings snapshot has no repository identity.')
    if any(not _action_is_allowed(action, snapshot)
           for action in plan.get('actions', [])):
        raise SetupError('The settings plan contains an unsupported action.')

    for action in plan.get('actions', []):
        kind = action['kind']
        target = action['target']
        expected = action['expected']
        value = action['value']
        if kind == 'environment':
            route = f'/repos/{repository}/environments/{target}'
            current = gh.get(route, optional=True)
            if current is not expected:
                raise SetupError('A target changed after planning; no further settings were changed.')
            gh.write('PUT', route, value)
            completed.append(f'environment:{target}')
            created = gh.get(route, optional=True)
            if created is None or not _environment_payload_matches(created, value):
                raise SetupError(
                    'The environment needs explicit manual repair before setup can continue.')
        elif kind == 'branch_rule':
            route = (f'/repos/{repository}/environments/{target}/'
                     'deployment-branch-policies')
            current = _branch_specs(gh.pages(route, 'branch_policies'))
            if current != expected:
                raise SetupError('A target changed after planning; no further settings were changed.')
            gh.write('POST', route, value)
            completed.append(f'branch-rule:{target}')
            current = _branch_specs(gh.pages(route, 'branch_policies'))
            if current != [value]:
                raise SetupError(
                    'The environment needs explicit manual repair before setup can continue.')
        elif kind == 'variable':
            current = _live_variable(gh, repository, target)
            if current != expected:
                raise SetupError('A target changed after planning; no further settings were changed.')
            gh.write('POST', f'/repos/{repository}/actions/variables',
                     {'name': target, 'value': value})
            completed.append(f'variable:{target}')
            if _live_variable(gh, repository, target) != value:
                raise SetupError('GitHub did not retain a repository variable setting.')
        else:
            route = f'/repos/{repository}/labels/{target}'
            current = gh.get(route, optional=True) is not None
            if current != expected:
                raise SetupError('A target changed after planning; no further settings were changed.')
            gh.write('POST', f'/repos/{repository}/labels', value)
            completed.append(f'label:{target}')
            if gh.get(route, optional=True) is None:
                raise SetupError('GitHub did not retain the repository label setting.')

    reviewer_logins = [reviewer['login'] for reviewer in snapshot.get('reviewers', [])]
    return inspect_settings(gh, repository, reviewer_logins)

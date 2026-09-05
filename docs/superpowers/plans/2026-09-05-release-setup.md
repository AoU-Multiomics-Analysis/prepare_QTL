# Release Setup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reusable, dry-run-first command that configures the existing GitHub release system and handles App credentials without exposing them.

**Architecture:** A pure planner consumes a sanitized GitHub settings snapshot. A guarded transport creates only absent settings after fresh reads. A separate App module handles browser approval, in-memory signing, installation verification, and secret upload; the CLI joins these stages without enabling releases.

**Tech Stack:** Python 3.11+, GitHub CLI, Python standard library, cryptography 50.0.1 for in-memory RSA signing, unittest. Use the existing MiniWDL/PyYAML environment for repository regression checks.

**Spec:** `docs/superpowers/specs/2026-09-05-release-setup-design.md`, approved by the user on 2026-09-05.

## Global Constraints

- The default run is read-only.
- Target GitHub.com initially.
- Do not infer reviewer identity or select a user silently.
- A dry run must not read private-key contents, open a browser, or start a local server.
- Create absent enable variables as `false`.
- Do not replace an existing App ID or secret.
- Do not change repository visibility or package visibility.
- Keep the returned private key in memory and send it to `gh secret set` through standard input.
- Do not log callback URLs, codes, keys, tokens, or raw API error bodies from sensitive requests.
- No live GitHub settings test will run without separate user approval.
- No local Docker build is needed. WDL and analysis scripts are unchanged.
- Do not enable releases, publish images, merge PRs, or submit Terra jobs.

## Workspace and files

Use the existing isolated worktree on `codex/release-setup`. It descends from
the PR #48 integration branch. Do not push this work or modify PR #48 without
user approval. Before opening a separate PR, check whether #48 has merged and
select its base accordingly.

Files and responsibilities:

| File | Responsibility |
| --- | --- |
| `ci/release_setup.py` | Settings snapshot, pure plan, guarded GitHub CLI transport, verification |
| `ci/release_app_setup.py` | App manifest callback, in-memory JWT, installation validation |
| `ci/setup_release.py` | CLI input validation, stage order, safe progress and recovery report |
| `ci/release-setup-requirements.txt` | Pinned App signing dependency |
| `tests/test_release_setup.py` | Planner, fake transport, and CLI regression tests |
| `tests/test_release_app_setup.py` | Callback, cryptography, App permissions, and secret redaction tests |
| `docs/image-release.md` | Setup examples, dependencies, partial-failure recovery, remaining approvals |
| `.github/workflows/image-release-plan.yml` | Run setup tests with no GitHub write token or private key |

## Task 1: Safe settings planner and transport

**Files:** Create `ci/release_setup.py` and `tests/test_release_setup.py`.

**Interfaces:**

```python
class SetupError(Exception):
    """A fixed, non-sensitive message suitable for the CLI."""

# JSON-compatible snapshot, with no credential values:
# repository: {full_name, id, owner_type, admin, main_sha}
# workflows_present: bool
# reviewers: [{type: 'User', id: int, login: str}]
# environments: {name: None | {environment: dict, branches: list}}
# variables: {name: str}
# secret_present: bool
# label_present: bool

def inspect_settings(gh, repository: str, reviewers: list[str]) -> dict: ...
def plan_settings(snapshot: dict, requested_app_id: int | None) -> dict: ...
# plan: {actions: list[dict], conflicts: list[str], notices: list[str]}
# Each action: {kind, target, expected, value}; never contains a PEM or token.
def apply_settings(gh, snapshot: dict, plan: dict, completed: list[str]) -> dict: ...
# Returns a freshly inspected snapshot; app credentials are handled in Task 2/3.

class GitHubCLI:
    def get(self, route: str, *, optional: bool = False): ...
    def pages(self, route: str, key: str) -> list: ...
    def write(self, method: str, route: str, body: dict): ...
    def upload_key(self, repository: str, pem: bytes) -> None: ...
```

- [ ] **Step 1: Write failing planner tests using a complete absent-settings fixture.**

```python
def empty_settings():
    return {
        'repository': {'full_name': 'owner/repo', 'id': 42,
                       'owner_type': 'Organization', 'admin': True,
                       'main_sha': 'a' * 40},
        'workflows_present': True,
        'reviewers': [{'type': 'User', 'id': 7, 'login': 'reviewer'}],
        'environments': {'release-publish': None, 'release-commit': None},
        'variables': {}, 'secret_present': False, 'label_present': False,
    }

def test_new_installation_never_enables_releases(self):
    plan = plan_settings(empty_settings(), None)
    values = {a['target']: a['value'] for a in plan['actions']
              if a['kind'] == 'variable'}
    self.assertEqual(values, {'RELEASE_ENABLED': 'false',
                             'RELEASE_AUTO_TRIGGER': 'false',
                             'RELEASE_AUTO_COMMIT': 'false'})
    self.assertEqual(plan['conflicts'], [])
```

Add cases for an existing true variable (no update action), a different App ID
(conflict), a secret without an App ID (conflict), absent administrator access,
missing main workflows, and an already-complete snapshot (no settings writes).
For environments, test no reviewers, additional branch/tag patterns, absent
main branch rule, unknown/custom protections, and preserved wait timers.

- [ ] **Step 2: Run the tests and confirm the new module is missing.**

```bash
python -m unittest discover -s tests -p test_release_setup.py -v
```

- [ ] **Step 3: Implement the pure planner and exact protection policy.**

Require one to six explicit reviewer users with repository read access. Check
the repository's canonical identity and administrator permission. Require
`image-release.yml` and `image-release-dispatch.yml` at `main`, including their
expected variable/environment names, before an apply can proceed. Resolve
main once and read both descriptors at that SHA.

New environments use this payload, with resolved reviewer IDs:

```python
payload = {
    'wait_timer': 0,
    'prevent_self_review': False,
    'reviewers': [{'type': 'User', 'id': 7}],
    'deployment_branch_policy': {
        'protected_branches': False, 'custom_branch_policies': True,
    },
}
branch_rule = {'name': 'main', 'type': 'branch'}
```

The explicit self-review setting permits the named administrator to approve
their own manual trial. Preserve stricter existing self-review controls.
Existing environments must already have the requested reviewer set and an
exact branch-only main restriction. Preserve other rules and timers without
PUT updates. Report incompatible existing settings as conflicts, not actions.
An existing secret cannot be inspected; report that limitation without
claiming its value or validity is known.

- [ ] **Step 4: Implement the transport with fake-process regression tests.**

Use `gh api --hostname github.com --method GET` explicitly for all reads;
body flags must not implicitly turn reads into POSTs. Use `--input -` for
JSON writes. Parse HTTP status separately from the response body. Only a
resource-specific 404 can mean absent after repository authorization has
been established. Treat 401, 403, rate limits, malformed JSON, and unsupported
environment controls as errors. Never include raw subprocess output in a
`SetupError`. Remove inherited `GH_DEBUG` from subprocess environments.

List variables, secrets, and deployment rules with `per_page=100` and page
iteration until exhausted; detect an unexpectedly repeated page. Read public
GHCR manifests from the repositories in `ci/image-stages.yml` when available.
Report inaccessible/private packages as a readiness issue; do not change them
or claim a manifest read proves package-write permission.

Before each write, re-read the target and compare it with `expected`. Create
the environment, then its main branch policy, then verify both before any
secret upload. Append only non-sensitive operation names to `completed`.
Do not delete partially created resources. A partial environment that fails
the policy check requires explicit manual repair rather than a broad PUT.

Secret upload must follow this process contract:

```python
subprocess.run(
    ['gh', 'secret', 'set', 'RELEASE_APP_PRIVATE_KEY',
     '--repo', repository, '--env', 'release-commit'],
    input=pem, capture_output=True, check=False, env=safe_environment,
)
```

Test the command vector has no key and errors cannot echo the key. The
`upload_key` method does not decide whether upload is allowed; Task 3 rechecks
the environment and secret absence immediately before invoking it.

- [ ] **Step 5: Run focused tests and commit the settings component.**

```bash
python -m unittest discover -s tests -p test_release_setup.py -v
git diff --check
git add ci/release_setup.py tests/test_release_setup.py
git commit -m "Add dry-run release settings planner"
```

## Task 2: Guided App setup and credential verification

**Files:** Create `ci/release_app_setup.py`, `ci/release-setup-requirements.txt`,
and `tests/test_release_app_setup.py`.

**Interfaces:** Consume `SetupError` from Task 1. Export:

```python
def app_manifest(repository: str, name: str, redirect_url: str) -> dict: ...
def validate_callback(path: str, query: str, expected_state: str) -> str: ...
# Returns exactly one validated temporary code; never logs it.
def register_app(repository: str, owner_type: str, *, timeout: int = 600) -> dict: ...
# Returns {id: int, slug: str, pem: bytes}; never print or serialize this object.
def app_jwt(app_id: int, pem: bytes, now: int) -> str: ...
def verify_installation(repository: str, app_id: int, pem: bytes) -> dict: ...
# Returns only {id, slug, repository, permissions}; no secret values.
```

- [ ] **Step 1: Add failing callback and manifest tests.**

```python
def test_callback_requires_exact_state_and_path(self):
    self.assertEqual(validate_callback('/callback', 'state=s&code=c', 's'), 'c')
    for path, query in [('/other', 'state=s&code=c'),
                        ('/callback', 'state=wrong&code=c'),
                        ('/callback', 'state=s&code=c&code=d'),
                        ('/callback', 'state=s')]:
        with self.assertRaises(SetupError):
            validate_callback(path, query, 's')

def test_manifest_has_only_required_permissions(self):
    manifest = app_manifest('owner/repo', 'Release helper',
                            'http://127.0.0.1:12345/callback')
    self.assertEqual(manifest['default_permissions'],
                     {'contents': 'write', 'pull_requests': 'read'})
    self.assertEqual(manifest['default_events'], [])
    self.assertFalse(manifest['hook_attributes']['active'])
```

- [ ] **Step 2: Run the tests and confirm the new module is missing.**

```bash
python -m unittest discover -s tests -p test_release_app_setup.py -v
```

- [ ] **Step 3: Implement the loopback manifest flow.**

Bind to `127.0.0.1` on an available port, never all interfaces. Generate state
with `secrets.token_urlsafe(32)`. Serve a no-cache HTML form on a separate
random start path; its explicit POST action is the official organization or
personal GitHub App registration URL. Escape all HTML attributes. The manifest
contains the repository URL, inactive webhook, exact permissions, and callback
URL. The start page requires a valid Host header and has no external assets.

Suppress HTTP access logs. Reject unknown paths, duplicate query parameters,
invalid state, and invalid codes. Set socket/request timeouts so a partial
request cannot hold the server indefinitely. Close the server after a valid
callback or timeout. Exchange the code with GitHub's manifest-conversion API
using an in-memory HTTPS request with a fixed GitHub API host; do not place the
code in a subprocess argument. Reject redirects. Surface only fixed error
messages. Do not persist the PEM, client secret, or webhook secret.

Display the returned public App ID and installation URL, not credentials.
Guide the user to select only the named repository. Organization approval is
not bypassed. The CLI waits for explicit user continuation after installation,
then calls `verify_installation`; it does not poll or register another App.

- [ ] **Step 4: Implement in-memory App authentication and tests.**

Pin `cryptography==50.0.1` in the new requirements file. Use its RSA PEM loader
and PKCS1v15/SHA256 signing, with JWT claims `iat=now-60`, `exp=now+540`, and
`iss=str(app_id)`. Use compact JSON and URL-safe base64 without padding.
Never use the repository administrator token as an App credential.

Use the App JWT for `GET /app` and
`GET /repos/OWNER/REPO/installation`. Verify App ID, owner, required permissions,
and non-suspended installation. This endpoint proves the App can access the
selected repository without minting an installation token. Do not request an
installation token or expand repository access for this verification. Report
an existing broader installation without changing it; reject unexpected extra
write permissions for a newly created App.

Generate a temporary RSA key in test memory. Verify the JWT signature using
its public key and assert claims/expiry. Fake API responses must cover wrong
App identity, owner mismatch, missing contents write, missing PR read,
suspended installation, failed conversion, rejected redirect, and timeout.
Insert a recognizable fake secret into errors and assert it never reaches
stdout, stderr, exception messages, or returned sanitized reports.

- [ ] **Step 5: Run tests and commit the App component.**

```bash
python -m pip install -r ci/release-setup-requirements.txt
python -m unittest discover -s tests -p test_release_app_setup.py -v
git diff --check
git add ci/release_app_setup.py ci/release-setup-requirements.txt tests/test_release_app_setup.py
git commit -m "Add guarded GitHub App setup flow"
```

## Task 3: CLI integration, recovery, and CI tests

**Files:** Create `ci/setup_release.py`; extend both test files; update
`docs/image-release.md` and `.github/workflows/image-release-plan.yml`.

**Interfaces:** Consume all public functions defined in Tasks 1 and 2. Export
`main(argv=None, *, gh=None, app_service=None) -> int` for fake-transport tests.
The default collaborators are constructed only inside `main`.

- [ ] **Step 1: Add failing CLI tests for read-only mode.**

Use a recording fake with the complete settings fixture from Task 1. Its
mutation methods raise `AssertionError`. Patch private-key reading and App
registration to raise `AssertionError`. A dry run with `--private-key-file`
must still succeed without touching that file:

```python
code = main(['--repo', 'owner/repo', '--reviewer', 'reviewer',
             '--app-id', '123', '--private-key-file', '/not/read/in/dry-run'],
            gh=recording_fake, app_service=forbidden_app_service)
self.assertEqual(code, 0)
self.assertEqual(recording_fake.writes, [])
```

Add argument tests for malformed repository names, missing reviewer, duplicate
reviewer, more than six reviewers, partial credential pairs, and create-App
combined with supplied credentials. Reject CLI key-value arguments: the only
key input is a path. Fail conflicts before any mutation, browser, or key read.

- [ ] **Step 2: Run tests and confirm the CLI is missing.**

```bash
python -m unittest discover -s tests -p test_release_setup.py -v
```

- [ ] **Step 3: Implement explicit staged execution.**

```text
parse → inspect → plan → print non-secret plan
  dry run: return with readiness notices
  apply: recheck → create absent protected settings → verify environments
       → create/reuse App credentials → guide installation → verify App
       → recheck secret absence and App ID → upload absent secret
       → create absent App ID variable → read back metadata → report
```

Print fixed `stage=... status=...` progress lines and a non-secret completed
operation list on failure. Catch known transport/credential exceptions without
chained tracebacks containing request data. Do not catch interrupts as success.
Keep the PEM out of snapshots, plans, dataclass reprs, and logging. Check an
existing private-key file is a regular file and not group/world-readable on
POSIX before reading it; report how to fix permissions without changing them.

Recheck both the environment and App ID/secret before upload. If another actor
creates a secret between reads, do not deliberately replace it. GitHub secret
PUT has no atomic create-only operation: document this remaining race and
require one setup operator at a time; do not claim atomic protection.

On an upload failure, preserve completed non-secret resources and explain how
to generate a replacement App key in GitHub and rerun the existing-App path.
If App ID metadata was not written, print the public App ID for recovery.
When an existing secret is reused, do not claim its contents were verified;
App-key verification covers the supplied key, not GitHub's unreadable stored
secret. A successful protected workflow trial is still required.

The exit codes are 0 for a conflict-free plan or verified configuration, 1 for
configuration/transport failures, and argparse's 2 for invalid CLI arguments.
Do not print a ready message when installation verification is unavailable.

- [ ] **Step 4: Test apply, repeat runs, and recovery with fake transports.**

Assert protected environments precede secret upload. Verify the second run
makes no settings or credential writes. Test stale variables/environments,
inaccessible package notices, secret-upload failure, partial environment
creation, missing installation, and preserved enable variables. Verify no
test starts a real browser, creates an App, modifies GitHub, or reads user keys.

- [ ] **Step 5: Add CI tests and document usage.**

In the existing read-only release-plan job, install the requirements and add:

```yaml
- name: Test release setup without GitHub writes
  run: |
    python -m pip install -r ci/release-setup-requirements.txt
    python -m unittest discover -s tests -p test_release_setup.py -v
    python -m unittest discover -s tests -p test_release_app_setup.py -v
```

Document Python/GitHub CLI prerequisites, dry-run and both apply examples,
reviewer selection, organization approval, default-disabled variables,
preservation of existing settings, partial-setup recovery, the unreadable-secret
limit, the secret-create race, and manual trial/activation steps. Do not change
the image build triggers or add a setup job with administrator credentials.

- [ ] **Step 6: Run complete relevant validation and request security review.**

```bash
python -m unittest discover -s tests -p 'test_*.py'
python ci/setup_release.py --help
python scripts/check_wdl_file_scope.py workflows
actionlint .github/workflows/image-release-plan.yml
git diff --check
```

The help command must not contact GitHub. Review the patch for credential
exposure, unwanted setting writes, API permission handling, callback origin/state
validation, and preservation of protections. Fix blocking findings and rerun
tests before asking permission to push/open a PR. Report that live setup,
GitHub App registration, and Terra execution remain untested until separately
approved trials occur.

- [ ] **Step 7: Commit the integrated CLI and documentation.**

```bash
git add ci/setup_release.py tests/test_release_setup.py tests/test_release_app_setup.py docs/image-release.md .github/workflows/image-release-plan.yml
git commit -m "Expose reusable release setup command and tests"
```

## Plan self-review

- Settings preservation, dry-run safety, pagination, and stale state: Task 1.
- Browser approval, callback safety, in-memory keys, and App verification: Task 2.
- Secret ordering, failure recovery, repeat runs, CI, documentation, and review: Task 3.
- No enablement or live configuration occurs during implementation.
- New dependency is isolated from analysis containers. Published version was
  checked on [PyPI](https://pypi.org/project/cryptography/50.0.1/).
- API authority: [App installation lookup](https://docs.github.com/en/rest/apps/apps),
  [environment settings](https://docs.github.com/en/rest/deployments/environments),
  [manifest registration](https://docs.github.com/en/apps/sharing-github-apps/registering-a-github-app-from-a-manifest).

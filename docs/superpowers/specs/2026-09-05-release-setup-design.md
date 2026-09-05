# Release setup tool design

## Goal

Provide a reusable local command that configures the GitHub release system.
The default run is read-only. An explicit apply run creates missing settings
and guides the user through GitHub App approval. It does not enable automatic
releases, publish images, merge PRs, or run Terra workflows.

## Recommended approach

Use a Python CLI and the user's authenticated GitHub CLI. Keep the setup tool
outside GitHub Actions: it needs administrator authority that the release
workflow must not hold. Target GitHub.com initially. Accept an explicit
`--repo OWNER/REPO` so the tool can be copied into a template repository.

Alternatives considered:

- A setup Actions workflow would require a new administrator secret before it
  could install the settings and credentials. This does not simplify setup.
- Terraform would add a separate provider, state, and credential workflow.
  That is unnecessary for these few settings.

## Command interface

The proposed entry point is `ci/setup_release.py`.

```bash
# Read settings and show proposed changes. No GitHub writes or App creation.
python ci/setup_release.py --repo OWNER/REPO --reviewer GITHUB_LOGIN

# Configure missing settings and guide new App registration/installation.
python ci/setup_release.py --repo OWNER/REPO --reviewer GITHUB_LOGIN \
  --create-app --apply

# Reuse an App and read its key from a local file, not a command-line value.
python ci/setup_release.py --repo OWNER/REPO --reviewer GITHUB_LOGIN \
  --app-id APP_ID --private-key-file /absolute/path/app.pem --apply
```

Allow repeated `--reviewer` inputs for up to six users. Do not infer reviewer
identity or select a user silently. Team reviewers can be added later.
`--create-app` and supplied App credentials are mutually exclusive. A dry run
must not read private-key contents, open a browser, or start a local server.

## Preflight and settings

Read the repository identity, administrator permissions, release workflow
presence on main, environment protections, deployment branch rules, variables,
label, and secret names. Treat authentication/permission errors as failures,
not as evidence that a setting is absent. Read all pages for list endpoints.

Use the existing release contract:

- Environments: `release-publish` and `release-commit`.
- Required reviewer(s) supplied by the user; deployment branch limited to
  `main` using a branch rule, not a tag rule.
- Label: `release-ready`.
- Variables: `RELEASE_APP_ID`, `RELEASE_ENABLED`, `RELEASE_AUTO_TRIGGER`, and
  `RELEASE_AUTO_COMMIT`.
- Environment secret: `release-commit/RELEASE_APP_PRIVATE_KEY`.

Create absent enable variables as `false`. Preserve existing values, including
an already-enabled installation, and report them prominently. Setup is not a
disable or activation command.

Reuse existing environments only when they meet the required protections.
Do not overwrite existing reviewers, wait timers, custom rules, branch rules,
or bypass settings. If existing protections need a change, report the exact
conflict and stop before writes. Do not silently remove protections when a
GitHub plan does not support the required controls.

Do not replace an existing App ID or secret. A repeated run reuses them and
reports any mismatch. Secret rotation is outside this tool's initial scope.
Do not change repository visibility or package visibility. Report missing
GHCR access or private packages as separate setup requirements.

## App registration and credentials

Support an existing App or the official GitHub App manifest flow. New App
registration requests Contents write and Pull requests read, with no event
subscriptions or active webhook. Registration belongs to the repository owner.
Organization policy can require an administrator to approve registration or
installation; the tool must not bypass that approval.

For guided creation, use a short-lived loopback-only callback with a random
state token, exact callback path, and timeout. Open the GitHub registration
page only in apply mode. Validate the callback state before exchanging the
temporary registration code. Do not log callback URLs, codes, keys, tokens,
or raw API error bodies from sensitive requests. Test callback validation
without contacting GitHub.

Keep the returned private key in memory and send it to `gh secret set` through
standard input. GitHub CLI encrypts secret values before upload. Do not write
the key into the repository or pass it through process arguments. If the key
upload fails, give recovery instructions without exposing the key. Do not
automatically delete the App or other resources to roll back a partial setup.

Guide the user to install the App on the selected repository only. Verify its
identity, required permissions, and access to the selected repository before
reporting setup ready. Any temporary installation token must be limited to
the selected repository. Do not silently expand an existing installation.

## Apply and recovery

Show a non-secret plan before apply operations. Re-read each target immediately
before changing it; stop on conflicts. Create environment restrictions before
uploading secrets. Verify environment settings after creation and upload the
key only after the commit environment is verified.

GitHub setup is not transactional. On failure, report which non-secret steps
completed and what remains. Do not automatically delete resources. Rerunning
the command must reuse completed steps without replacing credentials.

A ready report distinguishes configured settings from untested behavior. It
lists the manual publishing trial as the next step and states that automatic
release variables remain false on a new installation. Enabling automation
after a successful trial remains a separate, explicit administrator action.

## Files and tests

Separate the read-only planner, GitHub transport, and App registration flow:

- `ci/setup_release.py`: CLI and staged execution.
- `ci/release_setup.py`: settings inspection, plan, and verification.
- `ci/release_app_setup.py`: App registration and credential handling.
- `tests/test_release_setup.py`: pure plan and fake-transport tests.
- `tests/test_release_app_setup.py`: callback and secret handling tests.
- `docs/image-release.md`: setup commands, recovery, and remaining approvals.

Test missing/existing settings, pagination, permission errors, unsupported
protections, conflicting App IDs, repeat runs, stale state, failed secret
upload, callback state mismatch, timeout, and credential redaction. Dry-run
tests must prove that no mutating transport method or credential read occurs.
Run these tests in the existing read-only release-plan GitHub Actions job.

No local Docker build is needed. No live GitHub settings test will run without
separate user approval. WDL and analysis scripts are unchanged. This tool does
not provide evidence that a complete workflow or call caching works on Terra.

## Sources

- [GitHub environment API](https://docs.github.com/en/rest/deployments/environments)
- [Deployment branch policies](https://docs.github.com/en/rest/deployments/branch-policies)
- [GitHub App manifest registration](https://docs.github.com/en/apps/sharing-github-apps/registering-a-github-app-from-a-manifest)
- [GitHub CLI secret upload](https://cli.github.com/manual/gh_secret_set)

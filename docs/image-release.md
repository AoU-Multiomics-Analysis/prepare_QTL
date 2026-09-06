# Stage image releases

Keep the workflow source on `main`. Keep each image default on an immutable
digest. These are separate references: updating a script does not require you
to change a Terra image input by hand.

## Example: change the reference filter

1. Open a PR that changes `scripts/cell_type_specific_expression/filter_cell_type_beds.R`.
2. The release plan selects `cell_downstream`. The release workflow builds the
   cell-type image, or reuses an image with the same source fingerprint.
3. The test job applies the proposed downstream digest to an isolated copy of
   the PR. It keeps the estimation, fit, export, and QTL digests unchanged.
4. The tests run that image combination. If they pass, the optional commit job
   adds only the registered image-default edits to the same PR.
5. Review and merge the PR. Refresh the `main` workflow version in
   Dockstore/Terra as usual. New submissions use the new downstream default.

```text
script PR → source image → candidate WDL pins → tests → pin commit → merge
```

There is no automatic merge. Existing Terra submissions do not change.
Explicit image inputs saved in a Terra configuration override WDL defaults;
remove those overrides if you want to use the maintained defaults.

## Administrator setup

Use `ci/setup_release.py` from a local checkout to configure missing release
settings. The default command reads settings and prints a plan. Only `--apply`
permits changes. The command does not enable releases or start workflows.

Install Python 3.11 or later and the [GitHub CLI](https://cli.github.com/).
Authenticate `gh` to GitHub.com with repository administrator access. Merge the
release infrastructure into `main` before setup. The command checks both release
workflows at one `main` revision. It supports GitHub.com only.

```bash
python -m venv .venv-release-setup
. .venv-release-setup/bin/activate
python -m pip install -r ci/release-setup-requirements.txt
gh auth login --hostname github.com

# Read-only plan. No key file read, browser, local server, or GitHub write.
python ci/setup_release.py --repo OWNER/REPO --reviewer GITHUB_LOGIN

# Create missing settings and guide new App registration and installation.
python ci/setup_release.py --repo OWNER/REPO --reviewer GITHUB_LOGIN \
  --create-app --apply

# Use an existing App. The only private-key input is a file path.
chmod 600 /absolute/path/app.pem
python ci/setup_release.py --repo OWNER/REPO --reviewer GITHUB_LOGIN \
  --app-id APP_ID --private-key-file /absolute/path/app.pem --apply
```

Choose each reviewer explicitly. Repeat `--reviewer LOGIN` for up to six
different GitHub users with repository read access. Teams are not supported.
The tool creates `release-publish` and `release-commit` with those reviewers
and a deployment rule for the `main` branch. A tag rule does not meet this
requirement. Existing reviewer sets must match. Existing wait timers, custom
rules, self-review controls, and bypass settings are preserved. If the existing
protection is incompatible, setup stops before changes. Resolve that conflict
with a repository administrator. The tool does not remove protection when a
GitHub plan or organization policy restricts it.

The tool creates the `release-ready` label and absent release variables. It
preserves every existing variable value, including `true`. It is not a command
to disable an existing release installation. It also preserves an existing
`RELEASE_APP_ID` and `release-commit/RELEASE_APP_PRIVATE_KEY` secret. A different
requested App ID is a conflict. Secret rotation is outside this command.

New Apps request Contents write and Pull requests read, without event
subscriptions or an active webhook. GitHub hosts the registration approval.
The tool shows the public App ID and installation link. Select the target
repository only. Obtain organization approval when required. After installation,
type `continue` in the terminal. The tool then verifies App identity, owner,
permissions, and target repository access. It rejects a new App with broader
write permissions or all-repository access. For an existing App, it reports
broader access and does not change it. Selected-repository mode confirms target
access; it does not prove that no other repositories are selected.

`--create-app` cannot be combined with supplied credentials. It also stops if
the repository already has either App credential setting. Use the existing-App
command for a repeat run. Apply without a credential mode stops before changes
because installation verification needs a key. A dry run needs no credentials.
Even when a dry run includes `--private-key-file`, the file is not read.

The tool checks that a supplied key is a regular file. On POSIX systems, it
rejects group or world permissions and symbolic links. It does not change file
permissions. A generated App key stays in memory and goes to `gh secret set`
through standard input. It is not printed, saved, or placed in process arguments.
The local registration callback binds only to `127.0.0.1` and expires after
ten minutes. It validates the callback path and a random state token.

Only one operator must run setup at a time. The tool rechecks protection, App ID,
and secret names before upload. It stops if another operator changes the
settings. However, GitHub secret PUT has no atomic create-only option. Another
write between the last read and the upload can still be replaced. These checks
do not remove that race.

The tool cannot read GitHub's stored secret or compare it with a supplied key.
When a secret already exists, verification covers the supplied key only. The
stored secret remains unverified. A successful protected workflow trial is
required to test the stored credential.

Allow this repository's build job to publish the registered GHCR packages.
The registry reader requires public GHCR packages. Setup reports an unreadable
public manifest as a separate requirement. That read does not prove package
write access. The tool does not change repository or package visibility.

### Partial setup recovery

Setup is not a transaction. On failure, the tool lists completed operation
names without secret values. It preserves those resources. Exit code 0 means
a conflict-free dry-run plan or verified configuration. Exit code 1 means a
configuration or transport failure; invalid command arguments use exit code 2.
An unavailable installation check cannot produce a completion report.

If environment creation is incomplete, inspect its reviewers and branch rules
in GitHub. Repair the protection explicitly before rerunning setup. Do not
remove protection to make setup pass. Compatible completed resources are reused.

If registration succeeds but installation or key upload fails, keep the App.
Complete its installation and obtain any organization approval. Check the secret
names because an upload error can have an uncertain result. If the secret is
absent, generate a replacement private key in the GitHub App settings. Save the
file privately and use the existing-App command with the displayed public App ID.
Do not rerun new App creation to recover that key. The in-memory key is not saved.

If key upload succeeds but App ID creation fails, the tool prints the public
App ID. Confirm the App identity, then manually create the missing
`RELEASE_APP_ID` repository variable with that value. Preserve the stored secret.
Setup deliberately stops when a secret exists without App ID metadata. Once the
metadata is repaired, a repeat run with existing-App credentials can verify
the supplied key and reuse the stored secret. A protected workflow trial is
still required. No App, environment, variable, or secret is deleted on failure.

Live setup and GitHub App registration have not been tested for this change.
The CI tests use fake GitHub responses and test-only keys. They do not create
an App, start a real browser, or change GitHub settings. The setup command runs
locally; no Actions job receives administrator credentials.

### Manual trial and activation

Use these repository variables:

| Variable | Initial value | Effect when `true` |
| --- | --- | --- |
| `RELEASE_ENABLED` | `false` | Allows manual stage releases |
| `RELEASE_AUTO_TRIGGER` | `false` | Dispatches releases for PRs labeled `release-ready` |
| `RELEASE_AUTO_COMMIT` | `false` | Allows tested pin commits when dispatch also requests them |

The first trial needs only `RELEASE_ENABLED=true` and the publish environment.
Keep the other two variables false. In Actions, select **Stage Image Release**,
run it from `main`, enter a same-repository PR number, and leave `commit_pins`
false. This trial can publish images after approval, but cannot write to the PR.
Review the `release-test-results` artifact and its `candidate/pins.patch`.

Before a publishing trial, **Pinned Image Smoke (Manual)** can test all current
stage defaults. It pulls published images only. It does not build, publish, or
commit. Select the integration branch to test it before merge.

After a successful trial, configure the App and enable `RELEASE_AUTO_COMMIT`.
Run another manual trial with `commit_pins=true`. Confirm that the App commit
starts the normal PR checks. Enable `RELEASE_AUTO_TRIGGER` only after that works.
The automatic dispatcher reads PR metadata only; it does not check out PR code.

GitHub documents the [App token action](https://github.com/actions/create-github-app-token)
and the security constraints for
[`pull_request_target`](https://docs.github.com/en/actions/reference/workflows-and-actions/events-that-trigger-workflows#pull_request_target).

## Trust and stale results

The controller and test policy come from `main`, not from the PR. Releases
accept only open, same-repository PRs targeting current `main`. Update the PR
branch with current `main` before dispatch. A changed base or head stops the
release; dispatch again after review.

The initial policy rejects PRs that also change `ci/`, `.github/`, or `tests/`.
Merge reviewed policy and test changes first, then release the source change
in a separate PR. This is intentional: candidate code cannot replace its own
release tests. The integration PR itself cannot run this new release path
until its trusted policy has been merged.

Only the build job has package-write permission. It builds candidate source
after the publish-environment approval. Test jobs have no write credentials.
The commit job has the App token, but never runs candidate code. It independently
checks the source fingerprint and registry digest, recreates the allowed patch,
and uses an expected branch-head ID to prevent overwriting a newer commit.

The GitHub environments and App installation are part of the security boundary.
Review changes to Dockerfiles before approving publication. Do not make the
publish or commit jobs available to untrusted forks.

## Add scripts and stages

For an existing stage:

1. Put the new script in one of the stage's owned `script_roots` directories in
   `ci/image-stages.yml`.
2. In WDL, pass the matching registered stage image through each call. Use that
   input for the task's `runtime.docker` value.
3. Add or extend a runtime test that runs the script with that stage image.

A change to a shared script invalidates every declared consumer stage. Put a
helper in a shared directory only when each declared consumer image can run it.

A new stage also needs an entry in `ci/image-stages.yml`, registered immutable
defaults in `ci/release-pins.yml`, and a runtime gate in
`ci/test_release_images.py`. Add the stage to `SUPPORTED_STAGES` only after the
runtime gate exists.

The per-file `sources`, legacy runtime aliases, and legacy runtime roots are
transitional. Use owned directories and canonical runtime paths for new code.
This policy-only change does not split existing images, move source files, or
build or publish new stage images. Do not move sources until the trusted policy
and its tests are on `main`. Image splitting and publication need a separate
build and migration plan.

## Build reuse and image retention

### CI events and duplicate work

Source PRs run descriptor checks, release-policy checks, and relevant unit tests.
Cell-type and expression source unit tests pull the pinned dependency images;
they do not build images or run the complete workflow. A dependency change can
need a tested new digest before these source unit tests pass.

The `release-ready` label starts the stage release: build or reuse the selected
image, test the mixed image defaults, and commit the tested pins. The dispatcher
does not repeat a release for a single digest-only commit from
`aou-prepare-qtl-release[bot]`. It checks the parent SHA, author, message, file
statuses, and changed lines. Missing or uncertain metadata uses the normal
release path. Other WDL edits and source edits still dispatch. Adding the label
explicitly can request another trial.

After the pin commit, descriptor and policy checks still run. Source unit tests
compare the previous PR head with the new head and skip an update with no source
changes. After merge, lightweight descriptor and RNA-SeQC checks remain; images
are not rebuilt merely because the PR was merged.

The standard, cell-type, and Rust compatibility image builders are manual-only.
They no longer keep mutable `main` image tags current automatically. Use the
maintained WDL digest defaults. Manual compatibility builds remain available.
RNA-SeQC container CI tests without publishing; stage releases own publication.

Candidate changes to the cell-type or RNA-SeQC runtime tests still run their
container tests on PRs, because the trusted release path rejects candidate
changes to its own tests. These test-maintenance runs do not repeat after merge.
Cell-type test PRs pull published dependency images pinned by digest in the WDL.
The R tests use checked-out scripts mounted into those images. No images are
built in this job. Dependency changes need a tested image from the release process. The HSPE workflow,
precomputed workflow, and saved-model restart run only when smoke scripts,
fixtures, or the reference-fixture generator change, or on manual dispatch.
Full workflow tests use the existing pinned-image runner and all stage defaults.
A change to the CI YAML alone runs the image pulls and R tests; use manual dispatch
to check changes to the integration steps. Stage release tests still run the
complete integration suite.
An unlabelled source PR has not passed the stage release: do not merge runtime
changes until their release tests and final pin checks pass. This trigger change
does not add branch protection or automatic merging.

`ci/image-stages.yml` maps build inputs and stage consumers.
`ci/release-pins.yml` lists the exact WDL defaults that can change. See
[repository image routing](repository-image-routing.md) for input names.

Each release image gets a `release-src-<fingerprint>` tag and a source-fingerprint
label. The fingerprint includes tracked build-file contents, paths, and modes.
WDL pin commits do not change it. A repeat run reuses the matching image and
does not create a commit when its pins are already current. An automatic pin
commit causes another descriptor/policy validation run, but its recognized
digest-only synchronization does not dispatch another release.

Retain these tags and all digests referenced by released WDLs. Do not overwrite
release source tags. Removing an image can prevent a historical workflow from
running. Existing `:main` image build workflows remain available, but moving
that tag does not change a digest default.

Build-input patterns must include every file a Dockerfile copies. Add a new
stage, its pin targets, and a runtime test gate before adding a new workflow
family. External genotype images are not built or released by this system.

## Tests and limits

The release gate checks all WDL descriptors and the Terra file-creation rule.
It also checks candidate image forwarding, not just declared defaults.

| Family | Runtime gate |
| --- | --- |
| Cell type / expression / common | HSPE and supplied-proportion workflows, filtered QTL outputs, and saved-model restart with the actual stage defaults; expression scale and sample-list checks |
| Methylation R | Synthetic normalization, filtering, and output checks |
| Methylation Rust | Both bundled binaries on a small cohort, including known output values |
| Splicing | Bundled transformations and a synthetic BED-to-output run |
| Proteomics | Bundled helper functions and CLI/dependency startup |
| RNA-SeQC | Bundled script tests and candidate WDL task checks |

Proteomics BioMart mapping is not covered by an offline end-to-end test. The
workflow needs that external service. These tests do not prove scientific
accuracy on production data.

The complete workflows have **not been tested on Terra for this change**.
GitHub tests do not prove Terra call-cache reuse. A separate, approved Terra
trial must confirm cache behavior. Stable image digests remove one source of
cache changes; inputs, task commands, and other Cromwell settings can still
prevent a cache hit.

To roll back, disable automatic release variables and restore the required
older digest defaults in a reviewed PR. Do not delete the newer image or
rewrite the release branch history.

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

The new release workflow is disabled until `RELEASE_ENABLED` is `true`.
This PR does not create secrets, enable settings, or start a release.

Before you enable it:

1. Merge the release infrastructure into `main`.
2. Create GitHub environments named `release-publish` and `release-commit`.
   Configure required reviewers and restrict deployments to `main`. GitHub
   does not add protection rules merely because a workflow names an environment.
3. Allow this repository's build job to publish the four registered GHCR
   packages. The registry reader currently requires **public GHCR packages**.
4. For automatic pin commits, install a GitHub App on this repository only.
   Give it **Contents: read and write** and **Pull requests: read** permissions.
   It does not need workflow-write or administrator permissions.
5. Set repository variable `RELEASE_APP_ID` to the App ID. Store the private
   key as secret `RELEASE_APP_PRIVATE_KEY` in the `release-commit` environment.

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

## Build reuse and image retention

`ci/image-stages.yml` maps build inputs and stage consumers.
`ci/release-pins.yml` lists the exact WDL defaults that can change. See
[repository image routing](repository-image-routing.md) for input names.

Each release image gets a `release-src-<fingerprint>` tag and a source-fingerprint
label. The fingerprint includes tracked build-file contents, paths, and modes.
WDL pin commits do not change it. A repeat run reuses the matching image and
does not create a commit when its pins are already current. An automatic pin
commit can cause one more validation run, but it does not cause another image
build or another pin commit.

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

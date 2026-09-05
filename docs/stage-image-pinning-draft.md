# Draft: automatic stage-specific image pinning

Status: design only. No workflow, registry, repository permission, or Terra
configuration is changed by this document.

## Purpose

Keep using the WDL from `main` without manually copying image digests. A code
change should update the image used by the affected stages, not every stage.
This preserves one important condition for Cromwell call caching. It does not
guarantee a cache hit: commands, task inputs, cache settings, and output access
must also match.

Keep the existing Dockerfiles. Different stages may use different published
versions of the same image.

## Proposed workflow inputs

Both cell-type entry points expose the following String inputs. Each default
contains a verified `ghcr.io/...@sha256:...` reference, not a mutable tag.

| Input | Assigned work |
| --- | --- |
| `estimation_docker_image` | Input gene filtering, HSPE preparation/workers/merge, proportion validation and processing |
| `fit_docker_image` | TCA fitting and numerical model cleanup |
| `export_docker_image` | BED export and export QC |
| `downstream_docker_image` | Gene summaries, Haemopedia preparation, reference filtering, scatter preparation, and manifests |
| `qtl_docker_image` | The existing expression preparation workflow and its QTL tasks |

Each task receives only its assigned image through its existing `docker_image`
or `DockerImage` input. Do not pass the complete image map into every task;
that would change unrelated task inputs when one stage is updated.

Initially, the first four defaults may share one verified digest. QTL uses a
verified digest of the existing standard image. Select these digests from
published, tested images; do not invent them or resolve mutable tags at runtime.

Retire the shared `deconvolution_docker_image` input in the cell-type entry
points. Update the examples and smoke inputs at the same time. Users with a
saved override must replace it with the relevant stage inputs. An optional
legacy override is possible, but is not part of this draft: silently applying
one image to every stage would hide the new behavior.

## Identify affected stages

Add a reviewed YAML dependency map under `ci/`. This is CI configuration, not
a JSON wrapper passed from WDL to R.

- Assign each CLI script and R helper module to its consumers.
- Include shared helpers transitively, not just the script named by a task.
- Treat bootstrap changes and image dependency changes conservatively: update
  all stages that use that image.
- Fail the release check for a new or renamed source file without a mapping.
  Require review rather than silently treating it as unrelated.
- Test-only edits run relevant tests but do not update image defaults.
- WDL-only edits run validation and release-consistency checks, not automatic
  image builds. A changed CLI contract must reference an image that supports it;
  it cannot be approved solely because the WDL parses.

Example: a change confined to `filter_cell_type_beds.R` or its filtering-only
helpers updates the downstream pin. A change to a shared I/O helper may update
several pins. The current bootstrap loads all R modules, so this dependency
audit must be conservative and backed by tests.

The standard Dockerfile copies all of `scripts/`. It may still rebuild for
changes outside QTL preparation. Such a build must not automatically advance
the QTL pin if no QTL code or dependency changed.

## Release sequence for one code PR

1. Record the PR head commit and current base commit. Determine affected stages
   against the base. Serialize release operations per PR.
2. Build the candidate image from that exact source commit in GitHub Actions.
   Publish it under a retained build-specific tag and record its immutable digest.
   Do not make a new build merely to obtain a digest after testing.
3. Create a proposed pin update: new digest for affected stages; existing base
   digests for unaffected stages.
4. Test the actual mixed-image workflow, pulling images by digest. Do not override
   every stage with one newly built test image. Check stage interfaces, model/BED
   restart paths, downstream outputs, and WDL 1.0 compatibility.
5. A separate trusted updater checks the source commit and test result. If the PR
   head or relevant base has changed, discard the stale result and reevaluate.
6. The updater changes only the approved image-default fields in both WDL entry
   points and adds a commit to the same PR. It does not force-push.
7. Validate the final PR head. Require matching defaults in both entry points,
   immutable image references, and a release record tied to the tested source and
   digests. The final tree may differ from the tested tree only in the approved
   generated pin fields. Any other change needs fresh evaluation.
8. The user reviews and merges the PR. `main` then contains the compatible WDL
   and image references together.

## Avoid repeated builds

The updater must not build on its own pin-only commit. Use the dependency map
and a source fingerprint that excludes only the generated pin fields, not all
WDL content. Include relevant build inputs and base pins in the release record.
Reapplying an already-current pin must produce no commit.

Keep the existing WDL-only build exclusions. Add a lightweight release check
that runs on pin-only commits. Do not rely only on a bot commit message to skip
work: a human commit could use the same message.

If another PR changes relevant source or pins on `main`, reevaluate the candidate
against that base. Do not overwrite newer pins with an old build result.

## Permissions and repository setup

Use a repository-scoped GitHub App token for the trusted updater, with only the
permissions needed to commit approved pin changes and report PR status. Keep
that token out of build/test jobs that execute candidate code. Do not run PR
scripts or Dockerfiles in a privileged `pull_request_target` job.

Start with maintainer-approved, same-repository PRs. Fork PRs may run read-only
checks, but do not automatically receive publish or updater credentials. The
trusted updater treats build results as data: validate the PR, commit, repository,
digest format, and allowed file changes before writing anything.

A GitHub App also allows its pin commit to trigger the next validation run.
Ordinary `GITHUB_TOKEN` pushes generally do not start another workflow run;
explicit dispatch is an alternative if an App is not desired.

Before enforcement, a repository administrator must approve App installation,
required checks, and image retention. These settings are not changed by drafting
or implementing the WDL code alone.

References: [GitHub workflow triggering](https://docs.github.com/en/actions/how-tos/write-workflows/choose-when-workflows-run/trigger-a-workflow),
[privileged workflow safety](https://docs.github.com/en/actions/reference/security/securely-using-pull_request_target),
[Docker build digest output](https://github.com/docker/build-push-action/).

## What the user does

Push code, wait for build/test/pin checks, review the automatic commit, and merge.
There is no manual digest copying and no separate digest-update PR to track.

Keep the same WDL source URL on `main`. In Terra, select or refresh the workflow
version/configuration as needed to load new defaults. Existing submission
snapshots and explicit saved image overrides do not automatically change just
because the GitHub branch has advanced.

For an exact repeat, retain the WDL revision and resolved image inputs from the
original submission. Using the newest WDL on `main` is not an exact-repeat policy.

## Validation and rollout

- Unit-test dependency selection, default updates, unmapped files, stale builds,
  no-op updates, and consistent defaults across entry points.
- Validate WDL 1.0 with MiniWDL and Womtool, plus existing Terra file-path checks.
- Run mixed-image smoke tests in GitHub Actions. No local Docker builds.
- Verify a filtering-only update leaves estimation, fitting, and export image
  inputs unchanged. Verify shared code updates all affected stages.
- After approval, compare two otherwise identical Terra runs using fixed inputs
  and pins. Inspect cache results and metadata. Never claim Terra caching success
  based only on local or container tests.

The initial switch from a String `:main` task input to a digest may cause a
one-time cache miss. This design is primarily for stable future runs. It does
not recover a failed task midway through execution or replace an explicit model
or BED checkpoint.

Open setup choices: approve the GitHub App approach; choose initial verified
digests; audit the stage dependency map; and select image-retention and required
check policies. No secrets belong in this document.

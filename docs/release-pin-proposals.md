# Candidate image pin patches

This is the next implementation layer after the report-only planner and stage
image inputs. It creates a patch and a release record. It does not build images,
publish packages, edit working WDLs, commit, push, or approve a release.

## Inputs and output

`ci/release-pins.yml` maps stage names to exact workflow input defaults. The
current targets cover the four cell-type stages and the integrated QTL image.
Other workflow families still need explicit targets before they can use this
tool. A source change affecting an unsupported stage is rejected, not ignored.

The tool reads committed configuration and WDLs from the supplied head commit.
It uses the dependency planner to select stages. A filtering-only script change
can request a cell image build and a standard image build, but only the
downstream defaults change when only that stage consumes the changed code.

Example invocation after a build has produced a candidate digest:

```bash
python ci/propose_image_pins.py \
  --base origin/main \
  --head HEAD \
  --expected-base BASE_COMMIT_SHA \
  --expected-head BUILD_SOURCE_COMMIT_SHA \
  --candidate-image cell_type=ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression@sha256:BUILD_DIGEST \
  --output-dir release-candidate
```

Replace the capitalized values with actual full commit IDs and the 64-character
digest. `BUILD_DIGEST` alone is not a valid digest. The output directory must
not already exist. There is no overwrite option.

The directory contains:

- `pins.patch`: only changed workflow String defaults, suitable for review with
  `git apply --check release-candidate/pins.patch` in the matching source checkout.
- `release-candidate.json`: base/head commits, dependency plan, supplied image
  references, changed files, and the SHA-256 of the patch. Its status is always
  `candidate_not_validated`.

This JSON is a CI release artifact, not a WDL-to-script argument wrapper.

The tool locates inputs with the MiniWDL parser. It does not replace image strings
globally. Defaults must already be immutable image references. Duplicate targets,
missing inputs, nonliteral defaults, unexpected repositories, inconsistent entry
point defaults, and missing candidate images cause a failure before output.
An already-current digest produces an empty patch.

## What the guards prove

The selected local base/head refs must match the expected commit IDs, and base
must be an ancestor of head. The tool checks them again before emitting files.
It ignores uncommitted configuration and WDL edits. This protects candidate
generation against local ref changes; it does not prove a remote PR is current.

A supplied digest is checked for syntax and the expected repository only. This
tool does not contact the registry or certify that an image was built from the
recorded head. The release record and patch are not trusted attestations.

## Remaining integration

The complete automated release still needs:

1. A build job tied to the exact source commit that publishes and records the
   tested image digest, rather than accepting an arbitrary user assertion.
2. A test job that applies the candidate patch in an isolated checkout and runs
   the resulting mixture of old and new images, including restart paths.
3. A separate trusted updater that verifies the live PR head/base, image source,
   successful test run, patch contents, and allowed target fields before it commits.
4. A source fingerprint and pin-only commit policy to prevent repeated builds.
5. Approved repository permissions and image retention. The updater must not
   expose its credentials to PR scripts, Dockerfiles, or untrusted artifacts.

No GitHub App or write-enabled updater is installed by this change. Existing
smoke failures still block rollout. No Terra jobs or local container builds are
performed. The WDL source can remain on `main` once the compatible code and pins
are merged together; this candidate tool alone does not automate that merge.

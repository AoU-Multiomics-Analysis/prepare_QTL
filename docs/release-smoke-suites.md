# Release smoke suites

The image-release gate runs short tests for changed stages against immutable
image digests. Integration is manual. See [fast image releases](lean-release-tests.md)
for the default policy. The optional integration suites are described below.

| Check | Compact release suite | Full suite |
| --- | --- | --- |
| HSPE and multi-batch estimation | Deconvolution workflow, all ten groups; no QTL scatter | Complete QTL workflow, all ten groups |
| Precomputed proportions | Complete QTL workflow, three groups | Complete QTL workflow, ten groups |
| Saved-model restart | Deconvolution only; verify fit skips and identical BED contents | Complete QTL workflow again |
| Other affected images | Retain their targeted tests | Retain their targeted tests |

The compact precomputed case retains B cells, CD4 T cells, and monocyte/myeloid
using the existing synthetic fixture and a test-only mean-proportion threshold
of 0.12. It retains all 12 samples and all existing genes. An independent R test
checks the retained groups and weights. The manifest, sample order, coordinates,
INT/scaled matrices, PCs, covariates, and reference-filter assertions remain.
HSPE retains its original threshold and all reference cell types.

Compact CPU inputs are one core. Memory requests are 4 GB rather than 8 GB
for the small test data. These are fixture overrides, not Terra defaults or
claims about memory requirements for real cohorts. GitHub integration checks
must confirm that the new reservations suffice before this policy is merged.
Expression/common checks run once per distinct image digest, not twice when
both stages use the same image.

## Run on GitHub

The **Pinned Image Smoke** workflow is manual. It offers a `suite` choice and
defaults to `full`. The separate cell-type integration workflow is also manual.

On a GitHub runner with Docker available:

```sh
python ci/test_release_images.py --source . --all-stages --suite compact
python ci/test_release_images.py --source . --all-stages --suite full
```

Neither command builds images. Do not run local Docker builds for this change.
The compact suite has no Terra execution; a successful GitHub run is not a
complete Terra validation. No fixed speedup is promised until a GitHub run is timed.

These are test-policy changes. Merge them separately from source/image releases
and do not apply `release-ready` to their PR.

## Stage-specific release tests

`ci/image-stages.yml` lists `runtime_tests` for each cell-type stage. Each
selected stage runs those tests in its pinned image. The test harness links the
scripts bundled inside that image into a temporary test checkout. It does not
replace the image scripts with scripts from the mounted candidate checkout.
Missing HSPE, TCA, or edgeR dependencies fail this gate rather than skip it.

| Changed stage | Task-level coverage |
| --- | --- |
| Estimation | HSPE batching |
| Fit | Small TCA fit, model cleanup, and export gene/sample alignment |
| Export | Tensor extraction, BED coordinates, sample order, and QC |
| Downstream | Reference filtering |

These tests use small synthetic fixtures. For example, the export test builds
a tiny model to check extraction; it does not launch the complete FitTca WDL.
Other image families retain their existing targeted checks. Cheap repository
and WDL validation still runs for every release.

New scripts, shared helpers, and dependency changes use the affected-stage tests
without an automatic integration run. New stages still need a registered test.
For input/output or sample/gene-order changes, request manual integration when
needed. WDL-only changes retain static checks without runtime image tests.

The release log and `ci-runs/runtime-test-plan.json` show the selected stages
and test files. `--suite compact`, `--suite full`, or `--all-stages` explicitly
requests integration. Broad source tests are also manual.

## Reuse images within one job

The release runner and pinned smoke runner inspect each exact digest before
pulling it. If that digest is already present in the job's Docker daemon, the
runner logs `status=reused` and skips the pull. A different or missing digest
is pulled normally; a failed pull still fails the test job. Mutable tags are
not accepted. This does not retain images across fresh GitHub runners or
change MiniWDL's own Docker behavior.

# Release smoke suites

The image-release gate selects tests for changed stages. When an integration
check is required, it uses `compact`. It still tests the selected immutable
image digests before the release App can commit them. No production WDL,
analysis script, input default, or output format changes for these tests.

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

The **Pinned Image Smoke** workflow tests `compact` on pull requests that change
its release-test code. Manual runs offer a `suite` choice and default to `full`.
The separate cell-type integration workflow keeps its full-suite coverage.

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
| Estimation | Expression input, markers, HSPE, batching, and proportions |
| Fit | Small TCA fit, model cleanup, and model reuse |
| Export | Tensor extraction, BED coordinates, sample order, and QC |
| Downstream | Reference filtering, summaries, scatter inputs, and manifests |

These tests use small synthetic fixtures. For example, the export test builds
a tiny model to check extraction; it does not launch the complete FitTca WDL.
Other image families retain their existing targeted checks. Cheap repository
and WDL validation still runs for every release.

Only exact files in `stage_only_test_paths` can omit integration. The initial
list covers the fit entrypoint, HSPE estimator module, QC module, and gene
summary implementation. All other relevant code changes retain the compact
test. In particular, BED interfaces, shared helpers, dependency environments,
new scripts, and WDL changes keep that check. Missing change evidence also
keeps integration. The runner reads the exact PR base/head diff from the
trusted checkout, not from a candidate-provided test-selection file.

Path rules cannot detect a semantic interface change. If an exempt file changes
its inputs, outputs, units, or sample/gene ordering, remove its exemption in a
separate test-policy PR before releasing that change. Do not add directory-wide
exemptions for new scripts. New stages must have a runtime gate before release.

The release log and `ci-runs/runtime-test-plan.json` show the selected stages,
test files, and paths that require integration. `--suite compact` or
`--suite full` explicitly requests integration; `--all-stages` also keeps it.
Manual full tests and the separate source-unit workflow are unchanged.

For this test-policy PR, Pinned Image Smoke runs all stage gates against the
published images plus compact integration. This is broader than a normal
single-stage release so the new test harness is checked before use.

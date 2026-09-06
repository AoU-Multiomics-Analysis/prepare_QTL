# Release smoke suites

The image-release gate uses `compact`. It still tests the selected immutable
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

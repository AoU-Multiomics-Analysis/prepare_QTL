# Initial stage image defaults

Status: candidate defaults, not a completed runtime rollout. Automatic publishing
and pin-update commits are not enabled. Other repository entry points retain
their current image inputs; the first rollout covers the two cell-type entry points.

## Selected images

Verified on 2026-09-05 against public GHCR manifest bytes (SHA-256), manifest
headers, and successful GitHub build logs. Both images target Linux amd64 and
were built from repository commit `f40568d08d687b500ee9ae447c9869db092d4792`.

- Cell-type stages: `ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression@sha256:9f7af7c16fa3dc7a0b82c042a40145fa26afce4a96547791e0b29a9e8de4d754`.
  [Build log](https://github.com/AoU-Multiomics-Analysis/prepare_QTL/actions/runs/33983297365).
- QTL preparation: `ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:932f67a09f1635c22a8061a5c98c892393d321e7c17d0401531e7093c469c845`.
  [Build log](https://github.com/AoU-Multiomics-Analysis/prepare_QTL/actions/runs/33983297326).

The image source labels are inherited from micromamba and do not identify this
repository's source commit. The build logs provide that association instead.

## Input routing

- `estimation_docker_image`: gene filtering, HSPE batches, proportion processing.
- `fit_docker_image`: TCA fitting and model cleanup.
- `export_docker_image`: BED export and export QC.
- `downstream_docker_image`: summaries, reference preparation/filtering, scatter
  preparation, and manifests.
- `qtl_docker_image`: expression preparation and its QTL calls in the parent.

Replace saved `deconvolution_docker_image` overrides with the relevant new inputs.
The default values for estimation, fit, export, and downstream initially match.
They can now change independently. Updating a downstream default alone leaves
the other task image inputs unchanged. Shared source or environment changes can
still require several pins to change; use the dependency planner before release.

Each workflow returns a `stage_images` map. It is a workflow output only, not a
map passed into every task. The legacy manifest field `container_image` now means
the manifest task's downstream image. It must not be used as whole-run provenance.

## Validation and remaining gates

Run `python tests/test_stage_image_inputs.py` to evaluate the actual WDL call
expressions with distinct stage values. Dependency tests remain in
`tests/test_image_update_plan.py`.

The manual **Pinned Image Smoke (Manual)** GitHub Actions workflow pulls the
digest defaults, runs the synthetic HSPE and precomputed-proportion workflows,
and checks their outputs. It never builds or publishes images. It supports
different digest defaults for the different stages. Existing build smoke tests
still use explicit local test tags for every stage.

The published main build's end-to-end smoke run did not pass its QTL filtering
assertion (coordinates/sample order/CPM preservation). Its deconvolution assertions
passed. Do not remove that check to approve these pins. The genotype descriptor
also has a separate external import returning HTTP 404. Neither issue is fixed
by digest pinning. Resolve these validation failures before production rollout.

The pinned-image smoke has not yet been run for this branch. Different historical
cell-type digest combinations and model restart need runtime validation before
automated stage-selective releases. No complete Terra run or cache-hit comparison
has been performed. No Terra jobs were submitted.

The WDL source URL can remain on `main`; refresh the Terra workflow configuration
to load new defaults. Saved overrides and existing submissions do not change
automatically. The initial tag-to-digest input change can cause a one-time cache
miss. Stable images alone do not guarantee cache reuse.

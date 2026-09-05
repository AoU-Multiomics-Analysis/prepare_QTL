# Stage image defaults implementation plan

**Goal:** Give the cell-type workflow independent, immutable image defaults for estimation, fitting, export, downstream work, and QTL preparation.

**Architecture:** Keep the existing task image inputs and Dockerfiles. Each call receives only its stage image. Both entry points expose matching defaults; workflow outputs record the selected images without adding them to unrelated task inputs.

**Spec:** `docs/stage-image-pinning-draft.md`, limited to initial input routing and pinned defaults. Automated publishing and pin commits remain disabled.

**Constraints:** WDL 1.0; no workflow-scope file writing; no local container builds; no Terra jobs. Do not change the separate genotype workflow.

## Tasks

- [ ] Verify GHCR manifest digests and image source labels. Record the exact references and build evidence in documentation.
- [ ] Add `tests/test_stage_image_inputs.py`. Load real WDL ASTs with MiniWDL, evaluate image call-input expressions with distinct stage values, and assert that changing downstream leaves estimation/fit/export calls unchanged. Check immutable, matching defaults and parent forwarding.
- [ ] Run `/opt/homebrew/opt/python@3.11/bin/python3.11 tests/test_stage_image_inputs.py` and confirm failure before edits.
- [ ] Replace the shared input in both cell-type entry points with `estimation_docker_image`, `fit_docker_image`, `export_docker_image`, and `downstream_docker_image`. Pin `qtl_docker_image` in the parent. Add a `stage_images` workflow output map; retain the legacy manifest image field as the manifest-task image only.
- [ ] Update fixture inputs, smoke assertions, and user documentation for the retired input. Keep test image overrides explicit for all stages in existing build smoke tests.
- [ ] Add a manually dispatched, pull-only GitHub Actions smoke job using the actual pinned defaults. Reuse the synthetic fixture and assertions; do not build images. Keep regular WDL-only changes excluded from expensive automatic builds.
- [ ] Run image routing, dependency selection, WDL and Terra file-scope checks, relevant existing tests, and actionlint. Report remote smoke and Terra validation separately from local checks.

## Review gate

Do not claim runtime image compatibility or cache reuse from static tests. Published digest existence and successful builds are not a substitute for the pinned-image smoke run. Leave changes for review; no automatic updater or repository permission changes.

## Implementation checkpoint

The input routing, candidate defaults, workflow image outputs, fixture migration,
documentation, and manual pull-only smoke job are implemented locally.

Verification on 2026-09-05:

- 18 image/planner/CI-selection Python tests passed.
- 18 cell-type WDL regression tests passed; the optional Cromwell test was then
  run separately with Womtool 87 and Java 20 and passed for both entry points.
- The focused R end-to-end WDL and fixture-validation tests passed.
- Terra file-scope checks passed for all 32 WDL files. Actionlint and diff checks passed.
- The older documentation-contract suite fails on both this tree and an exported
  unchanged HEAD baseline. The broader logging check also reports existing gaps
  in other workflow families. These are not repaired by this image-routing change.
- Remote pinned smoke and full Terra execution remain untested. Existing main
  QTL smoke and genotype import failures remain rollout gates.
- Code review found a flat-output JSON handling error in the new smoke runner;
  it was corrected and covered by a failing-then-passing regression test.

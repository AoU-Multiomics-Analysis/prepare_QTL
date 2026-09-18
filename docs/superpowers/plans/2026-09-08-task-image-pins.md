# Task image pins implementation plan

Goal: pin each repository-managed WDL task independently, updating its digest only when its entrypoint, transitive shared dependencies, or image environment changes.

Design: retain shared Docker builds and the existing trusted release controller. Replace coarse release units with named task units. Register literal image defaults on each task and each calling workflow, forwarding only that task's image through nested calls. Preserve current digest values at migration. External images remain explicitly manual because this repository does not build them. Retain WDL 1.0 and File localization.

- [x] Inventory every task, entrypoint and transitive helper dependency; reject missing mappings.
- [x] Migrate pins and forwarding without changing commands or data inputs.
- [x] Make release planning, image validation and runtime selection work with task units.
- [x] Test isolated script edits, shared helper edits, environment edits, unchanged tasks, and complete WDL routing.
- [x] Run syntax/static checks and relevant local tests. Use GitHub for container checks; no Terra submission.
- [x] Open a separate policy migration PR. Do not label it release-ready: trusted policy changes must merge before a source release.

## Current-main reconciliation, 2026-09-18

- [x] Merge current main into #76 without rewriting branch history.
- [x] Carry forward #75's constrained fit, matching variance screening, cleanup,
  and published images. Compare every migrated default with the corresponding
  current-main stage default before running tests.
- [x] Retain #78's command-time inventory localization and regression tests.
- [ ] Run the release-policy Python suite, WDL checks, file-path tests, and
  focused R tests. Run GitHub validation and source checks before merging #76.
- [ ] Update #77 with the merged policy, target main, and verify scoped loading.
  Apply `release-ready` only to #77. Verify published image pins and validation
  before merging it. Do not submit Terra jobs or build local Docker images.

The migration cannot promise cache hits: inputs, task command/output definitions, and cached-output access still matter. Existing historical images are retained rather than claiming a reconstructed last-change digest. Subsequent releases maintain that invariant prospectively.

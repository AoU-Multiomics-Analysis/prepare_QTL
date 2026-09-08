# Task image pins implementation plan

Goal: pin each repository-managed WDL task independently, updating its digest only when its entrypoint, transitive shared dependencies, or image environment changes.

Design: retain shared Docker builds and the existing trusted release controller. Replace coarse release units with named task units. Register literal image defaults on each task and each calling workflow, forwarding only that task's image through nested calls. Preserve current digest values at migration. External images remain explicitly manual because this repository does not build them. Retain WDL 1.0 and File localization.

- [x] Inventory every task, entrypoint and transitive helper dependency; reject missing mappings.
- [x] Migrate pins and forwarding without changing commands or data inputs.
- [x] Make release planning, image validation and runtime selection work with task units.
- [x] Test isolated script edits, shared helper edits, environment edits, unchanged tasks, and complete WDL routing.
- [x] Run syntax/static checks and relevant local tests. Use GitHub for container checks; no Terra submission.
- [ ] Open a separate policy migration PR. Do not label it release-ready: trusted policy changes must merge before a source release.

The migration cannot promise cache hits: inputs, task command/output definitions, and cached-output access still matter. Existing historical images are retained rather than claiming a reconstructed last-change digest. Subsequent releases maintain that invariant prospectively.

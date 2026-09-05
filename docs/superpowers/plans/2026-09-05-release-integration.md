# Release integration implementation plan

**Goal:** Integrate source-specific image builds, stage pin proposals, mixed-image tests, and guarded pin commits in one PR on merged main.

**Architecture:** A trusted workflow runs from main, accepts a same-repository PR number, snapshots its base/head, builds only affected registered images, and produces a candidate patch. Test jobs have no write credentials. The optional final updater independently reconstructs the allowed patch using trusted policy and verifies the live PR before a non-forced update. Separate repository settings enable dispatch and commits; no settings are changed by implementation.

**Spec:** `docs/stage-image-pinning-draft.md` plus user approval for one integration PR with repository-wide image routing and validation fixes.

**Constraints:** WDL 1.0; no workflow-scope file writers; no local Docker builds; no Terra jobs; no automatic merge; no secrets available to PR test code. External genotype images remain manually maintained, not rebuilt by this repository.

## Task 1: Repository-wide image routing

Own WDL files outside cell-type, `ci/release-pins.yml`, `tests/test_repo_image_routing.py`, and `docs/repository-image-routing.md`. Extend image default targets to all consumers of the four repository images, including task-only descriptors and shared workflows. Preserve existing input names where possible, route each call to its stage image, and use verified published digests only. Coordinate any missing published image with the root agent. Do not edit the proposal parser; report the minimal interface extension it needs for task-scoped targets.

Test with MiniWDL ASTs that every repository-owned runtime image is a declared input and registered defaults are immutable. Check actual call forwarding, not merely text. Pin targets should identify workflow/task scope explicitly when necessary. Do not invent digests. Do not change external genotype image sources.

## Task 2: Release engine and workflows

Root owns `ci/release_integration.py`, new GitHub release workflows, `ci/propose_image_pins.py`, `ci/image-stages.yml`, and release engine tests. Implement pure helpers first with failing tests: source fingerprint, exact target patch, ref freshness, no-op behavior, build selection, trusted-policy checks. Record image digest/source fingerprints as build outputs. Use immutable digest refs for tests. Update proposal parser to accept named task scopes from Task 1 while preserving workflow-only compatibility.

A manually dispatched trusted workflow must run only from main and only for open same-repository PRs. Build jobs get package permission only after a protected environment gate; tests get read-only credentials. The optional updater uses a repository GitHub App token only in its own protected job. It must not execute candidate code. Recompute proposed target changes and verify current PR head/base before writing. Never force-push. Source-fingerprint image tags allow unchanged builds to be reused, and an already-current pin produces no commit.

Automatic event dispatch, if implemented, must be metadata-only and disabled by a repository variable. Never check out or execute PR source in pull_request_target. Retain manual dispatch regardless of automatic mode.

## Task 3: Runtime validation and existing failures

Inspect the existing QTL filtering assertion and genotype HTTP404 import failure before changes. Fix the actual cause with regression tests; do not weaken assertions to get a green run. Extend pinned smoke to the other registered image families using available synthetic fixtures; report missing test coverage explicitly. Cover model restart and actual per-stage digest routing. Root initially owns this task and may delegate a bounded diagnostic or test implementation after Task 1 implementation finishes.

## Task 4: Review and rollout documentation

Document GitHub App permissions/secrets, protected environments, enable variables, image retention, rollback, manual trial, update-loop prevention, stale PR handling, and the difference between CI validation and Terra cache validation. Run focused tests, repository descriptor checks, actionlint, and a whole-branch security review. Do not enable settings, publish images, or start Terra runs without approval.

## Progress and rulings

- Base: `4810de5` (PRs 45–47 merged). Existing linked worktree reused on `codex/release-integration`.
- Ruling: manual trusted dispatch is the first release path; automatic event dispatch and App commits remain disabled unless configured. This avoids adding write automation before administrator approval.
- Ruling: external genotype image builds are out of repository scope; their descriptors still receive validation.
- Implemented all four tasks. The reviewer found two candidate-test routing gaps; both were fixed and reviewed again. A Rust fixture boundary was also removed.
- Verification on 2026-09-05: 54 top-level Python tests passed; 18 cell-type WDL tests passed with one Cromwell test skipped in discovery, then the Cromwell test passed separately using Womtool 87. All 32 WDLs passed MiniWDL and the Terra file-scope check. Actionlint and diff whitespace checks passed.
- The local synthetic splicing run passed. Proteomics helper checks ran, but local CLI startup stopped because this host lacks OlinkAnalyze. The image includes that package; the GitHub container gate remains pending.
- No local container builds, image publications, repository setting changes, or Terra jobs were performed. GitHub runtime tests and the protected manual publishing trial are still required before enabling automation.

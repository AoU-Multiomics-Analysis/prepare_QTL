# Stage Directory Policy Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Register new scripts by stage directory and validate supported WDL task image wiring without per-task name lists.

**Architecture:** Extend the existing trusted registry with directory ownership and runtime path mappings. Resolve literal script commands to stages, then use MiniWDL call bindings to verify that a stage image reaches each task. Preserve current source paths, defaults, and release behavior during this policy-only phase.

**Tech Stack:** Python 3.11, MiniWDL 1.15.0, PyYAML 6.0.2, unittest, existing GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-09-05-stage-script-layout-design.md`.

## Global Constraints

- Preserve scientific calculations, file formats, WDL 1.0, and Terra file localization.
- Keep existing stage IDs, image repositories, and image-input names in this phase.
- Do not change existing digests without publication and runtime tests.
- Do not build Docker images locally.
- Do not enable releases, publish images, push branches, merge PRs, or submit Terra jobs without the applicable user approval.
- Do not infer dependencies from arbitrary shell or R code.
- Static validation does not prove scientific correctness or call-cache reuse.

## Phase boundary

This plan implements the first independently testable policy change. No scripts
move and no image is built here. Existing WDLs must pass unchanged. The approved
separate-stage/shared-base design needs a subsequent build-graph plan, trusted
policy merge, and source/image migration. Do not mark the whole migration done
when these tasks finish.

Use the existing clean linked worktree on `codex/stage-script-layout`. The prior
CI fix is already on main. Read AGENTS.md and the spec. Do not use the private
App key or invoke setup as part of implementation.

## File responsibilities

| File | Responsibility |
| --- | --- |
| `ci/image-stages.yml` | Stage roots and existing legacy classifications |
| `ci/script_stages.py` | Pure script ownership and runtime-path resolution |
| `ci/plan_image_updates.py` | Directory selection and registry validation |
| `ci/wdl_stage_routing.py` | Static command classification and image forwarding |
| `tests/test_script_stages.py` | Directory ownership and path failure cases |
| `tests/test_wdl_stage_routing.py` | Synthetic task and nested-call regression cases |
| `tests/test_image_update_plan.py` | Planner selection behavior |
| `tests/test_repo_image_routing.py` | Repository-wide image contract |
| `tests/test_stage_image_inputs.py` | Cell-stage defaults, outputs, and forwarding |
| `.github/workflows/image-release-plan.yml` | Read-only checks for the new policy |
| `ci/test_release_images.py` | Run trusted checks against candidate source |
| `docs/image-release.md` | Contributor steps and migration limits |

## Task 1: Directory ownership and planner selection

**Files:** Create `ci/script_stages.py`, `tests/test_script_stages.py`; modify
`ci/image-stages.yml`, `ci/plan_image_updates.py`, `tests/test_image_update_plan.py`.

**Interfaces:**

```python
def source_patterns(stage: dict) -> list[str]: ...
def stages_for_script(config: dict, source_path: str) -> set[str]: ...
def validate_script_roots(config: dict) -> list[str]: ...
def source_for_runtime(config: dict, runtime_path: str, root: Path) -> str: ...
```

`source_patterns` combines legacy `sources` with `script_roots` converted to
`root/**`. `stages_for_script` also includes declared shared consumers.
`source_for_runtime` returns one repository-relative path or raises ValueError;
it never reads script contents or executes a script.

- [ ] Add `script_roots` to the four cell stages: estimation, fit, export,
  downstream beneath `scripts/cell_type_specific_expression/`. Add shared and
  tools globs to the existing shared/excluded policies respectively. Add
  `scripts/expression/prepare/` and `scripts/expression/rnaseqc/` roots. Treat
  existing single-stage script directories and Rust crate roots as stage roots.
  Preserve every current per-file/shared mapping for compatibility.
- [ ] Write failing tests using new, nonexistent paths: ownership is a path
  policy, so a proposed new file must be classifiable before it is committed.

```python
plan = plan_changes(config, [
    'scripts/cell_type_specific_expression/export/new_tool.R'])
self.assertEqual(plan['stages'], ['cell_export'])
self.assertEqual(plan['unmapped'], [])
self.assertEqual(stages_for_script(config,
    'scripts/cell_type_specific_expression/shared/new_helper.R'),
    {'cell_estimation', 'cell_fit', 'cell_export', 'cell_downstream'})
```

- [ ] Add tests for unknown roots, path traversal, absolute source paths,
  overlapping owned roots, unknown shared consumers, roots outside image build
  coverage, and unchanged legacy stage selection. Include new scripts under
  each existing stage root, not just export. An owned root must belong to one
  stage; shared roots use their declared consumer set.
- [ ] Run focused tests and observe failure before implementation:
  `python -m unittest discover -s tests -p test_script_stages.py -v` and
  `python -m unittest discover -s tests -p test_image_update_plan.py -v`.
- [ ] Implement the pure helpers and use `source_patterns` in planner stage
  matching. Include `validate_script_roots` errors in registry validation.
  Use normalized POSIX paths; reject `..` and empty components instead of
  silently normalizing them into a different policy location.
- [ ] Add runtime prefix rules for `/opt/prepare_qtl/scripts/` to `scripts/`.
  Resolve legacy `/tmp/<basename>.R` against the existing standard-image copy
  roots only, requiring exactly one matching tracked source. Resolve current
  RNA-SeQC paths with the same canonical prefix rule. Unknown/ambiguous paths
  fail. Keep explicit legacy runtime aliases when basename resolution is not
  unique; do not select the first match. Test quoted paths with spaces.
- [ ] Re-run the two focused suites, inspect the diff, and commit only Task 1
  files with message `Register stage directories with legacy path support`.

## Task 2: Script-based WDL stage verification

**Files:** Create `ci/wdl_stage_routing.py`, `tests/test_wdl_stage_routing.py`;
modify `tests/test_repo_image_routing.py`, `tests/test_stage_image_inputs.py`.

**Consumes:** Task 1 helpers and registry. **Produces:**

```python
def task_stages(task, config: dict, source_root: Path) -> set[str]: ...
def validate_routing(source_root: Path, policy_root: Path) -> list[str]: ...
```

Errors name the WDL file, task or call, and unresolved path/input. The validator
uses trusted policy from `policy_root` and WDL/script existence from
`source_root`; candidate source must not replace the policy.

- [ ] Write small WDL 1.0 fixture trees in TemporaryDirectory. Include a new
  task running a literal `Rscript /opt/prepare_qtl/scripts/cell_type_specific_expression/export/new.R`
  and a workflow forwarding `export_docker_image`. No task-name registration is
  supplied. Use literal sentinel image values to distinguish stage inputs even
  when their real digests are identical.

```python
self.assertEqual(validate_routing(source_root, policy_root), [])
# Replace only the call binding with fit_docker_image in the fixture.
self.assertTrue(any('stage' in error.lower()
                    for error in validate_routing(source_root, policy_root)))
```

- [ ] Add fixtures for imported child workflows, scatters, conditionals,
  multiple scripts of one stage, hard-coded runtime images, omitted forwarding,
  a wrong intermediate binding, quoted script paths, unknown scripts, and a
  dynamic interpreter target. A dynamic path must fail rather than disappear
  from the scan. Check valid shell flags and continuation lines used by current
  repository commands. Do not execute task commands.
- [ ] Run `python -m unittest discover -s tests -p test_wdl_stage_routing.py -v`
  and observe the expected failures.
- [ ] Implement the literal-command subset explicitly: Rscript and python/
  python3 file invocations with repository runtime paths. Use MiniWDL command
  parts and shell tokenization; preserve unresolved placeholders so they cannot
  be misread as safe literals. Reject unsupported dynamic script targets with a
  specific diagnostic. Inline interpreter programs and external binaries need
  explicit existing contracts; they are not inferred from a filename.
- [ ] Derive image-input stage identity from the registered pin targets and
  recursively propagate sentinel identities through WDL call bindings to task
  runtimes. A shared helper can be used by one of its declared consumers, but
  a task executing distinct owned stages must fail for review. Avoid executing
  input file expressions or workflow file-writing functions during validation.
- [ ] Replace task-name mappings only for supported script-backed tasks. Keep
  a narrowly named explicit contract table for existing inline/external tasks.
  Do not delete immutable-default or unregistered-runtime checks. Extend the
  cell tests so new correctly wired script-backed calls do not require edits
  to their hard-coded call-name sets. Preserve stage_images output tests.
- [ ] Run the new suite and both existing routing suites on the actual repo.
  Every current WDL must pass unchanged. Check sentinel tests still fail for
  swapped fit/export inputs even though they share an image repository today.
- [ ] Commit Task 2 files as `Validate task image stages from script paths`.

## Task 3: Trusted CI integration and contributor instructions

**Files:** Modify `.github/workflows/image-release-plan.yml`,
`ci/test_release_images.py`, `docs/image-release.md`; extend
`tests/test_release_integration.py` where it covers trusted test invocation.

**Consumes:** `validate_routing(source_root, policy_root)` and the new tests.
No build-graph or Dockerfile changes are included in this task.

- [ ] Add regression coverage showing candidate policy cannot override trusted
  routing. Construct separate source/policy fixture roots, swap only candidate
  metadata, and require the bad task binding to remain rejected.

```python
errors = validate_routing(candidate_root, trusted_root)
self.assertTrue(errors)
# A matching candidate-owned registry must not turn the result into success.
```

- [ ] Run the focused integration/routing tests and observe the missing trust
  enforcement or integration before implementing it.
- [ ] Add explicit read-only CI commands for the new test files after existing
  pinned Python dependency installation. Keep job permissions, release flags,
  image build triggers, and approval environments unchanged.

```yaml
- name: Test stage directory and task routing policy
  run: |
    python -m unittest discover -s tests -p test_script_stages.py -v
    python -m unittest discover -s tests -p test_wdl_stage_routing.py -v
```

- [ ] Ensure `ci/test_release_images.py` invokes the routing checks from trusted
  main against candidate source, as it already does for image routing. Do not
  run candidate test files to approve candidate code.
- [ ] Document: place a new script in an owned directory; wire the stage image
  in WDL; add a runtime test. New stages still need registry/pin/runtime-gate
  registration. Shared code invalidates declared consumers. Mark legacy paths
  transitional and explain that this PR does not yet split or publish images.
- [ ] Run complete local verification:

```bash
python -m unittest discover -s tests -p 'test_*.py'
python ci/plan_image_updates.py --validate
python scripts/check_wdl_file_scope.py workflows
actionlint .github/workflows/image-release-plan.yml
git diff --check
```

Use the existing temporary Python runtime with MiniWDL/PyYAML/cryptography if
available. Callback tests may need approved local socket access; GitHub calls
remain fake. Do not build or run Docker locally for these tests.

- [ ] Request task and final whole-change review. Fix blocking findings and
  re-run covering tests. Confirm there are no changes to analysis scripts,
  WDL defaults, or release enable flags. Commit as
  `Run directory routing checks in trusted release validation`.
- [ ] Report this phase's results and ask permission to push/open its PR. The
  next work is a separate shared-base/stage-image build-graph plan; source moves
  must wait until their trusted policy and tests are on main.

## Plan self-review

Directory ownership, existing roots, shared consumers, and legacy path handling
are Task 1. Static task classification and correct image forwarding are Task 2.
Trusted candidate checks and documentation are Task 3. Runtime compatibility,
base dependency locking, base/child fingerprinting, image publication, source
moves, and pin-repository migration are deliberately assigned to the next
build/migration plan, not claimed by this policy-only implementation.

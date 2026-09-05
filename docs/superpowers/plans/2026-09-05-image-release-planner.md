# Repo-wide image release planner implementation plan

**Goal:** Add a report-only dependency planner covering repository workflows,
without changing runtime images or granting automated write access.

**Architecture:** A YAML registry maps local source changes to planned release
groups and image builds. A read-only Python command prints a Markdown report.
GitHub Actions runs tests and reports the PR change set; no builds or pin writes
are performed by the new workflow.

**Tech stack:** Python 3.11, PyYAML, unittest, GitHub Actions.

**Spec:** `docs/stage-image-pinning-draft.md`. This implements only its initial
dependency-selection component. The repo-wide extension includes four locally
built images and explicitly classifies genotype images as external.

## Constraints

- Preserve all current WDLs, image defaults, and existing build behavior.
- No Docker build locally; no Terra submissions; no registry/repository writes.
- Keep shared dependencies conservative. New unmapped source paths fail closed.
- A release-group report is not proof of cache reuse or of complete dependency
  analysis. Mixed-image execution tests remain required before pin automation.
- CI configuration is YAML; this is not a JSON argument wrapper for WDL tasks.

## Task 1: Dependency selection

Files: `ci/image-stages.yml`, `ci/plan_image_updates.py`,
`tests/test_image_update_plan.py`.

Interface: `plan_changes(config, changed_paths)` returns a dictionary with
`builds`, `stages`, `wdl_checks`, `test_changes`, and `unmapped` lists.

- [ ] Add table-driven unittest cases for a filtering-only edit, shared I/O edit,
  a standard expression edit, image environment edits, WDL-only edits, test-only
  edits, unknown source files, and changes spanning two image families.
- [ ] Run `python -m unittest discover -s tests -p test_image_update_plan.py` and
  confirm the missing planner fails before implementation.
- [ ] Implement selection from explicit paths and positive glob patterns. Every
  matching rule contributes stages. Image build dependencies are distinct from
  stage dependencies: the standard image copies all scripts but does not require
  every standard-image stage to adopt every rebuilt version.
- [ ] Add CLI `--base`, `--head`, `--config` and `--validate`. Read changed paths
  with `git diff --no-renames --name-only -z BASE...HEAD`. Decode each path without
  shell interpolation. Include both old and new paths for renames by disabling
  rename detection. Print a report; exit nonzero for unknown source files.
- [ ] Validate image/stage references, source coverage, and WDL coverage against
  tracked repository files. Untracked draft files must not change CI results.

## Task 2: Read-only GitHub report

File: `.github/workflows/image-release-plan.yml`.

- [ ] Run the planner tests on pull requests and manual dispatch.
- [ ] Check out full history with no persisted credentials. Use a read-only
  contents token; no package writes, secrets, bot commits, or PR comments.
- [ ] Compare the PR base SHA with the head SHA, using environment variables and
  quoted shell arguments. Write the Markdown report to the job summary.
- [ ] On manual dispatch validate coverage and compare the current commit with
  its parent. Clearly label this as a report, not an image publication.
- [ ] Run `actionlint`, planner tests, and `git diff --check`.

## Task 3: Explain the trial and remaining release work

File: `docs/image-release-planner.md`.

- [ ] Document a local filtering-only simulation and its expected groups.
- [ ] Explain that WDL-only and test-only edits do not request pin updates.
- [ ] Explain external genotype images, broad standard-image build dependencies,
  shared-module conservatism, and how to add a source path to the registry.
- [ ] List remaining work: verify initial digests, route stage inputs, test mixed
  images, add the trusted updater, configure App permissions, and enforce checks.
- [ ] Do not call the full pinning system complete after this report-only trial.

# Task-specific image pins

Each repository-built task has a release unit in `ci/image-stages.yml`. The
`stages` key is retained for compatibility with the release controller, but
its entries now identify individual tasks. `test_group` selects an existing
runtime test suite; it does not tie their image pins together.

`ci/release-pins.yml` records the task's literal image default and each workflow
input that forwards it. The updater changes all copies for an affected task
and leaves every other task's digest unchanged. Docker builds remain shared:
one candidate image per affected image family, not one build per task.

## Update rules

- An entrypoint edit updates the tasks that execute that entrypoint.
- A shared helper edit updates every task whose dependency list contains it.
- A Dockerfile or environment edit updates all tasks using that environment.
- Inline WDL command edits change the task command used by call caching. They
  do not require an image rebuild when no image contents changed.
- New runtime source files and new repository-managed tasks require explicit
  dependency records. Validation rejects missing records and missing helpers.
- Third-party genotype and tensorQTL images remain manual: this repository does
  not build their software. The three repository-built trans-LD tasks are managed.

The cell-type entrypoints have checked-in `modules/*.txt` helper lists. The
scoped loader (the dependent source PR) loads only those modules. A conservative
transitive symbol scan checks these lists, including constants and helper
functions. Literal R helper imports and relative Python imports are also checked.
Dynamic imports require review; a dependency list is not a proof of scientific
compatibility between versions.

## Migration

The policy PR preserves the current task command text and existing image digests.
It replaces workflow-wide/stage-wide image inputs with task-specific inputs such
as `expression__filter_expression_genes_image` and `hspe__run_hspe_batch_image`.
Users with image overrides in Terra must remove the old overrides or map them to
the applicable new inputs. Other workflow inputs retain their names and types.
The cell-type image metadata output is now `task_images`, keyed by task unit.

Merge the policy/WDL migration first, then retarget the dependent scoped-loader
PR to main and label that source PR `release-ready`. Do not label the policy PR:
the trusted release controller intentionally rejects mixed policy/source releases.
The loader change triggers one initial cell-type image refresh. Later unrelated
helper changes can preserve the other tasks' pins. Rebase this migration after
any concurrent image release (including the constrained TCA PR) to retain its
published digest values before merging.

The starting pins are a verified current baseline, not a reconstruction of the
historical last commit that changed each script. Subsequent successful releases
maintain task-specific pins prospectively. A digest change can prevent a cache
hit, but an unchanged digest does not guarantee one: inputs, commands, outputs,
and access to cached files must also match. No complete workflow was run on Terra
for this migration.

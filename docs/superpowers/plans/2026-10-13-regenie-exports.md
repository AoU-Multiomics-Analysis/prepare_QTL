# REGENIE export implementation plan

> **For agentic workers:** Use superpowers:executing-plans to implement this plan.

**Goal:** Add REGENIE phenotype and covariate outputs to the current QTL preparation workflows.

**Architecture:** One streaming exporter and one shared WDL task. Modality workflows
call the task after existing preparation and covariate merging. Public wrappers
pass through typed outputs.

**Tech Stack:** Python 3.11 standard library, WDL 1.0, MiniWDL, GitHub Actions,
the existing pinned REGENIE image.

**Spec:** [Design](../../regenie-exports-design.md)

## Global constraints

- WDL 1.0, Terra managed Cromwell, command logging, safe shell quoting.
- No workflow-scope file creation; preserve File types through localization.
- No local Docker build or Terra submission.
- Keep current outputs and source participant data out of Git.

## Review focus

- Different sample order and partial covariate overlap must align by ID.
- Unsafe trait IDs and duplicates must not mix traits or overwrite outputs.
- Missing, infinite, and malformed cells must have explicit behavior.
- Localized files with shell punctuation and absent optional files must work.
- Large feature counts must be processed in bounded batches.

### Task 1: Export and integrate all preparation outputs

**Files:** tools/regenie_export/export_regenie.py; workflows/common/ExportRegenie.wdl;
modality and public wrapper WDLs; tests/regenie_export/; docs/regenie-exports.md;
.github/workflows/regenie-exports.yml; .dockstore.yml; ci/image-stages.yml.

**Interface:** BED plus optional sample-row or TensorQTL covariate File, batch size,
output prefix. Outputs: Array[File] phenotype batches; optional File covariates;
File trait map, batch table, exclusions, log. The batch table is metadata; WDL
discovers phenotype Files with glob, not paths parsed from the metadata.

- [ ] Write CLI and localization regression tests. Assert exact transpose and
  numeric text, stable safe IDs, missing NA, reordered overlap, categorical
  covariates, duplicates, nonfinite values, empty input, and cloud URI failures.
- [ ] Run the tests and confirm failure because the exporter does not exist.
- [ ] Implement the exporter, literal WDL program, typed task, workflow calls,
  and pass-through outputs. Add common helper registration and routing policy.
- [ ] Run the exporter and localization tests, all descriptor checks, file-scope
  regression check, existing Python/R tests available locally, and the CI smoke.
- [ ] Document batch selection and output bindings. Commit; obtain a fresh code
  review, fix material findings with regression tests, push and open a PR.

## Execution record

Native execution. The user requested the code change; no further design or
worktree approval is needed for this reversible branch.

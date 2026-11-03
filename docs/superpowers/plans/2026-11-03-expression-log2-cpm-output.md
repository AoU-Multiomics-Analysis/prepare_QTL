# Expression Log2-CPM Output Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans in this session.

**Goal:** Export unscaled log2 expression with the same genes and samples as the scaled BED.

**Architecture:** Return retained sample IDs from the existing BED writer. Select
these columns from log2 expression and use the same BED writer without an
additional connectivity check. Expose a typed WDL File and retain existing PCs.

**Tech Stack:** R, tidyverse, WDL 1.0, miniwdl, existing GitHub Actions image tests.

**Spec:** `docs/superpowers/specs/2026-11-03-expression-log2-cpm-output-design.md`

## Global Constraints

- No workflow-scope file creation; preserve File typing and command logging.
- Keep all existing outputs and scaled PCA behavior.
- Do not build Docker locally or submit Terra jobs.
- Reuse this checkout on `codex/expression-log2-cpm-output`; preserve unrelated files.

## Review Focus

- Existing log2 input must not receive a second log transform.
- A connectivity outlier must be excluded from both scaled and log2 BEDs.
- Count-input annotation joins must give identical genes and coordinates.
- Numeric and shell-special sample IDs must retain their exact spelling and order.
- Signed input must retain the existing signed-log convention, without claiming ordinary CPM fold changes.

## Task 1: Expression output and workflow interface

**Files:** `scripts/expression/prepare/PrepareExpression.R`,
`workflows/expression/prepare_eQTL.wdl`,
`tests/test_prepare_expression_log2_cpm.R`, `tests/test_expression_input_wdl.py`.

**Interface:** `write_expression_bed()` returns retained sample IDs invisibly;
the task exposes `Log2CpmExpressionBed`, and the workflow exposes
`Log2CpmBedFile`.

- [x] Extend existing expression tests for the new BED, all three input modes,
  and a real connectivity outlier; run them before implementation and confirm
  the output is missing.
- [x] Test typed task/workflow output discovery and retained cloud File paths.
- [x] Add the unscaled log2 output with the scaled sample IDs. Update WDL outputs
  and lifecycle logs. Do not add a new PCA call.
- [x] Run expression tests, output/localization tests, WDL syntax and file-scope checks.

## Task 2: Input documentation and verification

**Files:** `README.md`, `docs/scripts.md`,
`docs/phenotype-normalization-filtering.md`, `docs/molecular-qtl-workflows.md`.

- [x] Document values, sample filtering, supported input modes, and the precise
  outlier manifest fields. Explain separate fixed covariates and all-PC input.
- [x] Verify that existing GitHub Actions expression tests include the extended
  smoke test and that source changes select the expression image stage.
- [x] Review the final diff, run the relevant checks, and record any baseline failures.

## Execution notes

The checkout was fast-forwarded from `582bdb6` to its existing `origin/main`
reference `fd100b0` before edits. Untracked files were preserved. The baseline
`tests/testthat` suite fails before this change in
`test_cell_type_specific_documentation.R` (stale input/manifest expectations)
and `test_eqtl_runtime_contract.R` (old `:main` image and forwarding expectations).
These unrelated checks are outside this change. Image builds and pin publication
remain release work; no new image digest will be invented.

Verification: the full expression CLI integration test, reference/signed CPM
integration, and sample-list integration all passed. The three expression WDL
tests and five CI path-selection tests passed. Targeted R PC/output tests passed.
MiniWDL checked expression and its cell-type wrapper; the workflow-scope check
passed for all 39 WDL files and the logging check passed. Dockstore validation
passed for 23 descriptors. The full R suite still has the same 14 baseline
assertion failures in the two files named above.

Read-only review confirmed the implementation and downstream mappings. Review
notes were addressed in the guides: fixed covariates are numeric, signed CPM is
an explicit exception, and existing PC ID rewriting is documented. Image pins
must be updated in all registered expression-stage targets before release.

The user then requested a GitHub push, Actions image build, and image-pin update.
The feature branch was rebased onto current main `c297869`. PR #90 contains
source commit `ffc70bf` and release-App pin commit `d664b65`. Stage Image Release
run `41726806674` passed. Its verified standard image is
`ghcr.io/aou-multiomics-analysis/prepare_qtl@sha256:f9329733aedfb159cf3d38e837ed195bd2db1de9bb56b6f46673814b33019765`.
The source fingerprint matches the pushed branch. All three registered expression
pins were updated: task, expression workflow, and cell-type wrapper. No other
stage pin changed. GitHub Actions source checks were requested for the pinned
revision. No Terra jobs were submitted.

The first source-check run stopped in three existing cell-type output-contract
assertions after main added REGENIE exports. The wrapper differs from main only
in its expression image pin. Updated the test to include the published REGENIE
file arrays and allow raw REGENIE exports while still rejecting raw BED and
residualized BED outputs. This test-only correction needs no image rebuild.

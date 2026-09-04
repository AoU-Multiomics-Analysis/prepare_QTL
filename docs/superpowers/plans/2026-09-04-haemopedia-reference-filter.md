# Haemopedia Reference Filter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Filter each exported TCA BED independently and compare mean log expression with Haemopedia before eQTL preparation.

**Architecture:** Prepare the optional raw-count reference once. A chunked post-export R task writes filtered BEDs, tables and plots; the existing QTL scatter consumes those BEDs. Original BEDs and the model remain unchanged.

**Tech Stack:** R, tidyverse, edgeR, ggplot2, WDL 1.0, MiniWDL, GitHub Actions.

**Spec:** `docs/superpowers/specs/2026-09-04-haemopedia-reference-filter-design.md`

## Global Constraints

- Remove a negative gene only from the affected cell type's BED.
- Mean expression is `mean(log2(CPM + 1))`; the threshold is strictly `> 0.01` in both datasets.
- Filtered BEDs remain linear CPM; never clip values or refit TCA.
- Equal weighting of all mapped reference sample columns; normalize the full raw count matrix before gene subset.
- Residual filtering is one pass and disabled when its optional cutoff is absent.
- WDL 1.0; no workflow-scope file-writing functions. Cloud-path manifest strings must not be localized.
- Existing global runtime retries; no local Docker builds or cloud submissions.
- R/tidyverse, minimal plots without titles or subtitles.
- Work only in the existing isolated `codex/haemopedia-reference-filter` worktree. Do not push or open a PR during implementation.

## Task 1: Reference and per-cell-type filtering engine

**Files:**
- Create `scripts/cell_type_specific_expression/R/reference_expression.R` (normalization and mapping).
- Create `scripts/cell_type_specific_expression/R/reference_filter.R` (statistics and filter decisions).
- Create `scripts/cell_type_specific_expression/R/reference_filter_io.R` (chunked BED IO and reports).
- Create `scripts/cell_type_specific_expression/R/reference_filter_plots.R` (plots).
- Create `scripts/cell_type_specific_expression/prepare_haemopedia.R` and `filter_cell_type_beds.R` (CLI).
- Create `tests/cell_type_specific_expression/testthat/test-reference-filter.R`.

**Interfaces:**
- `prepare_haemopedia.R COUNTS OUTPUT_DIR` writes `reference_summary.tsv.gz`, `reference_samples.tsv`, `reference_metadata.json`. Summary includes `gene_id`, `cell_type`, `n_samples`, `mean_log2_cpm1`, `median_log2_cpm1`. Metadata records mapping caveats, edgeR version, input provenance and full-matrix normalization; samples table holds library sizes and norm factors. Validate matching keys and subtype completeness.
- `filter_cell_type_beds.R CONFIG_JSON OUTPUT_DIR`. Config keys: `inventory` (localized TSV), `bed_paths` (array of localized BED paths), `reference_summary` (optional localized summary path or null), `min_mean_log2_cpm1` (number), `residual_cutoff` (number or null).
- Filter output files: `filtered_beds.txt` (task-relative paths, in inventory order), `filtered_inventory.tsv` (existing exact scatter inventory schema), `negative_summary.tsv.gz`, `gene_comparison.tsv.gz`, `filter_metrics.tsv`, and PDFs under `plots/`. BEDs live under `beds/` with names `<slug>.filtered.bed.gz`. Keep checksums and dimensions correct. Reports carry all original gene/cell pairs and explicit unavailable statuses.
- Reference matching failures cannot be treated as expression zero. No-reference cell types get negative filtering only. Missing entire reference gets negative filtering only.
- Per-cell regression uses OLS with intercept; `stats::rstandard` convention. Undefined residuals are NA, and cause an error if a removal cutoff was supplied. Save baseline and retained-set descriptive metrics separately.

- [ ] Write focused tests for the specified functions/CLIs, watching them fail before production code. A hand-checked scale fixture is:

```r
x <- matrix(c(0, 3, 15, -1e-12, 2, 4), nrow = 2, byrow = TRUE)
# First gene mean log2(CPM+1) = 2; second gene must be excluded.
# In another cell type with abs(x), both genes pass the negative check.
```

  Include full-matrix edgeR parity, weighted pooling, strict threshold equality,
  unsupported reference cell types, ID collisions, zero/nonfinite values,
  regression and one-pass cutoff, insufficient genes and perfect fit. Run real
  CLI tests with filenames containing spaces, chunk boundaries and different
  retained gene sets. Check BED coordinates/sample order/values and output plots.
- [ ] Run the focused test with `/private/tmp/celltype-deconvolution-r45/bin/Rscript -e 'testthat::test_dir("tests/cell_type_specific_expression/testthat", filter="reference-filter", reporter="summary", stop_on_failure=TRUE)'`; if edgeR is unavailable in that runtime, use the existing system R with installed edgeR and report the exact runtime.
- [ ] Implement the spec using small pure calculation functions plus two-pass chunked IO: first retain only per-gene summaries for regression, then stream unchanged retained expression rows into gzip BEDs. Validate inventory checksums and gene/sample counts, all IDs and coordinates. The original raw summary enforces identical coordinates across raw BEDs; new filtered outputs must permit unequal gene sets.
- [ ] Generate negative prevalence plots, top-100 negative-fraction heatmap with deterministic tie breaks, per-cell baseline scatter/residual plots and separate post-residual plots when enabled. Always create at least a negative overview PDF even without reference. Use explicit statuses for empty comparisons.
- [ ] Run the focused tests and existing gene-summary/integration tests, self-review, and commit only these owned files. Report RED and GREEN evidence.

## Task 2: WDL routing and cloud manifest columns

**Files:**
- Create `workflows/cell_type_specific_expression/tasks/reference_filter.wdl`.
- Modify `workflows/cell_type_specific_expression/deconvolution.wdl`, `prepare_cell_type_eQTL.wdl`, `tasks/integration.wdl`.
- Modify `scripts/cell_type_specific_expression/build_qtl_manifest.R`, `R/integration.R`.
- Create `tests/cell_type_specific_expression/test_reference_filter_wdl.py`; extend `testthat/test-qtl-manifest.R`.

**Interfaces:**
- Both workflows expose `File? haemopedia_counts`, `Float reference_min_mean_log2_cpm1 = 0.01`, `Float? reference_residual_cutoff`.
- Prepare reference conditionally; filter original BEDs unconditionally. Pass config through task-scoped `write_json(object {...})`. CLI signatures/output paths are defined in Task 1.
- New workflow outputs: `filtered_cell_type_beds`, `filtered_cell_type_bed_inventory`, `negative_expression_summary`, `reference_gene_comparison`, `reference_filter_metrics`, `reference_filter_plots`, `reference_filter_log`; optional `haemopedia_reference_summary`, `haemopedia_reference_samples`, `haemopedia_reference_metadata`.
- Preserve all old public outputs. Existing `PrepareScatterInputs` consumes filtered BEDs and filtered inventory. No other changes to PrepareExpression transforms.
- Final QTL manifest adds `source_cpm_bed`, `filtered_cpm_bed`, `negative_expression_summary`, `reference_gene_comparison`, `reference_filter_metrics`. Paths are full upstream cloud URLs. Per-cell paths must be aligned by slug rather than trusting unrelated array orders. Shared report paths may repeat.
- Maintain old direct manifest API/CLI callers by making new metadata optional and validating it only when provided. No need to alter the older deconvolution manifest schema.

- [ ] Add failing typed-WDL tests verifying dependency routing and task boundary types, and a real R manifest test with two cells and shared gs:// reports. Example input paths:

```r
source_beds <- c("gs://bucket/raw/b_cells.bed.gz", "gs://bucket/raw/cd4_t_cells.bed.gz")
filtered_beds <- c("gs://bucket/filtered/b_cells.filtered.bed.gz", "gs://bucket/filtered/cd4_t_cells.filtered.bed.gz")
# Returned rows retain all previous QTL paths and both new paths, with a shared
# gs://bucket/reports/filter_metrics.tsv allowed in both rows.
```

- [ ] Run targeted tests to observe missing routing/columns, implement task commands with start/failure/completion logging and safe task-scoped serialization, and expose outputs. Reuse gene-summary runtime settings for the streaming task and explicit bounded reference preparation memory.
- [ ] Run typed-WDL checks, `scripts/check_wdl_file_scope.py workflows`, logging checks, and R manifest tests. Update existing contract tests only where the approved routing intentionally changes.
- [ ] Self-review and commit owned files. Provide RED/GREEN evidence and identify any CI-facing schema changes to Task 3.

## Task 3: CI fixture, user guide, and full validation

**Files:**
- Modify `.github/workflows/cell-type-specific-expression-ci.yml` and test smoke scripts under `tests/cell_type_specific_expression/smoke/`.
- Add a tiny raw-reference fixture generator under `tests/cell_type_specific_expression/` and use its output in one existing end-to-end smoke mode.
- Update cell-type workflow documentation and input examples in their existing directories.
- Modify `envs/CellTypeSpecificExpression/environment.yml` only if an explicit pinned edgeR dependency is absent (inspect actual filename first).

**Interfaces:**
- Use Task 1's CLI outputs and Task 2's workflow output names. One existing e2e mode exercises the optional reference; the other exercises no reference. Do not create another costly full-workflow mode solely for this feature.
- Fixture raw counts cover all mapped prefixes and synthetic gene IDs, with deterministic positive counts; include boundary cases in focused R tests rather than making the e2e TCA fit unstable.

- [ ] Add assertions that actual filtered BEDs are nonnegative, samples/coordinates/values match retained source rows, manifest paths identify the correct raw and filtered cell, and optional reference metrics exist. Update old exact-schema smoke expectations to preserve original columns while allowing approved additions.
- [ ] Wire focused tests and fixture creation into GitHub Actions before MiniWDL runs, using the built container. Continue building containers only in GitHub Actions.
- [ ] Document raw counts versus CPM input, full-matrix TMM normalization, sample mapping, mean-of-logs convention, strict per-cell negatives, residual cutoff off by default, one-pass removal, unsupported reference populations and original-output preservation. Add sample JSON keys:

```json
{
  "PrepareCellTypeEqtlWorkflow.haemopedia_counts": "gs://YOUR_BUCKET/GSE115736_Haemopedia-Human-RNASeq_raw.txt.gz",
  "PrepareCellTypeEqtlWorkflow.reference_min_mean_log2_cpm1": 0.01
}
```

- [ ] Run the available R suite, MiniWDL type checks, static Terra/logging regressions, and lint. Inspect generated PDFs. Do not claim a Terra or GitHub Actions run occurred locally.
- [ ] Review all changes against the spec, record limitations, and commit. No push/PR without new user authorization.

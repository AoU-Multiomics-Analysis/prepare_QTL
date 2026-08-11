# Optional methylation-depth correlation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the PacBio methylation-versus-sequencing-depth correlation diagnostic optional while preserving the existing metadata schema and all other methylation QC behavior.

**Architecture:** Add one workflow-level Boolean defaulting to `true`, propagate it through every methylation cohort entry point, and gate the diagnostic in both the Rust chromosome merger and the R manifest merger. Disabled runs will write stable-schema sentinel values (`NA` correlation and `0` contributing samples) without changing filtering or phenotype generation.

**Tech Stack:** WDL, Rust 2021 with clap/csv, R with data.table/optparse, Cargo tests.

## Global Constraints

- The default behavior must remain unchanged.
- Only the methylation-versus-sequencing-depth diagnostic is optional.
- Minimum coverage, extreme-coverage filtering, local CpG correlation, connectivity filtering, site filters, and phenotype outputs must not change.
- Output metadata column names and order must remain unchanged.

---

### Task 1: Add failing regression coverage for the Rust diagnostic toggle

**Files:**
- Modify: `rust/methylation_merge/src/main.rs`
- Test: `rust/methylation_merge/src/main.rs` inline unit tests

**Interfaces:**
- The test will exercise a small helper that summarizes optional correlation vectors.
- The helper will return `(0, None)` when disabled and the existing count/correlation when enabled.

- [x] **Step 1: Write the failing tests**

Add an inline `#[cfg(test)]` module with tests named `disabled_depth_correlation_returns_no_observations` and `enabled_depth_correlation_returns_spearman_summary`. The first should call the helper with `false` and assert count `0` and `None`; the second should call it with `true` and three nonconstant paired values and assert count `3` plus a correlation near `1.0`.

- [x] **Step 2: Run the focused tests and verify the expected failure**

Run:

```bash
cargo test --manifest-path rust/methylation_merge/Cargo.toml disabled_depth_correlation_returns_no_observations -- --nocapture
```

Expected result: compilation/test failure because the helper does not yet exist.

### Task 2: Implement the Rust toggle and wire the current WDL path

**Files:**
- Modify: `rust/methylation_merge/src/main.rs`
- Modify: `workflows/methylation/cohort_aggregation.wdl`
- Modify: `workflows/methylation/AggregateMethylationCohort.wdl`

**Interfaces:**
- Rust CLI adds `--skip-coverage-methylation-correlation`.
- `process_site` receives a Boolean `compute_coverage_methylation_correlation`.
- Current Terra cohort WDLs expose `ComputeCoverageMethylationCorrelation = true` and pass the skip switch only when false.

- [x] **Step 1: Add the minimal Rust implementation**

Add the helper used by Task 1, add the clap flag with default false for the skip switch, and pass `!args.skip_coverage_methylation_correlation` into `process_site`. Allocate and populate correlation vectors only when computation is enabled. Write `0` and `NA` in the metadata fields when disabled.

- [x] **Step 2: Run the focused tests and verify they pass**

Run:

```bash
cargo test --manifest-path rust/methylation_merge/Cargo.toml -- --nocapture
```

Expected result: all Rust tests pass.

- [x] **Step 3: Add the WDL Boolean and conditional command argument**

Add `Boolean ComputeCoverageMethylationCorrelation = true` to `AggregateMethylationCohort.wdl`, pass it into `CohortAggregation.AggregateMethylationData`, add it to `cohort_aggregation.wdl`'s `MergeMethylationChromosome` input, and append:

```wdl
~{if !ComputeCoverageMethylationCorrelation then "--skip-coverage-methylation-correlation" else ""}
```

to the Rust command line.

### Task 3: Implement the R-path toggle and wire legacy WDL entry points

**Files:**
- Modify: `scripts/methylation/MergeMethylationCohort.R`
- Modify: `workflows/methylation/merge_methylation.wdl`
- Modify: `workflows/methylation/AggregateMethylationCohortArrays.wdl`

**Interfaces:**
- R CLI adds `--SkipCoverageMethylationCorrelation` defaulting to `FALSE`.
- WDL Boolean remains named `ComputeCoverageMethylationCorrelation` and is translated to the R skip switch.

- [x] **Step 1: Add the R skip option**

Add the optparse Boolean switch. When set, do not create normalized log-coverage vectors for the diagnostic and write `0`/`NA` for the two metadata fields. When unset, preserve the current Spearman calculation exactly.

- [x] **Step 2: Wire the manifest/shard WDL**

Add the Boolean input to `merge_methylation.wdl`, pass it through the nested cohort workflow, and conditionally append `--SkipCoverageMethylationCorrelation` to the R merger command.

- [x] **Step 3: Wire the legacy array WDL**

Add the same Boolean input and pass it to each R chromosome-merge task, preserving the existing default behavior.

### Task 4: Document and validate the behavior

**Files:**
- Modify: `docs/methylation-qtl.md`
- Modify: `docs/scripts.md`

- [x] **Step 1: Document the new input and sentinel outputs**

State that the option defaults to enabled, affects only the diagnostic correlation, and produces `NA`/`0` in the metadata columns when disabled.

- [x] **Step 2: Validate WDL syntax and source behavior**

Run the repository's available WDL validation command for all modified WDL files. Run:

```bash
cargo fmt --manifest-path rust/methylation_merge/Cargo.toml -- --check
cargo test --manifest-path rust/methylation_merge/Cargo.toml -- --nocapture
git diff --check
```

Confirm the diff contains no changes to coverage thresholds, extreme-coverage logic, site filtering, connectivity, or phenotype generation.

### Task 5: Commit and push the implementation

**Files:**
- Commit all changes in the isolated branch.

- [x] **Step 1: Review status and diff**

Run:

```bash
git status --short
git diff --stat
git diff -- workflows/methylation rust/methylation_merge/src/main.rs scripts/methylation/MergeMethylationCohort.R docs/methylation-qtl.md docs/scripts.md
```

- [ ] **Step 2: Commit**

```bash
git add docs scripts/methylation workflows/methylation rust/methylation_merge/src/main.rs
git commit -m "feat: make methylation depth correlation optional"
```

- [ ] **Step 3: Push the branch**

```bash
git push -u origin codex/optional-depth-correlation
```

Report the pushed branch only after the command succeeds.

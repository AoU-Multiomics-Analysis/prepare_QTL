# Cell-Type-Specific eQTL Integration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `prepare_QTL` the canonical home of the whole-blood dtangle/TCA pipeline and add one WDL 1.0 workflow that scatters QTL preparation across retained cell-type expression BEDs.

**Architecture:** Migrate the reviewed deconvolution implementation at `CellTypeDeconvolution` commit `9a98d2e` into isolated workflow and script directories. Keep `prepare_eQTL.wdl` as the single-matrix implementation, call it once per validated cell type, and collect aligned outputs into a portable basename manifest plus authoritative WDL file arrays.

**Tech Stack:** WDL 1.0, R 4.5, tidyverse, dtangle 2.0.10, TCA 1.2.1, MiniWDL 1.15.0, Micromamba, Docker, GitHub Actions, Dockstore, Terra/Cromwell

**Spec:** `docs/superpowers/specs/2026-09-02-cell-type-specific-eqtl-integration-design.md`

## Global Constraints

- `prepare_QTL` is the only maintained source after migration.
- Preserve the reviewed deconvolution behavior from commit `9a98d2e`; do not redesign dtangle, LM22 grouping, TCA, or BED reconstruction.
- All new or modified WDL sources must use `version 1.0`.
- The expression input is a coordinate-preserving BED with finite, nonnegative, linear CPM values.
- Use `Float log2_pseudocount = 0.0`; reject negative or non-finite values.
- When the pseudocount is zero, reject zero CPM. When it is positive, allow zero CPM. Always reject negative CPM.
- Apply `log2(CPM + log2_pseudocount)` exactly once to both the dtangle and TCA views.
- Do not apply `log2(CPM + 1)` again in scattered eQTL preparation.
- Keep `AdditionalCovariates` required in the end-to-end workflow and optional in the existing single-matrix workflow.
- Remove connectivity outliers independently for each cell type and independently for INT and scaled branches.
- Do not add residualized BED files to the cell-type manifest.
- Manifest file fields use stable basenames; ordered `Array[File]` outputs are authoritative.
- Use one global `preemptible_attempts` value and one global `max_retries` value across the end-to-end path while retaining task-specific CPU, memory, and disk settings.
- Every new or modified WDL command block must log `stage`, `start_time`, `completion_time`, `dimensions`, and `outputs`.
- R additions must use tidyverse syntax.
- The dedicated container must use a pinned Micromamba base and pinned conda-forge/bioconda dependencies.
- Do not build Docker locally. Build and smoke test containers in GitHub Actions.
- Implement on `codex/cell-type-specific-expression`, based on `feat/log2-cpm-bed-input`; both branches are intended to merge into `main`.

---

## File structure

### Migrated deconvolution implementation

- Create `scripts/cell_type_specific_expression/R/*.R` for the nine reviewed R modules plus workflow validation.
- Create `scripts/cell_type_specific_expression/*.R` for the seven reviewed command-line programs and fixture generator.
- Create `workflows/cell_type_specific_expression/tasks/*.wdl` for the four reviewed task modules.
- Create `workflows/cell_type_specific_expression/deconvolution.wdl` as the standalone public workflow.
- Create `tests/cell_type_specific_expression/testthat/` and `tests/cell_type_specific_expression/smoke/` for migrated and new tests.

### New integration implementation

- Create `scripts/cell_type_specific_expression/scatter_contract.R` for inventory validation and scatter metadata.
- Create `scripts/cell_type_specific_expression/build_qtl_manifest.R` for final manifest validation and writing.
- Create `workflows/cell_type_specific_expression/tasks/integration.wdl` for validation and manifest tasks.
- Create `workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl` for the public end-to-end workflow.

### Existing reusable eQTL implementation

- Modify `workflows/expression/prepare_eQTL.wdl` to accept runtime policy and image inputs without changing its default single-matrix behavior.
- Modify `workflows/common/calculate_phenotypePCs.wdl`, `workflows/common/MergeCovariates.wdl`, and `workflows/common/ResidualizePhenotypes.wdl` to accept those runtime inputs and emit required logs.

### Container, CI, and docs

- Create `envs/CellTypeSpecificExpression/Dockerfile` and `envs/CellTypeSpecificExpression/environment.yml`.
- Create `.github/workflows/cell-type-specific-expression-ci.yml`.
- Modify `.dockstore.yml`, `README.md`, `workflows/README.md`, and create `docs/cell-type-specific-expression.md`.

---

### Task 1: Migrate the reviewed deconvolution R implementation

**Files:**
- Create: `scripts/cell_type_specific_expression/R/bed_outputs.R`
- Create: `scripts/cell_type_specific_expression/R/constants.R`
- Create: `scripts/cell_type_specific_expression/R/dtangle_stage.R`
- Create: `scripts/cell_type_specific_expression/R/expression.R`
- Create: `scripts/cell_type_specific_expression/R/expression_bed.R`
- Create: `scripts/cell_type_specific_expression/R/io.R`
- Create: `scripts/cell_type_specific_expression/R/proportions.R`
- Create: `scripts/cell_type_specific_expression/R/qc.R`
- Create: `scripts/cell_type_specific_expression/R/tca_stage.R`
- Create: `scripts/cell_type_specific_expression/R/workflow_validation.R`
- Create: `scripts/cell_type_specific_expression/bootstrap.R`
- Create: `scripts/cell_type_specific_expression/build_deconvolution_manifest.R`
- Create: `scripts/cell_type_specific_expression/export_tca_beds.R`
- Create: `scripts/cell_type_specific_expression/fit_tca.R`
- Create: `scripts/cell_type_specific_expression/process_proportions.R`
- Create: `scripts/cell_type_specific_expression/run_dtangle.R`
- Create: `scripts/cell_type_specific_expression/validate_proportion_mode.R`
- Create: `scripts/cell_type_specific_expression/lint_r.R`
- Create: `tests/cell_type_specific_expression/testthat.R`
- Create: `tests/cell_type_specific_expression/testthat/helper-load.R`
- Create: `tests/cell_type_specific_expression/testthat/test-bed-outputs.R`
- Create: `tests/cell_type_specific_expression/testthat/test-dtangle-stage.R`
- Create: `tests/cell_type_specific_expression/testthat/test-dtangle.R`
- Create: `tests/cell_type_specific_expression/testthat/test-expression.R`
- Create: `tests/cell_type_specific_expression/testthat/test-io.R`
- Create: `tests/cell_type_specific_expression/testthat/test-proportions.R`
- Create: `tests/cell_type_specific_expression/testthat/test-tca-stage.R`
- Create: `tests/cell_type_specific_expression/testthat/test-workflow-validation.R`

**Interfaces:**
- Consumes: the exact reviewed source files from `/Users/evinmpadhi/Documents/trans sqtls/CellTypeDeconvolution/.worktrees/dtangle-tca-pipeline` at `9a98d2e`.
- Produces: the same R functions and CLI contracts under `/opt/prepare_qtl/scripts/cell_type_specific_expression`, before the pseudocount extension in Task 2.

- [ ] **Step 1: Add a failing migration-boundary test**

Create `tests/cell_type_specific_expression/testthat/test-migration-boundary.R`:

```r
testthat::test_that("the migrated modules use only their canonical directory", {
  root <- normalizePath(testthat::test_path("..", "..", ".."))
  bootstrap <- readLines(
    file.path(root, "scripts", "cell_type_specific_expression", "bootstrap.R"),
    warn = FALSE
  )
  testthat::expect_match(paste(bootstrap, collapse = "\n"),
    "scripts/cell_type_specific_expression/R", fixed = TRUE)
  testthat::expect_false(any(grepl("/opt/celltype", bootstrap, fixed = TRUE)))
})
```

- [ ] **Step 2: Run the test to verify it fails**

Run:

```bash
Rscript -e 'testthat::test_file("tests/cell_type_specific_expression/testthat/test-migration-boundary.R")'
```

Expected: FAIL because the canonical migrated bootstrap file does not exist.

- [ ] **Step 3: Add the reviewed source with path-only adaptations**

Use the contents at commit `9a98d2e` as the source of truth. Map `R/<name>.R` to
`scripts/cell_type_specific_expression/R/<name>.R`. Map each root `scripts/*.R`
CLI to `scripts/cell_type_specific_expression/`, and rename only
`build_manifest.R` to `build_deconvolution_manifest.R` to avoid collision with
the new QTL manifest.

Set the migrated bootstrap root as:

```r
script_root <- Sys.getenv(
  "CELL_TYPE_SPECIFIC_EXPRESSION_ROOT",
  unset = "/opt/prepare_qtl/scripts/cell_type_specific_expression"
)
module_root <- file.path(script_root, "R")
```

For repository tests, `helper-load.R` must set
`CELL_TYPE_SPECIFIC_EXPRESSION_ROOT` to the repository script directory before
it sources the modules.

- [ ] **Step 4: Migrate the tests and update only canonical paths**

Bring over the reviewed R-function tests in the eight files listed above. Keep
all scientific assertions and change source paths to the new directories. WDL,
documentation, and manifest-boundary tests belong to later tasks. Do not bring
over historical tests for deleted shards, HDF5, or expression preparation.
Migrate `tools/lint_r.R` to the listed canonical script path and restrict it to
the migrated R and CLI directories.

- [ ] **Step 5: Run the migrated R suite**

Run:

```bash
Rscript tests/cell_type_specific_expression/testthat.R
```

Expected: PASS with only explicit skips for unavailable local `TCA` 1.2.1 or
`dtangle` 2.0.10.

- [ ] **Step 6: Check the migration diff**

Run:

```bash
git diff --check
rg -n "/opt/celltype|CellTypeDeconvolution" scripts/cell_type_specific_expression tests/cell_type_specific_expression
```

Expected: no obsolete runtime paths. A provenance comment that names the old
repository is permitted only in migration documentation, not runtime source.

- [ ] **Step 7: Commit**

```bash
git add scripts/cell_type_specific_expression tests/cell_type_specific_expression
git commit -m "refactor: migrate cell type deconvolution modules"
```

---

### Task 2: Add the shared log2 pseudocount contract

**Files:**
- Modify: `scripts/cell_type_specific_expression/R/expression.R`
- Modify: `scripts/cell_type_specific_expression/R/expression_bed.R`
- Modify: `scripts/cell_type_specific_expression/run_dtangle.R`
- Modify: `scripts/cell_type_specific_expression/fit_tca.R`
- Modify: `scripts/cell_type_specific_expression/export_tca_beds.R`
- Modify: `scripts/cell_type_specific_expression/build_deconvolution_manifest.R`
- Test: `tests/cell_type_specific_expression/testthat/test-expression.R`
- Test: `tests/cell_type_specific_expression/testthat/test-bed-outputs.R`

**Interfaces:**
- Consumes: the migrated direct-CPM readers and `make_tca_expression()` / `make_dtangle_expression()` from Task 1.
- Produces: `validate_log2_pseudocount(x) -> numeric(1)`, `read_expression_bed(path, log2_pseudocount = 0)`, `make_tca_expression(expression, log2_pseudocount = 0)`, and `make_dtangle_expression(expression, annotation, log2_pseudocount = 0)`.

- [ ] **Step 1: Write failing pseudocount tests**

Add tests that require:

```r
testthat::expect_equal(validate_log2_pseudocount(0), 0)
testthat::expect_equal(validate_log2_pseudocount(0.5), 0.5)
testthat::expect_error(validate_log2_pseudocount(-1), "non-negative")
testthat::expect_error(validate_log2_pseudocount(Inf), "finite")

expression <- list(
  coordinates = tibble::tibble(
    `#chr` = "chr1", start = 0L, end = 1L, gene_id = "g1"
  ),
  cpm = matrix(c(0, 3), nrow = 1,
    dimnames = list("g1", c("s1", "s2")))
)
testthat::expect_error(make_tca_expression(expression, 0), "strictly positive")
testthat::expect_equal(
  unname(make_tca_expression(expression, 1)),
  unname(log2(expression$cpm + 1))
)
```

Add a duplicate-symbol dtangle test that proves the pseudocount is added after
linear-CPM aggregation:

```r
result <- make_dtangle_expression(expression, annotation, 1)
testthat::expect_equal(result$log_expression["GENE1", ],
  log2(colSums(expression$cpm) + 1))
```

- [ ] **Step 2: Run focused tests to verify they fail**

Run:

```bash
Rscript tests/cell_type_specific_expression/testthat.R \
  tests/cell_type_specific_expression/testthat/test-expression.R
```

Expected: FAIL because the functions do not accept or validate a pseudocount.

- [ ] **Step 3: Implement the minimal shared validation and transform**

Use:

```r
validate_log2_pseudocount <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop("log2_pseudocount must be one finite numeric value", call. = FALSE)
  }
  if (x < 0) {
    stop("log2_pseudocount must be non-negative", call. = FALSE)
  }
  as.numeric(x)
}
```

Validate the BED matrix as nonnegative. Reject zero values when the validated
pseudocount is zero. Transform the TCA matrix with
`log2(cpm + log2_pseudocount)`. Aggregate duplicate symbols in linear CPM and
then transform the symbol matrix with the same expression.

- [ ] **Step 4: Pass the option through every expression reader**

Add this `optparse` option to dtangle, TCA fitting, and TCA export:

```r
optparse::make_option(
  "--log2-pseudocount",
  dest = "log2_pseudocount",
  type = "double",
  default = 0,
  help = "Non-negative pseudocount added before the one log2 transform."
)
```

Pass it to every `read_expression_bed()` call and to the matching expression
view constructor. Record `log2_pseudocount` in logs and effective parameters.

- [ ] **Step 5: Run the focused and full migrated suites**

```bash
Rscript tests/cell_type_specific_expression/testthat.R \
  tests/cell_type_specific_expression/testthat/test-expression.R \
  tests/cell_type_specific_expression/testthat/test-bed-outputs.R
Rscript tests/cell_type_specific_expression/testthat.R
```

Expected: PASS with only the recorded local dependency skips.

- [ ] **Step 6: Commit**

```bash
git add scripts/cell_type_specific_expression tests/cell_type_specific_expression
git commit -m "feat: add a configurable deconvolution pseudocount"
```

---

### Task 3: Migrate the standalone WDL and pinned container

**Files:**
- Create: `workflows/cell_type_specific_expression/deconvolution.wdl`
- Create: `workflows/cell_type_specific_expression/tasks/dtangle.wdl`
- Create: `workflows/cell_type_specific_expression/tasks/proportions.wdl`
- Create: `workflows/cell_type_specific_expression/tasks/tca.wdl`
- Create: `workflows/cell_type_specific_expression/tasks/qc.wdl`
- Create: `envs/CellTypeSpecificExpression/Dockerfile`
- Create: `envs/CellTypeSpecificExpression/environment.yml`
- Create: `scripts/check_wdl_logging.py`
- Test: `tests/cell_type_specific_expression/testthat/test-wdl-contract.R`

**Interfaces:**
- Consumes: migrated CLI paths and pseudocount option from Tasks 1-2.
- Produces: public `workflow CellTypeDeconvolution` with the reviewed output contract and `Float log2_pseudocount = 0.0`.

- [ ] **Step 1: Write failing WDL contract tests**

Require all five WDL files to start with `version 1.0`, require the top-level
input below, and require all three expression-consuming calls to receive it:

```wdl
Float log2_pseudocount = 0.0
```

Also require migrated runtime paths to begin with
`/opt/prepare_qtl/scripts/cell_type_specific_expression/` and reject
`/opt/celltype`.

- [ ] **Step 2: Run the WDL tests to verify they fail**

```bash
Rscript tests/cell_type_specific_expression/testthat.R \
  tests/cell_type_specific_expression/testthat/test-wdl-contract.R
```

Expected: FAIL because the migrated WDL files and container do not exist.

- [ ] **Step 3: Add the standalone WDL sources**

Start from the exact reviewed WDL sources at `9a98d2e`. Update import paths,
workflow name, CLI paths, and default image only. Add the pseudocount input to
`RunDtangle`, `FitTca`, and `ExportTcaBeds`, and pass
`--log2-pseudocount '~{log2_pseudocount}'` in each command.

Use this default:

```wdl
String deconvolution_docker_image = "ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression:main"
```

Retain one `preemptible_attempts = 2` input and one `max_retries = 2` input,
passed to every deconvolution task. Retain task-specific resources.

Migrate the reviewed WDL logging checker from `CellTypeDeconvolution` commit
`9a98d2e` to `scripts/check_wdl_logging.py`. Use it for every new or modified
WDL source.

- [ ] **Step 4: Add the pinned environment and Dockerfile**

Use:

```yaml
name: cell-type-specific-expression
channels:
  - conda-forge
  - bioconda
  - nodefaults
dependencies:
  - r-base=4.5.3
  - r-dtangle=2.0.10
  - r-tca=1.2.1
  - r-tidyverse=2.0.0
  - r-optparse=1.8.2
  - r-jsonlite=2.0.0
  - r-digest=0.6.39
  - r-testthat=3.3.1
  - r-lintr=3.4.0
  - bioconductor-limma=3.66.0
```

Use `mambaorg/micromamba:2.9.0-ubuntu22.04`, install the environment into
`base`, copy only `scripts/cell_type_specific_expression`, and verify dtangle
and TCA versions during the build.

- [ ] **Step 5: Run static checks without building Docker**

```bash
miniwdl check workflows/cell_type_specific_expression/deconvolution.wdl
python3 scripts/validate_dockstore.py
git diff --check
```

Expected: MiniWDL passes. Dockstore still passes for existing descriptors. Do
not build the new image locally.

- [ ] **Step 6: Commit**

```bash
git add workflows/cell_type_specific_expression envs/CellTypeSpecificExpression \
  tests/cell_type_specific_expression
git commit -m "feat: migrate the standalone deconvolution workflow"
```

---

### Task 4: Make eQTL preparation reusable inside the scatter

**Files:**
- Modify: `workflows/expression/prepare_eQTL.wdl`
- Modify: `workflows/common/calculate_phenotypePCs.wdl`
- Modify: `workflows/common/MergeCovariates.wdl`
- Modify: `workflows/common/ResidualizePhenotypes.wdl`
- Test: `tests/testthat/test_expression_log2_cpm_wdl.R`
- Test: `tests/testthat/test_phenotype_pc_workflow_outputs.R`
- Create: `tests/testthat/test_eqtl_runtime_contract.R`

**Interfaces:**
- Consumes: existing `eQTLPrepareData` inputs and outputs.
- Produces: backward-compatible optional runtime inputs `DockerImage`, `preemptible_attempts`, and `max_retries` on the eQTL and imported common workflows.

- [ ] **Step 1: Write failing runtime and logging tests**

Require these defaults on `eQTLPrepareData` and every invoked task:

```wdl
String DockerImage = "ghcr.io/aou-multiomics-analysis/prepare_qtl:main"
Int preemptible_attempts = 2
Int max_retries = 2
```

Require each modified command block to contain `stage=`, `start_time=`,
`completion_time=`, `dimensions=`, and `outputs=`. Require each runtime block to
use `docker: DockerImage`, `preemptible: preemptible_attempts`, and
`maxRetries: max_retries`.

- [ ] **Step 2: Run focused tests to verify they fail**

```bash
Rscript tests/testthat/test_eqtl_runtime_contract.R
```

Expected: FAIL because common tasks use hard-coded images and omit retry policy
and complete logging.

- [ ] **Step 3: Add backward-compatible runtime inputs**

Add the three inputs with defaults to each workflow and task. Pass them through
all calls in `eQTLPrepareData`, including residualization even though the new
end-to-end workflow sets `ResidualizeNormalizedInputs = false`.

Do not change the mutually exclusive `CountGCT` / `Log2CpmBed` behavior or any
existing output type.

- [ ] **Step 4: Add complete WDL command logging**

Use this command structure in each modified task:

```bash
set -euo pipefail
stage="<task_stage>"
printf 'stage=%s start_time=%s dimensions=pending outputs=%s\n' \
  "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "<declared outputs>"
<existing command>
printf 'stage=%s completion_time=%s dimensions=complete outputs=%s\n' \
  "$stage" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "<declared outputs>"
```

- [ ] **Step 5: Run existing and new tests**

```bash
Rscript tests/testthat/test_expression_log2_cpm_wdl.R
Rscript tests/testthat/test_phenotype_pc_workflow_outputs.R
Rscript tests/testthat/test_eqtl_runtime_contract.R
miniwdl check workflows/expression/prepare_eQTL.wdl
```

Expected: PASS. The container-dependent expression integration remains a
GitHub Actions check because local `biomaRt` is unavailable.

- [ ] **Step 6: Commit**

```bash
git add workflows/expression workflows/common tests/testthat
git commit -m "refactor: expose reusable eQTL runtime policy"
```

---

### Task 5: Build the validated scatter contract

**Files:**
- Create: `scripts/cell_type_specific_expression/R/integration.R`
- Create: `scripts/cell_type_specific_expression/prepare_scatter_inputs.R`
- Create: `workflows/cell_type_specific_expression/tasks/integration.wdl`
- Create: `tests/cell_type_specific_expression/testthat/test-integration.R`

**Interfaces:**
- Consumes: `File cell_type_bed_inventory` and `Array[File] cell_type_beds`.
- Produces: aligned `Array[String] cell_types`, `Array[String] cell_type_slugs`, `Array[File] expression_beds`, and `Array[String] output_prefixes`.

- [ ] **Step 1: Write failing inventory validation tests**

Define and test:

```r
prepare_scatter_contract <- function(inventory, bed_paths, output_prefix)
```

The result must contain `cell_type`, `cell_type_slug`, `expression_bed`, and
`output_prefix`. Add failures for missing columns, duplicate cell types,
duplicate slugs, duplicate paths, wrong scale, nonpositive dimensions,
basename mismatch, count mismatch, and task-local directory components in the
public path column.

Require exact prefix construction:

```r
testthat::expect_equal(
  contract$output_prefix,
  c("cohort.cd4_t_cells", "cohort.monocytes")
)
```

- [ ] **Step 2: Run focused tests to verify they fail**

```bash
Rscript tests/cell_type_specific_expression/testthat.R \
  tests/cell_type_specific_expression/testthat/test-integration.R
```

Expected: FAIL because `prepare_scatter_contract()` does not exist.

- [ ] **Step 3: Implement the pure R contract**

Validate this exact inventory schema:

```r
c("logical_name", "path", "n_genes", "n_samples", "scale", "cell_group", "slug")
```

Match `basename(bed_paths)` to `inventory$path` without sorting. Return a tibble
in inventory order and preserve display names separately from slugs.

- [ ] **Step 4: Add the CLI and WDL task**

The CLI writes four newline-delimited files. The task outputs them as arrays:

```wdl
Array[String] cell_types = read_lines("scatter/cell_types.txt")
Array[String] cell_type_slugs = read_lines("scatter/cell_type_slugs.txt")
Array[File] expression_beds = cell_type_beds
Array[String] output_prefixes = read_lines("scatter/output_prefixes.txt")
```

The task command must log validated cell count and output paths. It must not
copy the BEDs.

- [ ] **Step 5: Run R and WDL checks**

```bash
Rscript tests/cell_type_specific_expression/testthat.R \
  tests/cell_type_specific_expression/testthat/test-integration.R
miniwdl check workflows/cell_type_specific_expression/tasks/integration.wdl
```

Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add scripts/cell_type_specific_expression workflows/cell_type_specific_expression/tasks/integration.wdl tests/cell_type_specific_expression
git commit -m "feat: validate cell type scatter inputs"
```

---

### Task 6: Add the end-to-end scattered WDL

**Files:**
- Create: `workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl`
- Create: `tests/cell_type_specific_expression/testthat/test-end-to-end-wdl.R`

**Interfaces:**
- Consumes: standalone `CellTypeDeconvolution`, `PrepareScatterInputs`, and imported `expression/prepare_eQTL.wdl`.
- Produces: one `PrepareCellTypeEqtl` subworkflow call per retained cell type and aligned arrays of every requested QTL-ready output.

- [ ] **Step 1: Write a failing orchestration contract test**

Require imports for the standalone workflow, integration tasks, and existing
eQTL workflow. Require `File AdditionalCovariates`, not `File?`. Require:

```wdl
scatter (index in range(length(PrepareScatterInputs.cell_types)))
```

Require each eQTL call to receive its aligned BED and prefix, the shared sample
list and covariates, `ResidualizeNormalizedInputs = false`, the standard QTL
image, and the global retry policy.

- [ ] **Step 2: Run the test to verify it fails**

```bash
Rscript tests/cell_type_specific_expression/testthat.R \
  tests/cell_type_specific_expression/testthat/test-end-to-end-wdl.R
```

Expected: FAIL because the public orchestrator does not exist.

- [ ] **Step 3: Implement the orchestration**

Inside the scatter, call:

```wdl
call eqtl.eQTLPrepareData as PrepareCellTypeEqtl {
  input:
    OutputPrefix = PrepareScatterInputs.output_prefixes[index],
    Log2CpmBed = PrepareScatterInputs.expression_beds[index],
    SampleList = SampleList,
    AdditionalCovariates = AdditionalCovariates,
    ResidualizeNormalizedInputs = false,
    DockerImage = qtl_docker_image,
    preemptible_attempts = preemptible_attempts,
    max_retries = max_retries,
    memory = eqtl_memory,
    disk_space = eqtl_disk_gb,
    num_threads = eqtl_cpu
}
```

Convert the optional merged-covariate outputs to required files inside the
scatter:

```wdl
File int_merged_covariates = select_first([PrepareCellTypeEqtl.IntQtlCovariates])
File scaled_merged_covariates = select_first([PrepareCellTypeEqtl.ScaledQtlCovariates])
```

Expose aligned arrays for INT/scale BEDs, selected/all PCs, merged covariates,
and outlier reports. Do not expose residualized BED arrays.

- [ ] **Step 4: Run MiniWDL and contract tests**

```bash
miniwdl check workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl
Rscript tests/cell_type_specific_expression/testthat.R \
  tests/cell_type_specific_expression/testthat/test-end-to-end-wdl.R
```

Expected: PASS and resolve all nested imports under WDL 1.0.

- [ ] **Step 5: Commit**

```bash
git add workflows/cell_type_specific_expression tests/cell_type_specific_expression
git commit -m "feat: scatter eQTL preparation by cell type"
```

---

### Task 7: Create the per-cell-type QTL manifest

**Files:**
- Modify: `scripts/cell_type_specific_expression/R/integration.R`
- Create: `scripts/cell_type_specific_expression/build_qtl_manifest.R`
- Modify: `workflows/cell_type_specific_expression/tasks/integration.wdl`
- Modify: `workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl`
- Test: `tests/cell_type_specific_expression/testthat/test-qtl-manifest.R`

**Interfaces:**
- Consumes: cell names/slugs and 10 aligned `Array[File]` values from Task 6.
- Produces: `cell_type_qtl_manifest.tsv` with one row per retained cell type and stable basenames.

- [ ] **Step 1: Write failing manifest tests**

Require this exact schema and order:

```r
c(
  "cell_type", "cell_type_slug", "int_bed", "scaled_bed",
  "int_phenotype_pcs", "int_phenotype_pcs_all",
  "scaled_phenotype_pcs", "scaled_phenotype_pcs_all",
  "int_merged_covariates", "scaled_merged_covariates",
  "int_connectivity_outliers", "scaled_connectivity_outliers"
)
```

Test mismatched array lengths, duplicate basenames within one file category,
empty names/slugs, missing files, task-local paths in the written manifest, and
preservation of scatter order.

- [ ] **Step 2: Run tests to verify they fail**

```bash
Rscript tests/cell_type_specific_expression/testthat.R \
  tests/cell_type_specific_expression/testthat/test-qtl-manifest.R
```

Expected: FAIL because `build_cell_type_qtl_manifest()` does not exist.

- [ ] **Step 3: Implement the pure manifest function**

Define:

```r
build_cell_type_qtl_manifest <- function(
    cell_types,
    cell_type_slugs,
    int_beds,
    scaled_beds,
    int_pcs,
    int_pcs_all,
    scaled_pcs,
    scaled_pcs_all,
    int_covariates,
    scaled_covariates,
    int_outliers,
    scaled_outliers) {
  file_lists <- list(
    int_bed = int_beds,
    scaled_bed = scaled_beds,
    int_phenotype_pcs = int_pcs,
    int_phenotype_pcs_all = int_pcs_all,
    scaled_phenotype_pcs = scaled_pcs,
    scaled_phenotype_pcs_all = scaled_pcs_all,
    int_merged_covariates = int_covariates,
    scaled_merged_covariates = scaled_covariates,
    int_connectivity_outliers = int_outliers,
    scaled_connectivity_outliers = scaled_outliers
  )
  expected_length <- length(cell_types)
  if (expected_length == 0L || length(cell_type_slugs) != expected_length ||
      any(lengths(file_lists) != expected_length)) {
    stop("Cell-type manifest arrays must have one equal nonzero length",
      call. = FALSE)
  }
  if (anyNA(cell_types) || any(!nzchar(cell_types)) ||
      anyDuplicated(cell_types) > 0L || anyNA(cell_type_slugs) ||
      any(!nzchar(cell_type_slugs)) || anyDuplicated(cell_type_slugs) > 0L) {
    stop("Cell types and slugs must be nonempty and unique", call. = FALSE)
  }
  if (any(!file.exists(unlist(file_lists, use.names = FALSE)))) {
    stop("Every cell-type manifest file must exist", call. = FALSE)
  }
  if (any(vapply(file_lists, anyDuplicated, integer(1)) > 0L)) {
    stop("Manifest file paths must be unique within each output category",
      call. = FALSE)
  }
  tibble::tibble(
    cell_type = cell_types,
    cell_type_slug = cell_type_slugs
  ) |>
    dplyr::bind_cols(purrr::map_dfc(file_lists, basename))
}
```

The returned tibble must apply `basename()` to every file field. It must never
sort rows or infer identity from output filenames.

- [ ] **Step 4: Add the CLI and BuildQtlManifest task**

Pass each WDL file array through a newline-delimited argument file to avoid
shell-length and quoting errors. The task must validate the array lengths in R,
write `outputs/cell_type_qtl_manifest.tsv`, and log cell count, manifest path,
and completion.

- [ ] **Step 5: Wire the manifest and authoritative arrays into public outputs**

The top-level workflow outputs the manifest plus one ordered `Array[File]` for
every file column. Also expose the standalone deconvolution proportions, TCA
weights/model, filter report, original cell-type BEDs/inventory,
reconstruction QC, effective parameters, logs, and deconvolution manifest.

- [ ] **Step 6: Run focused, full, and static tests**

```bash
Rscript tests/cell_type_specific_expression/testthat.R \
  tests/cell_type_specific_expression/testthat/test-qtl-manifest.R \
  tests/cell_type_specific_expression/testthat/test-end-to-end-wdl.R
Rscript tests/cell_type_specific_expression/testthat.R
miniwdl check workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl
git diff --check
```

Expected: PASS with only the explicit local TCA/dtangle skips.

- [ ] **Step 7: Commit**

```bash
git add scripts/cell_type_specific_expression workflows/cell_type_specific_expression tests/cell_type_specific_expression
git commit -m "feat: build the cell type QTL manifest"
```

---

### Task 8: Add deterministic fixtures and GitHub Actions smoke tests

**Files:**
- Modify: `scripts/cell_type_specific_expression/generate_synthetic_fixture.R`
- Create: `tests/cell_type_specific_expression/fixtures/synthetic_expression.bed`
- Create: `tests/cell_type_specific_expression/fixtures/synthetic.gtf`
- Create: `tests/cell_type_specific_expression/fixtures/synthetic_signature.tsv`
- Create: `tests/cell_type_specific_expression/fixtures/precomputed_proportions.tsv`
- Create: `tests/cell_type_specific_expression/fixtures/samples.tsv`
- Create: `tests/cell_type_specific_expression/fixtures/additional_covariates.tsv`
- Create: `tests/cell_type_specific_expression/fixtures/expected_groups.txt`
- Create: `tests/cell_type_specific_expression/fixtures/dtangle-e2e.inputs.json`
- Create: `tests/cell_type_specific_expression/fixtures/precomputed-e2e.inputs.json`
- Create: `tests/cell_type_specific_expression/smoke/assert_deconvolution_outputs.R`
- Create: `tests/cell_type_specific_expression/smoke/assert_qtl_outputs.R`
- Create: `.github/workflows/cell-type-specific-expression-ci.yml`
- Modify: `.github/workflows/docker-image.yml`

**Interfaces:**
- Consumes: both public WDLs and both container definitions.
- Produces: deterministic dtangle and precomputed-proportion fixture inputs and GitHub evidence for the full scatter.

- [ ] **Step 1: Extend the fixture contract before changing the generator**

Add tests that require two end-to-end input JSON files. One must use LM22 and
one must use precomputed proportions. Use `log2_pseudocount = 0.0` in one and a
positive value with at least one zero CPM in the other. Both must require
`AdditionalCovariates`, omit residualization, and point to local CI image tags.

- [ ] **Step 2: Run fixture tests to verify they fail**

```bash
Rscript tests/cell_type_specific_expression/testthat.R \
  tests/cell_type_specific_expression/testthat/test-workflow-validation.R
```

Expected: FAIL because integrated fixtures do not exist.

- [ ] **Step 3: Generate deterministic fixtures**

Extend the migrated generator to write:

```text
dtangle-e2e.inputs.json
precomputed-e2e.inputs.json
samples.tsv
additional_covariates.tsv
expected_groups.txt
```

Generate twice into separate temporary directories and require `diff -rq` to
produce no output. Compare one generated directory to the checked-in fixtures.

- [ ] **Step 4: Add complete smoke assertions**

The QTL assertion must verify one row per expected group, the exact 12-column
manifest schema, stable basenames, matching authoritative arrays, readable
INT/scaled BEDs, selected/all PC files, merged covariates, outlier reports, and
important upstream deconvolution outputs. It must reject residualized manifest
columns and confirm that the positive-pseudocount mode accepts the zero input.

- [ ] **Step 5: Add GitHub Actions CI**

The workflow must:

1. Install `miniwdl==1.15.0`.
2. Run MiniWDL checks and the WDL logging checker.
3. Build the standard image from `envs/PhenotypePCs/Dockerfile` as
   `prepare-qtl:test`.
4. Build the dedicated image from
   `envs/CellTypeSpecificExpression/Dockerfile` as
   `cell-type-specific-expression:test`.
5. Run the migrated R suite inside the dedicated image.
6. Run the existing pre-normalized expression test inside the standard image.
7. Run both end-to-end WDL fixtures.
8. Assert both outputs inside the appropriate image.
9. Upload logs, `error.json`, `outputs.json`, manifests, and QC files only when
   a step fails.

Use Docker only in GitHub Actions. Do not run either Docker build locally.

Add a second build-and-push job to `.github/workflows/docker-image.yml` for
`envs/CellTypeSpecificExpression/Dockerfile`. Publish it as
`ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression`
with the same event-specific metadata tags as the standard image. The WDL
default `:main` tag must therefore exist after merge to `main`.

- [ ] **Step 6: Run all local non-container checks**

```bash
Rscript tests/cell_type_specific_expression/testthat.R
Rscript tests/testthat/test_expression_log2_cpm_wdl.R
Rscript tests/testthat/test_phenotype_pc_workflow_outputs.R
miniwdl check workflows/cell_type_specific_expression/deconvolution.wdl
miniwdl check workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl
python3 scripts/validate_dockstore.py
git diff --check
```

Expected: PASS with only explicit local package skips. Do not claim the
container smoke passes until GitHub Actions runs it.

- [ ] **Step 7: Commit**

```bash
git add .github/workflows/cell-type-specific-expression-ci.yml \
  .github/workflows/docker-image.yml \
  scripts/cell_type_specific_expression tests/cell_type_specific_expression
git commit -m "ci: smoke test cell type specific eQTL preparation"
```

---

### Task 9: Register and document the canonical workflows

**Files:**
- Modify: `.dockstore.yml`
- Modify: `README.md`
- Modify: `workflows/README.md`
- Create: `docs/cell-type-specific-expression.md`
- Create: `tests/testthat/test_cell_type_specific_documentation.R`

**Interfaces:**
- Consumes: final public input/output names from Tasks 3 and 6-7.
- Produces: Dockstore registrations and user guidance for both public workflows.

- [ ] **Step 1: Write failing registration and documentation tests**

Require `.dockstore.yml` to register:

```text
/workflows/cell_type_specific_expression/deconvolution.wdl
/workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl
```

Require the docs to state the input scale, pseudocount rules, one-time log2
transform, no second QTL log transform, independent connectivity filtering,
required additional covariates, manifest columns, basename rule, authoritative
arrays, and two-branch merge strategy.

- [ ] **Step 2: Run tests to verify they fail**

```bash
Rscript tests/testthat/test_cell_type_specific_documentation.R
```

Expected: FAIL because registrations and guidance are absent.

- [ ] **Step 3: Add Dockstore entries and documentation**

Add distinct workflow names `cell_type_deconvolution` and
`prepare_cell_type_eQTL`. Document example inputs for both proportion modes and
the exact 12-column manifest. State that `prepare_QTL` is canonical and the old
repository is deprecated after this integration reaches `main`.

- [ ] **Step 4: Run the documentation and descriptor checks**

```bash
Rscript tests/testthat/test_cell_type_specific_documentation.R
python3 scripts/validate_dockstore.py
miniwdl check workflows/cell_type_specific_expression/deconvolution.wdl
miniwdl check workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl
git diff --check
```

Expected: PASS and list both new descriptors.

- [ ] **Step 5: Commit**

```bash
git add .dockstore.yml README.md workflows/README.md \
  docs/cell-type-specific-expression.md tests/testthat/test_cell_type_specific_documentation.R
git commit -m "docs: publish cell type specific eQTL workflows"
```

---

## Final verification and integration sequence

After all nine tasks and per-task reviews:

1. Run one whole-branch code review against `feat/log2-cpm-bed-input`.
2. Permit one final fix wave and one scoped re-review.
3. Run fresh local verification:

```bash
Rscript tests/cell_type_specific_expression/testthat.R
Rscript tests/testthat/test_expression_log2_cpm_wdl.R
Rscript tests/testthat/test_phenotype_pc_workflow_outputs.R
Rscript tests/testthat/test_eqtl_runtime_contract.R
Rscript tests/testthat/test_cell_type_specific_documentation.R
Rscript tests/test_prepare_methylation.R
miniwdl check workflows/cell_type_specific_expression/deconvolution.wdl
miniwdl check workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl
python3 scripts/validate_dockstore.py
git diff --check
git status --short
```

4. Push `codex/cell-type-specific-expression`.
5. Open the cell-type-specific pull request against `main` and state that it
   depends on the pre-normalized-expression pull request. The branch contains
   that dependency's commits until they reach `main`; targeting `main` ensures
   normal pull-request CI runs.
6. Monitor GitHub Actions until both end-to-end smoke modes pass.
7. Merge `feat/log2-cpm-bed-input` into `main` first.
8. Confirm that the cell-type-specific pull-request diff shrinks after the
   dependency merges. Rebase onto `main` only if GitHub does not remove the
   already merged commits cleanly. Rerun all checks and merge after review.
9. Only after the canonical workflow is present on `prepare_QTL/main`, add a
   deprecation notice to the old `CellTypeDeconvolution` README that links to
   the two canonical workflow paths. Do not delete the old repository.

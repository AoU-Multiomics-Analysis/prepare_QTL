end_to_end_wdl_path <- function() {
  testthat::test_path(
    "..",
    "..",
    "..",
    "workflows",
    "cell_type_specific_expression",
    "prepare_cell_type_eQTL.wdl"
  )
}

read_end_to_end_wdl <- function() {
  path <- end_to_end_wdl_path()
  testthat::expect_true(file.exists(path), info = path)
  if (!file.exists(path)) {
    return("")
  }
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

testthat::test_that("end-to-end WDL imports and calls the reusable workflows", {
  text <- read_end_to_end_wdl()
  if (!nzchar(text)) {
    return(invisible())
  }

  testthat::expect_identical(
    readLines(end_to_end_wdl_path(), n = 1L, warn = FALSE),
    "version 1.0"
  )
  testthat::expect_match(
    text,
    'import "deconvolution.wdl" as deconvolution',
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    'import "tasks/integration.wdl" as integration',
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    'import "../expression/prepare_eQTL.wdl" as eqtl',
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "call deconvolution.CellTypeDeconvolution as CellTypeDeconvolution",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "call integration.PrepareScatterInputs as PrepareScatterInputs",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "log2_pseudocount = log2_pseudocount",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "deconvolution_docker_image = deconvolution_docker_image",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "docker_image = deconvolution_docker_image",
    fixed = TRUE
  )
})

testthat::test_that("end-to-end WDL requires QTL covariates and scatters by aligned index", {
  text <- read_end_to_end_wdl()
  if (!nzchar(text)) {
    return(invisible())
  }

  testthat::expect_match(text, "File AdditionalCovariates", fixed = TRUE)
  testthat::expect_false(grepl("File? AdditionalCovariates", text, fixed = TRUE))
  testthat::expect_match(
    text,
    "scatter (index in range(length(PrepareScatterInputs.cell_types)))",
    fixed = TRUE
  )

  assignments <- c(
    "OutputPrefix = PrepareScatterInputs.output_prefixes[index]",
    "Log2CpmBed = PrepareScatterInputs.expression_beds[index]",
    "SampleList = SampleList",
    "AdditionalCovariates = AdditionalCovariates",
    "ResidualizeNormalizedInputs = false",
    "DockerImage = qtl_docker_image",
    "preemptible_attempts = preemptible_attempts",
    "max_retries = max_retries",
    "memory = eqtl_memory",
    "disk_space = eqtl_disk_gb",
    "num_threads = eqtl_cpu"
  )
  purrr::walk(assignments, ~ testthat::expect_match(text, .x, fixed = TRUE))
  testthat::expect_identical(
    stringr::str_count(
      text,
      stringr::fixed("preemptible_attempts = preemptible_attempts")
    ),
    3L
  )
  testthat::expect_identical(
    stringr::str_count(
      text,
      stringr::fixed("max_retries = max_retries")
    ),
    3L
  )

  testthat::expect_match(
    text,
    "File int_merged_covariate = select_first([PrepareCellTypeEqtl.IntQtlCovariates])",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "File scaled_merged_covariate = select_first([PrepareCellTypeEqtl.ScaledQtlCovariates])",
    fixed = TRUE
  )
})

testthat::test_that("end-to-end WDL exposes aligned QTL-ready file arrays", {
  text <- read_end_to_end_wdl()
  if (!nzchar(text)) {
    return(invisible())
  }

  expected_outputs <- c(
    "Array[File] int_beds = PrepareCellTypeEqtl.IntBedFile",
    "Array[File] scaled_beds = PrepareCellTypeEqtl.ScaledBedFile",
    "Array[File] int_phenotype_pcs = PrepareCellTypeEqtl.IntPhenotypePCsOut",
    "Array[File] int_phenotype_pcs_all = PrepareCellTypeEqtl.IntPhenotypePCsAllOut",
    "Array[File] scaled_phenotype_pcs = PrepareCellTypeEqtl.ScaledPhenotypePCsOut",
    "Array[File] scaled_phenotype_pcs_all = PrepareCellTypeEqtl.ScaledPhenotypePCsAllOut",
    "Array[File] int_merged_covariates = int_merged_covariate",
    "Array[File] scaled_merged_covariates = scaled_merged_covariate",
    "Array[File] int_connectivity_outliers = PrepareCellTypeEqtl.IntConnectivityOutliers",
    "Array[File] scaled_connectivity_outliers = PrepareCellTypeEqtl.ScaledConnectivityOutliers"
  )
  purrr::walk(expected_outputs, ~ testthat::expect_match(text, .x, fixed = TRUE))
  testthat::expect_false(grepl("ResidualizedBedFile", text, fixed = TRUE))
})

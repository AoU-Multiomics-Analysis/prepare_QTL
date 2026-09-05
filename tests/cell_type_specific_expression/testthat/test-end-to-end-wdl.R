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

end_to_end_output_declarations <- function(text) {
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  output_start <- which(lines == "  output {")
  testthat::expect_length(output_start, 1L)
  if (length(output_start) != 1L) {
    return(character())
  }
  relative_end <- which(lines[(output_start + 1L):length(lines)] == "  }")
  testthat::expect_true(length(relative_end) > 0L)
  if (length(relative_end) == 0L) {
    return(character())
  }
  output_end <- output_start + relative_end[[1L]]
  trimws(lines[(output_start + 1L):(output_end - 1L)])
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
    "call integration.BuildQtlManifest as BuildQtlManifest",
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
    "CpmBed = PrepareScatterInputs.expression_beds[index]",
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
    4L
  )
  testthat::expect_identical(
    stringr::str_count(
      text,
      stringr::fixed("max_retries = max_retries")
    ),
    4L
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

  output_declarations <- end_to_end_output_declarations(text)
  array_matches <- stringr::str_match(
    output_declarations,
    "^Array\\[File\\]\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*(.+)$"
  )
  qtl_array_names <- array_matches[
    !is.na(array_matches[, 1L]) &
      grepl("PrepareCellTypeEqtl|merged_covariate", array_matches[, 3L]),
    2L
  ]
  expected_qtl_array_names <- c(
    "int_beds", "scaled_beds", "int_phenotype_pcs",
    "int_phenotype_pcs_all", "scaled_phenotype_pcs",
    "scaled_phenotype_pcs_all", "int_merged_covariates",
    "scaled_merged_covariates", "int_connectivity_outliers",
    "scaled_connectivity_outliers"
  )
  testthat::expect_identical(qtl_array_names, expected_qtl_array_names)
  all_array_names <- array_matches[!is.na(array_matches[, 1L]), 2L]
  testthat::expect_identical(
    all_array_names,
    c("cell_type_beds", "filtered_cell_type_beds", "reference_filter_plots",
      expected_qtl_array_names)
  )

  public_output_names <- stringr::str_match(
    output_declarations,
    "^(?:Array\\[File\\]|File\\??)\\s+([A-Za-z_][A-Za-z0-9_]*)\\s*="
  )[, 2L]
  public_output_names <- public_output_names[!is.na(public_output_names)]
  testthat::expect_false(any(grepl(
    "raw|residualized",
    public_output_names,
    ignore.case = TRUE
  )))
})

testthat::test_that("end-to-end WDL builds the manifest from authoritative arrays", {
  text <- read_end_to_end_wdl()
  if (!nzchar(text)) {
    return(invisible())
  }
  call_block <- stringr::str_extract(
    text,
    "call integration[.]BuildQtlManifest as BuildQtlManifest \\{[\\s\\S]*?\\n  \\}"
  )
  testthat::expect_false(is.na(call_block))
  expected_assignments <- c(
    "cell_types = PrepareScatterInputs.cell_types",
    "cell_type_slugs = PrepareScatterInputs.cell_type_slugs",
    "int_beds = PrepareCellTypeEqtl.IntBedFile",
    "scaled_beds = PrepareCellTypeEqtl.ScaledBedFile",
    "int_phenotype_pcs = PrepareCellTypeEqtl.IntPhenotypePCsOut",
    "int_phenotype_pcs_all = PrepareCellTypeEqtl.IntPhenotypePCsAllOut",
    "scaled_phenotype_pcs = PrepareCellTypeEqtl.ScaledPhenotypePCsOut",
    "scaled_phenotype_pcs_all = PrepareCellTypeEqtl.ScaledPhenotypePCsAllOut",
    "int_merged_covariates = int_merged_covariate",
    "scaled_merged_covariates = scaled_merged_covariate",
    "int_connectivity_outliers = PrepareCellTypeEqtl.IntConnectivityOutliers",
    "scaled_connectivity_outliers = PrepareCellTypeEqtl.ScaledConnectivityOutliers",
    "source_beds = CellTypeDeconvolution.cell_type_beds",
    "filtered_beds = CellTypeDeconvolution.filtered_cell_type_beds",
    "negative_summary = CellTypeDeconvolution.negative_expression_summary",
    "gene_comparison = CellTypeDeconvolution.reference_gene_comparison",
    "filter_metrics = CellTypeDeconvolution.reference_filter_metrics"
  )
  purrr::walk(
    expected_assignments,
    ~ testthat::expect_match(call_block, .x, fixed = TRUE)
  )

  output_declarations <- end_to_end_output_declarations(text)
  testthat::expect_true(
    "File cell_type_qtl_manifest = BuildQtlManifest.manifest" %in%
      output_declarations
  )
  testthat::expect_true(
    "File cell_type_qtl_manifest_log = BuildQtlManifest.log" %in%
      output_declarations
  )
})

testthat::test_that("end-to-end WDL exposes the full deconvolution provenance", {
  text <- read_end_to_end_wdl()
  if (!nzchar(text)) {
    return(invisible())
  }
  output_declarations <- end_to_end_output_declarations(text)
  expected_outputs <- c(
    "File? proportion_mode_validation_log = CellTypeDeconvolution.proportion_mode_validation_log",
    "File? estimated_proportions = CellTypeDeconvolution.estimated_proportions",
    "File? hspe_markers = CellTypeDeconvolution.hspe_markers",
    "File? hspe_metadata = CellTypeDeconvolution.hspe_metadata",
    "File? hspe_overlap_report = CellTypeDeconvolution.hspe_overlap_report",
    "File? transformed_lm22 = CellTypeDeconvolution.transformed_lm22",
    "File? hspe_log = CellTypeDeconvolution.hspe_log",
    "File? proportions_lm22 = CellTypeDeconvolution.proportions_lm22",
    "File? proportions_combined = CellTypeDeconvolution.proportions_combined",
    "File tca_weights = CellTypeDeconvolution.tca_weights",
    "File? cell_group_filter_report = CellTypeDeconvolution.cell_group_filter_report",
    "File? proportions_log = CellTypeDeconvolution.proportions_log",
    "File tca_model = CellTypeDeconvolution.tca_model",
    "File? tca_model_log = CellTypeDeconvolution.tca_model_log",
    "File? tca_excluded_genes = CellTypeDeconvolution.tca_excluded_genes",
    "File? fit_tca_log = CellTypeDeconvolution.fit_tca_log",
    "Array[File] cell_type_beds = CellTypeDeconvolution.cell_type_beds",
    "File cell_type_bed_inventory = CellTypeDeconvolution.cell_type_bed_inventory",
    "File reconstruction_by_sample = CellTypeDeconvolution.reconstruction_by_sample",
    "File qc_summary = CellTypeDeconvolution.qc_summary",
    "File qc_plots = CellTypeDeconvolution.qc_plots",
    "File export_log = CellTypeDeconvolution.export_log",
    "File export_detail_log = CellTypeDeconvolution.export_detail_log",
    "File output_manifest = CellTypeDeconvolution.output_manifest",
    "File output_inventory = CellTypeDeconvolution.output_inventory",
    "File manifest_log = CellTypeDeconvolution.manifest_log",
    "File effective_parameters_file = CellTypeDeconvolution.effective_parameters_file"
  )

  purrr::walk(
    expected_outputs,
    ~ testthat::expect_true(.x %in% output_declarations, info = .x)
  )
  testthat::expect_false(any(grepl("hspe_shared_bulk", output_declarations)))
})

qtl_manifest_columns <- c(
  "entity:cell_type_id", "cell_type", "cell_type_slug", "int_bed", "scaled_bed",
  "int_phenotype_pcs", "int_phenotype_pcs_all",
  "scaled_phenotype_pcs", "scaled_phenotype_pcs_all",
  "int_merged_covariates", "scaled_merged_covariates",
  "int_connectivity_outliers", "scaled_connectivity_outliers", "source_cpm_bed",
  "filtered_cpm_bed", "negative_expression_summary", "reference_gene_comparison",
  "reference_filter_metrics"
)

make_qtl_manifest_inputs <- function() {
  arguments <- c(
    "int_beds", "scaled_beds", "int_pcs", "int_pcs_all",
    "scaled_pcs", "scaled_pcs_all", "int_covariates",
    "scaled_covariates", "int_outliers", "scaled_outliers"
  )
  files <- purrr::map(arguments, function(category) {
    paste0(
      "gs://test-bucket/submission/call-eqtl/shard-", c(1, 0),
      "/outputs/", category, ".tsv"
    )
  })
  names(files) <- arguments
  c(list(
    cell_types = c("Monocytes", "CD4 T cells"),
    cell_type_slugs = c("monocytes", "cd4_t_cells")
  ), files)
}

testthat::test_that("QTL manifest preserves cloud URLs and scatter identity without localization", {
  inputs <- make_qtl_manifest_inputs()
  manifest <- do.call(build_cell_type_qtl_manifest, inputs)
  testthat::expect_identical(names(manifest), qtl_manifest_columns)
  testthat::expect_identical(manifest[["entity:cell_type_id"]], c("monocytes", "cd4_t_cells"))
  testthat::expect_identical(manifest$cell_type, c("Monocytes", "CD4 T cells"))
  testthat::expect_identical(manifest$cell_type_slug, c("monocytes", "cd4_t_cells"))
  testthat::expect_identical(manifest$int_bed, c(
    "gs://test-bucket/submission/call-eqtl/shard-1/outputs/int_beds.tsv",
    "gs://test-bucket/submission/call-eqtl/shard-0/outputs/int_beds.tsv"
  ))
  testthat::expect_identical(unname(as.list(manifest[4:13])), unname(inputs[-(1:2)]))
  testthat::expect_true(all(is.na(manifest[14:18])))
})

testthat::test_that("QTL manifest aligns filter BED metadata by slug and permits shared reports", {
  inputs <- make_qtl_manifest_inputs()
  inputs$source_beds <- c(
    "gs://bucket/raw/b_cells.bed.gz", "gs://bucket/raw/cd4_t_cells.bed.gz"
  )
  inputs$source_bed_slugs <- c("b_cells", "cd4_t_cells")
  inputs$filtered_beds <- c(
    "gs://bucket/filtered/cd4_t_cells.filtered.bed.gz",
    "gs://bucket/filtered/b_cells.filtered.bed.gz"
  )
  inputs$filtered_bed_slugs <- c("cd4_t_cells", "b_cells")
  inputs$cell_types <- c("B cells", "CD4 T cells")
  inputs$cell_type_slugs <- c("b_cells", "cd4_t_cells")
  inputs$negative_summary <- "gs://bucket/reports/negative_summary.tsv.gz"
  inputs$gene_comparison <- "gs://bucket/reports/gene_comparison.tsv.gz"
  inputs$filter_metrics <- "gs://bucket/reports/filter_metrics.tsv"

  manifest <- do.call(build_cell_type_qtl_manifest, inputs)
  testthat::expect_identical(manifest$source_cpm_bed, inputs$source_beds)
  testthat::expect_identical(manifest$filtered_cpm_bed, rev(inputs$filtered_beds))
  testthat::expect_identical(
    manifest$reference_filter_metrics,
    rep(inputs$filter_metrics, 2L)
  )
})

testthat::test_that("QTL manifest requires equal nonzero array lengths", {
  inputs <- make_qtl_manifest_inputs()
  inputs$scaled_outliers <- inputs$scaled_outliers[-1L]
  testthat::expect_error(do.call(build_cell_type_qtl_manifest, inputs), "one equal nonzero length")
  testthat::expect_error(
    do.call(build_cell_type_qtl_manifest, purrr::map(inputs, ~ .x[0])),
    "one equal nonzero length"
  )
})

testthat::test_that("QTL manifest rejects missing and duplicate cell identities", {
  for (field in c("cell_types", "cell_type_slugs")) {
    for (bad_value in c("", NA_character_, "Monocytes", "monocytes")) {
      inputs <- make_qtl_manifest_inputs()
      inputs[[field]][[1L]] <- bad_value
      inputs[[field]][[2L]] <- bad_value
      testthat::expect_error(do.call(build_cell_type_qtl_manifest, inputs), "nonempty and unique")
    }
  }
})

testthat::test_that("QTL manifest rejects unsafe IDs and multiline labels", {
  inputs <- make_qtl_manifest_inputs()
  inputs$cell_type_slugs[[1L]] <- "monocytes/other"
  testthat::expect_error(do.call(build_cell_type_qtl_manifest, inputs), "safe")
  inputs <- make_qtl_manifest_inputs()
  inputs$cell_types[[1L]] <- "Monocytes\nother row"
  testthat::expect_error(do.call(build_cell_type_qtl_manifest, inputs), "control")
})

testthat::test_that("QTL manifest rejects incomplete paths and duplicate full paths", {
  for (bad_value in c("", NA_character_, "int.bed.gz", "gs://bucket", "gs://bucket/",
                      "gs://bad bucket/file", "gs://bucket/file\nother", "https://example.org/file")) {
    inputs <- make_qtl_manifest_inputs()
    inputs$int_beds[[1L]] <- bad_value
    testthat::expect_error(do.call(build_cell_type_qtl_manifest, inputs), "paths")
  }
  inputs <- make_qtl_manifest_inputs()
  inputs$int_beds[[2L]] <- inputs$int_beds[[1L]]
  testthat::expect_error(do.call(build_cell_type_qtl_manifest, inputs), "paths must be unique")
})

testthat::test_that("local smoke-test paths remain metadata without container localization", {
  inputs <- make_qtl_manifest_inputs()
  directory <- withr::local_tempdir()
  paths <- file.path(directory, c("monocytes.bed.gz", "cd4.bed.gz"))
  purrr::walk(paths, ~ writeLines("fixture", .x))
  inputs$int_beds <- paths
  manifest <- do.call(build_cell_type_qtl_manifest, inputs)
  testthat::expect_identical(manifest$int_bed, paths)
  inputs$int_beds[[1L]] <- file.path(directory, "absent.bed.gz")
  manifest <- do.call(build_cell_type_qtl_manifest, inputs)
  testthat::expect_identical(manifest$int_bed, inputs$int_beds)
})

testthat::test_that("QTL manifest CLI preserves cloud URLs from safely serialized arrays", {
  root <- normalizePath(testthat::test_path("..", "..", ".."))
  withr::local_dir(root)
  withr::local_envvar(CELL_TYPE_SPECIFIC_EXPRESSION_ROOT = file.path(
    root, "scripts", "cell_type_specific_expression"
  ))
  inputs <- make_qtl_manifest_inputs()
  inputs$int_beds[[1L]] <- "gs://test-bucket/path with spaces/it's-a-file.tsv"
  directory <- withr::local_tempdir()
  cli_names <- c(
    cell_types = "--cell-types", cell_type_slugs = "--cell-type-slugs",
    int_beds = "--int-beds", scaled_beds = "--scaled-beds",
    int_pcs = "--int-pcs", int_pcs_all = "--int-pcs-all",
    scaled_pcs = "--scaled-pcs", scaled_pcs_all = "--scaled-pcs-all",
    int_covariates = "--int-covariates", scaled_covariates = "--scaled-covariates",
    int_outliers = "--int-outliers", scaled_outliers = "--scaled-outliers"
  )
  arguments <- purrr::imap(inputs, function(values, argument) {
    path <- file.path(directory, paste0(argument, ".json"))
    jsonlite::write_json(values, path, auto_unbox = FALSE)
    c(cli_names[[argument]], path)
  }) |>
    unlist(use.names = FALSE)
  output <- file.path(directory, "manifest.tsv")
  status <- suppressWarnings(system2(
    "Rscript", shQuote(c("scripts/cell_type_specific_expression/build_qtl_manifest.R",
                        arguments, "--output", output))
  ))
  testthat::expect_identical(status, 0L)
  if (status != 0L) return(invisible())
  manifest <- readr::read_tsv(output, show_col_types = FALSE, name_repair = "minimal")
  testthat::expect_identical(names(manifest), qtl_manifest_columns)
  testthat::expect_identical(manifest[["entity:cell_type_id"]], c("monocytes", "cd4_t_cells"))
  testthat::expect_identical(unname(as.list(manifest[4:13])), unname(inputs[-(1:2)]))
  testthat::expect_true(all(is.na(manifest[14:18])))
})

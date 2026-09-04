source(testthat::test_path("..", "smoke", "reference_filter_metrics.R"), local = TRUE)
filter_script_root <- testthat::test_path("..", "..", "..", "scripts", "cell_type_specific_expression")
filter_test_environment <- environment()
invisible(lapply(file.path(filter_script_root, "R", c(
  "expression_bed.R", "bed_outputs.R", "reference_expression.R", "reference_filter.R",
  "reference_filter_plots.R", "reference_filter_io.R"
)), function(path) source(path, local = filter_test_environment)))

metric_row <- function(cell_type, status, reference_n_samples = 2L) {
  tibble::tibble(
    cell_type = cell_type, slug = "b_cells", comparison_status = status,
    deconvolution_n_samples = 3L, reference_n_samples = reference_n_samples,
    n_original = 5L, n_negative_excluded = 1L, n_reference_excluded = 1L,
    n_residual_excluded = 0L, n_retained = 3L, metric_set = "baseline",
    n_genes = 3L, pearson_r = 0.9, spearman_rho = 0.8,
    r_squared = 0.7, intercept = 0.1, slope = 1.1
  )
}

testthat::test_that("reference-enabled smoke rejects malformed and negative-only metrics", {
  malformed <- dplyr::select(metric_row("B cells", "available"), -n_retained)
  testthat::expect_error(
    validate_reference_filter_metrics(malformed, "B cells", "B cells", TRUE),
    "documented columns"
  )

  negative_only <- metric_row("B cells", "reference_not_provided", NA_integer_)
  negative_only$n_reference_excluded <- 0L
  negative_only$n_retained <- 4L
  testthat::expect_error(
    validate_reference_filter_metrics(negative_only, "B cells", "B cells", TRUE),
    "reference-enabled"
  )
})

testthat::test_that("smoke accepts generated reference metrics and verifies negative-only status", {
  tmp <- tempfile("smoke metrics ")
  dir.create(tmp)
  bed <- file.path(tmp, "b_cells.bed.gz")
  bed_table <- tibble::tibble(
    `#chr` = "1", start = 0:3, end = 1:4, gene_id = paste0("g", 1:4),
    s1 = c(1, 2, 4, 8), s2 = c(2, 4, 8, 16), s3 = c(3, 6, 12, 24)
  )
  readr::write_tsv(bed_table, bed)
  inventory <- tibble::tibble(
    logical_name = "cell_type_bed", path = basename(bed),
    sha256 = digest::digest(file = bed, algo = "sha256", serialize = FALSE),
    n_genes = 4L, n_samples = 3L, scale = "cpm",
    cell_group = "B cells", slug = "b_cells"
  )
  reference <- tibble::tibble(
    gene_id = paste0("g", 1:4), cell_type = "B cells", n_samples = 2L,
    mean_log2_cpm1 = c(1, 2, 3, 4), median_log2_cpm1 = c(1, 2, 3, 4)
  )
  generated <- filter_cell_type_beds(
    inventory, bed, file.path(tmp, "output"), reference_summary = reference
  )
  testthat::expect_silent(
    validate_reference_filter_metrics(generated$metrics, "B cells", "B cells", TRUE)
  )

  negative_only <- metric_row("B cells", "reference_not_provided", NA_integer_)
  negative_only$n_reference_excluded <- 0L
  negative_only$n_retained <- 4L
  negative_only[c("pearson_r", "spearman_rho", "r_squared", "intercept", "slope")] <- NA_real_
  negative_path <- tempfile(fileext = ".tsv")
  readr::write_tsv(negative_only, negative_path, na = "NA")
  negative_only <- read_reference_filter_metrics(negative_path)
  testthat::expect_silent(
    validate_reference_filter_metrics(negative_only, "B cells", "B cells", FALSE)
  )
})

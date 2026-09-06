testthat::test_that("the compact fixture retains three ordered groups with nonzero weights", {
  source(file.path(script_root, "bootstrap.R"), local = .GlobalEnv)
  proportions <- readr::read_tsv(
    file.path(pipeline_root, "tests/cell_type_specific_expression/fixtures/precomputed_proportions.tsv"),
    show_col_types = FALSE
  ) |>
    tibble::column_to_rownames("sample_id") |>
    as.matrix()
  combined <- combine_lm22_proportions(proportions)
  retained <- filter_and_adjust_groups(combined, mean_threshold = .12, zero_floor = 1e-6)
  testthat::expect_identical(colnames(retained$weights),
    c("B cells", "CD4 T cells", "Monocyte/myeloid"))
  testthat::expect_true(all(retained$weights > 0))
  testthat::expect_equal(unname(rowSums(retained$weights)), rep(1, 12))
})

testthat::test_that("manifest smoke validation accepts the configured compact threshold", {
  expressions <- parse(file.path(pipeline_root,
    "tests/cell_type_specific_expression/smoke/assert_deconvolution_outputs.R"))
  checks <- Filter(function(expression) {
    is.call(expression) && identical(expression[[1L]], as.name("require_true")) &&
      grepl('numeric_parameter("group_mean_threshold")', paste(deparse(expression), collapse = " "), fixed = TRUE)
  }, as.list(expressions))
  testthat::expect_length(checks, 1L)
  parameters <- list(min_lm22_overlap = .8, hspe_marker_fraction = .1,
    hspe_marker_method = "ratio", hspe_quantile_normalize = FALSE,
    group_mean_threshold = .12, zero_floor = 1e-6, tca_max_iters = 10,
    tca_parallel = FALSE, gene_type = c("protein_coding", "lncRNA"),
    random_seed = 20260901, scale = "cpm")
  environment <- list2env(list(parameters = parameters, mean_threshold = .12,
    expected_tca_parallel = FALSE, manifest = list(tca_version = "1.2.1"),
    numeric_parameter = function(name) as.numeric(parameters[[name]]),
    input_value = function(name) parameters[[name]],
    require_true = function(condition, message) if (!isTRUE(condition)) stop(message)))
  testthat::expect_no_error(eval(checks[[1L]], environment))
})

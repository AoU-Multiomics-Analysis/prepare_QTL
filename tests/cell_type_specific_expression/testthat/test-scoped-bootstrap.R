source(testthat::test_path("helper-load.R"), local = .GlobalEnv)

testthat::test_that("filter entrypoint loads its helpers without loading HSPE", {
  directory <- withr::local_tempdir()
  script <- file.path(directory, "filter_expression_genes.R")
  bootstrap <- file.path(script_root, "bootstrap.R")
  writeLines(c(
    paste0("source(", deparse(bootstrap), ")"),
    'stopifnot(exists("filter_expression_by_gene_types"))',
    'stopifnot(!exists("prepare_hspe_batches"))',
    'stopifnot(!exists("build_tca_fit_arguments"))'
  ), script)
  result <- system2(file.path(R.home("bin"), "Rscript"), shQuote(script),
                    stdout = TRUE, stderr = TRUE)
  testthat::expect_true(is.null(attr(result, "status")), info = paste(result, collapse = "\n"))
})

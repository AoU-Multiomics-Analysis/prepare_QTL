#!/usr/bin/env Rscript
# Run trusted stage tests against the selected image's bundled scripts, not the
# candidate checkout mounted by Docker. Fixtures and WDLs come from the checkout.
suppressPackageStartupMessages(library(tidyverse))

prepare_stage_workspace <- function(source_root, bundled, directory) {
  if (!file.exists(file.path(bundled, "bootstrap.R"))) stop("Missing bundled cell-type scripts")
  dir.create(directory)
  purrr::walk(c("tests", "workflows"), function(name) {
    if (!file.copy(file.path(source_root, name), directory, recursive = TRUE)) {
      stop("Cannot copy test inputs: ", name)
    }
  })
  dir.create(file.path(directory, "scripts"))
  if (!file.symlink(bundled, file.path(directory, "scripts", "cell_type_specific_expression"))) {
    stop("Cannot link bundled scripts")
  }
}

main <- function() {
  tests <- commandArgs(trailingOnly = TRUE)
  if (length(tests) == 0L || any(basename(tests) != tests)) {
    stop("Supply one or more test filenames")
  }
  # A missing scientific dependency must fail the image gate, not skip tests.
  purrr::walk(c("hspe", "TCA", "edgeR", "testthat"), function(package) {
    if (!requireNamespace(package, quietly = TRUE)) stop("Missing image dependency: ", package)
  })
  bundled <- "/opt/prepare_qtl/scripts/cell_type_specific_expression"
  source_root <- normalizePath(".")
  directory <- tempfile("cell-stage-tests-")
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  prepare_stage_workspace(source_root, bundled, directory)
  setwd(directory)
  on.exit(setwd(source_root), add = TRUE, after = FALSE)
  purrr::walk(tests, function(test) {
    path <- file.path("tests", "cell_type_specific_expression", "testthat", test)
    if (!file.exists(path)) stop("Missing registered test: ", test)
    testthat::test_file(path, reporter = "summary", stop_on_failure = TRUE)
  })
}

main()

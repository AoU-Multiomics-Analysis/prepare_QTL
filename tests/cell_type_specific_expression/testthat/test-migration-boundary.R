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

testthat::test_that("hspe CLI reads direct CPM BED and GTF", {
  script_path <- testthat::test_path("..", "..", "..", "scripts", "cell_type_specific_expression", "run_hspe.R")
  text <- paste(readLines(script_path, warn = FALSE), collapse = "\n")
  testthat::expect_match(text, '"--expression"', fixed = TRUE)
  testthat::expect_match(text, '"--gtf"', fixed = TRUE)
})

testthat::test_that("hspe CLI rejects incomplete direct CPM input mode", {
  script_path <- testthat::test_path("..", "..", "..", "scripts", "cell_type_specific_expression", "run_hspe.R")

  result <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      script_path,
      "--expression", "expression.bed",
      "--lm22", "lm22.tsv",
      "--output-dir", "outputs"
    ),
    stdout = TRUE,
    stderr = TRUE
  ))

  testthat::expect_true(!is.null(attr(result, "status")))
  testthat::expect_match(
    paste(result, collapse = "\n"),
    "Missing required options: gtf",
    fixed = TRUE
  )
})

testthat::test_that("hspe CLI validates LM22 before the large inputs", {
  script_path <- testthat::test_path("..", "..", "..", "scripts", "cell_type_specific_expression", "run_hspe.R")
  lm22_path <- tempfile(fileext = ".tsv")
  writeLines("Gene symbol\tB cells naive\nABCB4\t1", lm22_path)

  result <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      script_path,
      "--expression", "missing-expression.bed",
      "--gtf", "missing-annotation.gtf",
      "--lm22", lm22_path,
      "--output-dir", "outputs"
    ),
    stdout = TRUE,
    stderr = TRUE
  ))

  testthat::expect_true(!is.null(attr(result, "status")))
  testthat::expect_match(
    paste(result, collapse = "\n"),
    "LM22 must contain exactly the 22 standard LM22 columns",
    fixed = TRUE
  )
})

testthat::test_that("hspe CLI rejects the retired prepared-log option", {
  script_path <- testthat::test_path("..", "..", "..", "scripts", "cell_type_specific_expression", "run_hspe.R")

  result <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      script_path,
      "--bulk-log", "bulk.tsv"
    ),
    stdout = TRUE,
    stderr = TRUE
  ))

  testthat::expect_true(!is.null(attr(result, "status")))
  testthat::expect_match(
    paste(result, collapse = "\n"),
    'long flag "bulk-log" is invalid',
    fixed = TRUE
  )
})

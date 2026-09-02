args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0L) {
  testthat::test_dir(
    "tests/cell_type_specific_expression/testthat",
    reporter = "summary",
    stop_on_failure = TRUE
  )
} else {
  purrr::walk(
    args,
    ~ testthat::test_file(.x, reporter = "summary", stop_on_failure = TRUE)
  )
}

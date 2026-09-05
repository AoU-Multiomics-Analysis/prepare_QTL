# Exercise the smoke reader without running the whole workflow assertion script.
smoke_expressions <- parse(testthat::test_path("..", "smoke", "assert_qtl_outputs.R"))
reader_definition <- Filter(function(expression) {
  is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
    identical(expression[[2L]], as.name("read_qtl_bed"))
}, as.list(smoke_expressions))
stopifnot(length(reader_definition) == 1L)
eval(reader_definition[[1L]])
require_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

testthat::test_that("smoke reader preserves the numeric inputs used by PrepareExpression", {
  expected_samples <- paste0("s", seq_len(12L))
  environment(read_qtl_bed) <- environment()
  values <- c(118.013452573863, 118.013452575535, 118.013452569707,
    118.013452572676, 118.013452571262, 118.013452571942,
    118.013452574756, 118.013452577389, 118.013452576283,
    118.013452573087, 118.013452573520, 118.013452576932)
  bed <- tibble::tibble(`#chr` = "chr1", start = 0L, end = 1L, gene_id = "g1") |>
    dplyr::bind_cols(tibble::as_tibble(as.list(stats::setNames(values, expected_samples))))
  path <- tempfile(fileext = ".bed")
  on.exit(unlink(path), add = TRUE)
  data.table::fwrite(bed, path, sep = "\t")
  # PrepareExpression reads CPM with fread. Near-constant genes amplify even
  # one-bit differences between numeric parsers during standardization.
  production_input <- data.table::fread(path, check.names = FALSE) |>
    dplyr::select(dplyr::all_of(expected_samples)) |>
    as.matrix()
  smoke_input <- read_qtl_bed(path, "CPM") |>
    dplyr::select(dplyr::all_of(expected_samples)) |>
    as.matrix()
  testthat::expect_equal(unname(smoke_input), unname(production_input), tolerance = 0)
  testthat::expect_equal(
    as.numeric(scale(log2(smoke_input[1L, ] + 1))),
    as.numeric(scale(log2(production_input[1L, ] + 1))), tolerance = 1e-7
  )
})

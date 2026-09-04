testthat::test_that("reference smoke fixture supplies raw counts and all mapped populations", {
  directory <- tempfile("reference fixture ")
  dir.create(directory)
  output <- file.path(directory, "counts.tsv.gz")
  generator <- file.path(pipeline_root, "tests/cell_type_specific_expression/generate_reference_fixture.R")
  bed <- file.path(pipeline_root, "tests/cell_type_specific_expression/fixtures/synthetic_expression.bed")
  logs <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(generator, bed, output)), stdout = TRUE, stderr = TRUE))
  testthat::expect_null(attr(logs, "status"), info = paste(logs, collapse = "\n"))
  testthat::expect_true(file.exists(output))
  if (!file.exists(output)) return(invisible(NULL))
  counts <- readr::read_tsv(output, show_col_types = FALSE)
  original <- readr::read_tsv(bed, show_col_types = FALSE)
  testthat::expect_identical(counts$gene_id, original$gene_id)
  testthat::expect_setequal(sub("\\.[0-9]+$", "", names(counts)[-1]),
    c("NveB", "MemB", "CD4T", "CD8T", "NK", "Mono", "MonoNonClassical",
      "Neut", "Eo", "myDC", "myDC123", "pDC"))
  values <- as.matrix(counts[-1])
  testthat::expect_true(all(is.finite(values) & values > 0 & values == floor(values)))
  testthat::expect_true(all(apply(values, 2, stats::sd) > 0))
})

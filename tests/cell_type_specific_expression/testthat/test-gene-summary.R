source(testthat::test_path("helper-load.R"), local = .GlobalEnv)
source(file.path(script_root, "bootstrap.R"), local = .GlobalEnv)

make_summary_fixture <- function(single_sample = FALSE) {
  directory <- tempfile("summary inputs ")
  dir.create(directory)
  coordinates <- tibble::tibble(
    `#chr` = c("chr2", "chr1", "chrX"), start = c(20L, 10L, 30L),
    end = c(21L, 11L, 31L), gene_id = c("g2.1", "g1", "g3")
  )
  values <- matrix(c(0, 2, 4, 6, 3, 3, 3, 3, -2, 0, 2, 8),
    nrow = 3, byrow = TRUE,
    dimnames = list(coordinates$gene_id, c("S4", "S1", "S3", "S2")))
  if (single_sample) values <- values[, 1, drop = FALSE]
  result <- write_cell_type_beds(
    list("CD8 T cells" = values, "B cells" = values + 10), coordinates, directory
  )
  result$directory <- directory
  result
}

testthat::test_that("summaries retain gene/cell identity and compute across samples in CPM", {
  fixture <- make_summary_fixture()
  output <- file.path(fixture$directory, "summary.tsv.gz")
  summarize_cell_type_beds(fixture$inventory, rev(unname(fixture$paths)), output, chunk_size = 1L)
  result <- readr::read_tsv(output, show_col_types = FALSE)
  testthat::expect_identical(result$cell_type, rep(c("CD8 T cells", "B cells"), each = 3))
  testthat::expect_identical(result$gene_id, rep(c("g2.1", "g1", "g3"), 2))
  testthat::expect_identical(result[["#chr"]], rep(c("chr2", "chr1", "chrX"), 2))
  testthat::expect_equal(result$start, rep(c(20, 10, 30), 2))
  testthat::expect_equal(result$end, rep(c(21, 11, 31), 2))
  testthat::expect_true(all(result$n_samples == 4 & result$scale == "cpm"))
  testthat::expect_equal(result$mean_cpm, c(3, 3, 2, 13, 13, 12))
  testthat::expect_equal(result$median_cpm, c(3, 3, 1, 13, 13, 11))
  testthat::expect_equal(result$sd_cpm, rep(c(sqrt(20 / 3), 0, sqrt(56 / 3)), 2))
  testthat::expect_equal(result$se_mean_cpm, result$sd_cpm / 2)
  testthat::expect_equal(result$q1_cpm, c(1.5, 3, -0.5, 11.5, 13, 9.5))
  testthat::expect_equal(result$q3_cpm, c(4.5, 3, 3.5, 14.5, 13, 13.5))
  testthat::expect_equal(result$iqr_cpm, c(3, 0, 4, 3, 0, 4))
  second_output <- file.path(fixture$directory, "larger_chunks.tsv.gz")
  summarize_cell_type_beds(fixture$inventory, unname(fixture$paths), second_output, chunk_size = 10L)
  testthat::expect_identical(readr::read_file(output), readr::read_file(second_output))
})

testthat::test_that("one sample has undefined SD and SE, not zero uncertainty", {
  fixture <- make_summary_fixture(single_sample = TRUE)
  output <- file.path(fixture$directory, "summary.tsv.gz")
  summarize_cell_type_beds(fixture$inventory, unname(fixture$paths), output)
  result <- readr::read_tsv(output, show_col_types = FALSE)
  testthat::expect_true(all(is.na(result$sd_cpm) & is.na(result$se_mean_cpm)))
  testthat::expect_true(all(result$iqr_cpm == 0 & result$n_samples == 1))
  testthat::expect_equal(result$mean_cpm, c(0, 3, -2, 10, 13, 8))
})

testthat::test_that("summary rejects invalid BED data and inventory mismatches", {
  fixture <- make_summary_fixture()
  output <- file.path(fixture$directory, "summary.tsv.gz")
  inventory <- fixture$inventory
  inventory$scale <- "tpm"
  testthat::expect_error(summarize_cell_type_beds(inventory, unname(fixture$paths), output), "scale")
  inventory <- fixture$inventory
  inventory$n_samples <- 5
  testthat::expect_error(summarize_cell_type_beds(inventory, unname(fixture$paths), output), "sample")
  inventory <- fixture$inventory
  inventory$n_genes <- 4
  testthat::expect_error(summarize_cell_type_beds(inventory, unname(fixture$paths), output), "gene")
  testthat::expect_error(summarize_cell_type_beds(fixture$inventory, unname(fixture$paths)[1], output), "count")
  original <- readr::read_tsv(fixture$paths[[1]], show_col_types = FALSE)
  for (bad_value in c(NA_real_, Inf, NaN)) {
    bad <- original
    bad$S4[1] <- bad_value
    readr::write_tsv(bad, fixture$paths[[1]])
    testthat::expect_error(summarize_cell_type_beds(fixture$inventory, unname(fixture$paths), output), "finite|missing|parse")
  }
  bad <- original
  bad$gene_id[3] <- bad$gene_id[1]
  readr::write_tsv(bad, fixture$paths[[1]])
  testthat::expect_error(summarize_cell_type_beds(fixture$inventory, unname(fixture$paths), output, chunk_size = 1L), "unique|duplicate")
  bad <- original
  names(bad)[6] <- names(bad)[5]
  readr::write_tsv(bad, fixture$paths[[1]])
  testthat::expect_error(summarize_cell_type_beds(fixture$inventory, unname(fixture$paths), output), "unique|duplicate")
})

testthat::test_that("summary CLI writes a compressed output and logs completion", {
  fixture <- make_summary_fixture()
  inventory_path <- file.path(fixture$directory, "inventory.tsv")
  paths_path <- file.path(fixture$directory, "beds.txt")
  output <- file.path(fixture$directory, "summary.tsv.gz")
  readr::write_tsv(fixture$inventory, inventory_path)
  writeLines(rev(unname(fixture$paths)), paths_path)
  cli <- file.path(script_root, "summarize_cell_type_beds.R")
  logs <- system2(file.path(R.home("bin"), "Rscript"),
    shQuote(c(cli, inventory_path, paths_path, output)), stdout = TRUE, stderr = TRUE)
  testthat::expect_null(attr(logs, "status"))
  testthat::expect_true(file.exists(output))
  testthat::expect_true(any(grepl("stage=summarize_cell_type_beds.*completion_time", logs)))
  if (file.exists(output)) {
    testthat::expect_equal(nrow(readr::read_tsv(output, show_col_types = FALSE)), 6)
    connection <- file(output, "rb")
    testthat::expect_identical(readBin(connection, "raw", 2), as.raw(c(31, 139)))
    close(connection)
  }
})

source(testthat::test_path("helper-load.R"), local = .GlobalEnv)
source(file.path(script_root, "bootstrap.R"), local = .GlobalEnv)

batch_test_inputs <- function() {
  reference <- rbind(
    type_C = c(200, 100, 10, 10, 20, 20),
    type_A = c(20, 20, 200, 100, 10, 10),
    type_B = c(10, 10, 20, 20, 200, 100)
  )
  colnames(reference) <- paste0("G", 1:6)
  weights <- rbind(S3 = c(.6, .3, .1), S1 = c(.15, .25, .6), S2 = c(.2, .5, .3))
  list(Y = log2(weights %*% reference), references = log2(reference),
       overlap_report = tibble::tibble(gene_symbol = colnames(reference)),
       transformed_lm22 = t(log2(reference)), overlap_count = 6L,
       overlap_fraction = 1, quantile_normalize = FALSE, lm22_qc = list())
}

testthat::test_that("preparation selects markers before splitting samples", {
  testthat::skip_if_not_installed("hspe")
  prepared <- prepare_hspe_batches(batch_test_inputs(), batch_size = 2, marker_fraction = .5)
  testthat::expect_identical(colnames(prepared$batches[[1]]$Y), c("G1", "G3", "G5"))
  testthat::expect_identical(rownames(prepared$batches[[1]]$Y), c("S3", "S1"))
  testthat::expect_identical(rownames(prepared$batches[[2]]$Y), "S2")
  testthat::expect_equal(unname(prepared$shared$markers), list(1L, 2L, 3L))
  testthat::expect_equal(dim(prepared$shared$references), c(3L, 3L))
  testthat::expect_null(prepared$shared$Y)
  testthat::expect_identical(prepared$shared$samples$sample_id, c("S3", "S1", "S2"))
  for (size in list(0, -1, 1.5, NA, Inf)) {
    testthat::expect_error(prepare_hspe_batches(batch_test_inputs(), size), "batch_size")
  }
})

testthat::test_that("seeds are stable across batching and sample ordering", {
  ids <- c("S3", "S1", "S2")
  seeds <- hspe_sample_seeds(ids, 123L)
  testthat::expect_identical(seeds, rev(hspe_sample_seeds(rev(ids), 123L)))
  testthat::expect_identical(seeds[2], hspe_sample_seeds("S1", 123L))
  testthat::expect_true(all(is.finite(seeds) & seeds > 0 & seeds <= .Machine$integer.max))
  testthat::expect_false(identical(seeds, hspe_sample_seeds(ids, 124L)))
  testthat::expect_error(hspe_sample_seeds(ids, 0), "random_seed")
})

testthat::test_that("batched fits match a single batch and restore original sample order", {
  testthat::skip_if_not_installed("hspe")
  inputs <- batch_test_inputs()
  split <- prepare_hspe_batches(inputs, batch_size = 2, marker_fraction = .5, random_seed = 123)
  together <- prepare_hspe_batches(inputs, batch_size = 10, marker_fraction = .5, random_seed = 123)
  results <- purrr::map(split$batches, ~ fit_hspe_batch(.x, split$shared))
  merged <- merge_hspe_batches(rev(results), split$shared)
  single <- fit_hspe_batch(together$batches[[1]], together$shared)
  testthat::expect_identical(merged$proportions, single$proportions)
  upstream <- hspe::hspe(Y = inputs$Y[1, , drop = FALSE], references = inputs$references,
    n_markers = .5, seed = split$shared$samples$random_seed[1], sto = TRUE)
  testthat::expect_identical(as.numeric(merged$proportions[1, ]), as.numeric(upstream$estimates))
  testthat::expect_equal(unname(merged$proportions),
    rbind(c(.6, .3, .1), c(.15, .25, .6), c(.2, .5, .3)), tolerance = 1e-4)
  testthat::expect_identical(merged$diagnostics$sample_id, c("S3", "S1", "S2"))
  testthat::expect_identical(merged$diagnostics$random_seed, split$shared$samples$random_seed)
  testthat::expect_equal(merged$metadata$batch_count, 2)
  testthat::expect_equal(merged$metadata$sample_count, 3)
  testthat::expect_error(merge_hspe_batches(results[1], split$shared), "missing|sample")
  testthat::expect_error(merge_hspe_batches(c(results, results[1]), split$shared), "duplicate")
  wrong <- results
  colnames(wrong[[1]]$proportions) <- rev(colnames(wrong[[1]]$proportions))
  testthat::expect_error(merge_hspe_batches(wrong, split$shared), "cell.type|column")
  wrong <- results
  wrong[[1]]$proportions[1, 1] <- NA_real_
  testthat::expect_error(merge_hspe_batches(wrong, split$shared), "finite")
  wrong <- results
  wrong[[1]]$diagnostics$random_seed[1] <- 1L
  testthat::expect_error(merge_hspe_batches(wrong, split$shared), "seed")
  wrong <- results
  wrong[[1]]$proportions[1, ] <- c(.8, .8, .8)
  testthat::expect_error(merge_hspe_batches(wrong, split$shared), "sum to one")
  wrong <- results
  wrong[[1]]$optimizer_version <- "different"
  testthat::expect_error(merge_hspe_batches(wrong, split$shared), "versions")
})

testthat::test_that("preparation preserves cohort-normalized marker values", {
  testthat::skip_if_not_installed("hspe")
  inputs <- batch_test_inputs()
  joint <- limma::normalizeBetweenArrays(cbind(t(inputs$references), t(inputs$Y)))
  inputs$references <- t(joint[, 1:3])
  inputs$Y <- t(joint[, 4:6])
  inputs$quantile_normalize <- TRUE
  prepared <- prepare_hspe_batches(inputs, batch_size = 2, marker_fraction = .5)
  assembled <- do.call(rbind, purrr::map(prepared$batches, "Y"))
  testthat::expect_identical(assembled, inputs$Y[, colnames(assembled), drop = FALSE])
  testthat::expect_identical(prepared$shared$references,
    inputs$references[, colnames(assembled), drop = FALSE])
  testthat::expect_true(prepared$shared$metadata$quantile_normalize)
})

testthat::test_that("batch commands write localized inputs and merged outputs", {
  testthat::skip_if_not_installed("hspe")
  directory <- tempfile("hspe batch cli ")
  dir.create(directory)
  cli <- function(script, args) {
    output <- system2(file.path(R.home("bin"), "Rscript"),
      shQuote(c(file.path(script_root, script), args)), stdout = TRUE, stderr = TRUE)
    testthat::expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  }
  fixture <- file.path(pipeline_root, "tests/cell_type_specific_expression/fixtures")
  cli("prepare_hspe_batches.R", c(
    "--expression", file.path(fixture, "synthetic_expression.bed"),
    "--gtf", file.path(fixture, "synthetic.gtf"),
    "--lm22", file.path(fixture, "synthetic_signature.tsv"),
    "--batch-size", "5", "--output-dir", directory))
  testthat::expect_true(file.exists(file.path(directory, "hspe_prepared.rds")))
  paths <- list.files(directory, "^batch_.*[.]rds$", full.names = TRUE)
  testthat::expect_length(paths, 3)
  if (length(paths) == 3L) {
    shared <- readRDS(file.path(directory, "hspe_prepared.rds"))
    batches <- purrr::map(paths, readRDS)
    testthat::expect_identical(purrr::map_int(batches, ~ nrow(.x$Y)), c(5L, 5L, 2L))
    testthat::expect_identical(colnames(batches[[1]]$Y), colnames(shared$references))
    testthat::expect_true(ncol(batches[[1]]$Y) < shared$metadata$overlap_count)
  }
  small <- prepare_hspe_batches(batch_test_inputs(), batch_size = 2, marker_fraction = .5)
  shared_path <- file.path(directory, "small shared.rds")
  saveRDS(small$shared, shared_path)
  result_paths <- purrr::imap_chr(small$batches, function(batch, index) {
    batch_path <- file.path(directory, paste0("small batch ", index, ".rds"))
    out <- file.path(directory, paste0("result ", index))
    saveRDS(batch, batch_path)
    cli("run_hspe_batch.R", c(shared_path, batch_path, out))
    file.path(out, "hspe_batch_result.rds")
  })
  result_list <- file.path(directory, "result paths.txt")
  writeLines(rev(result_paths), result_list)
  merged_dir <- file.path(directory, "merged outputs")
  cli("merge_hspe_batches.R", c(shared_path, result_list, merged_dir))
  testthat::expect_true(file.exists(file.path(merged_dir, "hspe_metadata.json")))
  p <- read_numeric_matrix(file.path(merged_dir, "hspe_proportions.tsv"), "sample_id")
  testthat::expect_identical(rownames(p), c("S3", "S1", "S2"))
  testthat::expect_equal(unname(p),
    rbind(c(.6, .3, .1), c(.15, .25, .6), c(.2, .5, .3)), tolerance = 1e-4)
  diagnostics <- readr::read_tsv(file.path(merged_dir, "hspe_sample_diagnostics.tsv"),
                                show_col_types = FALSE)
  testthat::expect_identical(diagnostics$sample_id, rownames(p))
})

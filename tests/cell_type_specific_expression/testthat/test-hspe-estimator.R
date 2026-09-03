source(testthat::test_path("helper-load.R"), local = .GlobalEnv)
source(file.path(script_root, "bootstrap.R"), local = .GlobalEnv)

make_hspe_test_inputs <- function() {
  reference <- matrix(
    c(200, 10, 20, 20, 200, 10, 10, 20, 200), nrow = 3, byrow = TRUE,
    dimnames = list(c("type_C", "type_A", "type_B"), c("G1", "G2", "G3"))
  )
  weights <- matrix(c(0.6, 0.3, 0.1, 0.15, 0.25, 0.6), nrow = 2, byrow = TRUE)
  bulk <- weights %*% reference
  rownames(bulk) <- c("S2", "S1")
  list(
    Y = log2(bulk), references = log2(reference),
    overlap_report = tibble::tibble(gene_symbol = colnames(reference))
  )
}

testthat::test_that("HSPE recovers known mixtures without changing sample or cell order", {
  testthat::skip_if_not_installed("hspe")
  inputs <- make_hspe_test_inputs()
  fit <- estimate_hspe(inputs, marker_fraction = 1, random_seed = 123L)

  testthat::expect_equal(
    unname(fit$proportions),
    matrix(c(0.6, 0.3, 0.1, 0.15, 0.25, 0.6), nrow = 2, byrow = TRUE),
    tolerance = 1e-4
  )
  testthat::expect_identical(rownames(fit$proportions), c("S2", "S1"))
  testthat::expect_identical(colnames(fit$proportions), c("type_C", "type_A", "type_B"))
  testthat::expect_equal(unname(rowSums(fit$proportions)), c(1, 1), tolerance = 1e-8)
  testthat::expect_identical(fit$metadata$hspe_version, "0.1")
  testthat::expect_identical(fit$metadata$optimizer, "DEoptimR")
  testthat::expect_identical(fit$metadata$random_seed, 123L)
  testthat::expect_identical(fit$markers$gene_symbol, c("G1", "G2", "G3"))

  repeated <- estimate_hspe(inputs, marker_fraction = 1, random_seed = 123L)
  testthat::expect_identical(repeated$proportions, fit$proportions)
})

testthat::test_that("HSPE preserves one-sample matrices and rejects invalid seeds", {
  testthat::skip_if_not_installed("hspe")
  inputs <- make_hspe_test_inputs()
  inputs$Y <- inputs$Y[1, , drop = FALSE]
  fit <- estimate_hspe(inputs, marker_fraction = 1, random_seed = 123L)
  testthat::expect_identical(dim(fit$proportions), c(1L, 3L))
  testthat::expect_identical(rownames(fit$proportions), "S2")
  for (seed in list(NA, Inf, 0, -1, 1.5, numeric(), c(1, 2))) {
    testthat::expect_error(estimate_hspe(inputs, random_seed = seed), "random_seed")
  }
})

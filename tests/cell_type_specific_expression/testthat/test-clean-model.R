source(testthat::test_path("helper-load.R"), local = .GlobalEnv)
source(file.path(script_root, "bootstrap.R"), local = .GlobalEnv)

cleanup_fixture <- function() {
  genes <- c("keep1", "bad", "keep2")
  samples <- c("s1", "s2", "s3")
  groups <- c("A", "B")
  zero_columns <- matrix(numeric(), 3L, 0L, dimnames = list(genes, NULL))
  list(
    W = matrix(c(.25, .75, .5, .5, .75, .25), 3L, byrow = TRUE,
               dimnames = list(samples, groups)),
    mus_hat = matrix(c(1, 2, 3, 4, 5, 6), 3L, byrow = TRUE,
                     dimnames = list(genes, groups)),
    sigmas_hat = matrix(c(1, 2, .01, 699957.8, 2, 3), 3L, byrow = TRUE,
                        dimnames = list(genes, groups)),
    tau_hat = 1,
    deltas_hat = zero_columns, gammas_hat = zero_columns,
    deltas_hat_pvals = zero_columns, gammas_hat_pvals = zero_columns,
    gammas_hat_pvals.joint = zero_columns,
    C1 = matrix(numeric(), 3L, 0L, dimnames = list(samples, NULL)),
    C2 = matrix(numeric(), 3L, 0L, dimnames = list(samples, NULL)),
    expression_scale = "cpm", tca_parallel = FALSE
  )
}

testthat::test_that("restart accepts fitted and cleaned CPM models without changing weights", {
  original <- cleanup_fixture()
  cleaned <- prepare_restart_model(original)
  testthat::expect_identical(cleaned$model$W, original$W)
  testthat::expect_identical(cleaned$model$C2, original$C2)
  again <- prepare_restart_model(cleaned$model)
  testthat::expect_identical(again$model, cleaned$model)
  testthat::expect_identical(again$report, cleaned$report)
  # TCA 1.2.1 creates unnamed zero-column matrices when covariates are absent.
  no_covariates <- cleanup_fixture()
  rownames(no_covariates$C1) <- NULL
  rownames(no_covariates$C2) <- NULL
  testthat::expect_identical(prepare_restart_model(no_covariates)$model$C2, no_covariates$C2)
  original$expression_scale <- "log2_cpm"
  testthat::expect_error(prepare_restart_model(original), "linear CPM")
  original$expression_scale <- NULL
  testthat::expect_error(prepare_restart_model(original), "linear CPM")
  invalid <- cleanup_fixture()
  invalid$W[1, 1] <- NA_real_
  testthat::expect_error(prepare_restart_model(invalid), "weights")
})

testthat::test_that("restart aligns BED genes but does not silently drop missing model genes", {
  model <- prepare_restart_model(cleanup_fixture())$model
  X <- matrix(1:9, 3L, dimnames = list(c("keep2", "extra", "keep1"), c("s1", "s2", "s3")))
  testthat::expect_identical(align_restart_expression(X, model), X[c(3, 1), , drop = FALSE])
  testthat::expect_error(align_restart_expression(X[-1, ], model), "missing.*keep2")
  testthat::expect_error(align_restart_expression(X[, 3:1], model), "sample order")
})

testthat::test_that("restart cleanup CLI accepts a cleaned model and writes its stored weights", {
  withr::local_dir(pipeline_root)
  work <- withr::local_tempdir()
  model <- clean_tca_model(cleanup_fixture())$model
  path <- file.path(work, "model.rds")
  saveRDS(model, path)
  output <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(
    "scripts/cell_type_specific_expression/clean_tca_model.R", "--reuse-model", "--model", path,
    "--output-dir", file.path(work, "output")
  )), stdout = TRUE, stderr = TRUE)
  testthat::expect_null(attr(output, "status"), info = paste(output, collapse = "\n"))
  testthat::expect_identical(readRDS(file.path(work, "output", "tca_model_cleaned.rds")), model)
  testthat::expect_equal(read_numeric_matrix(file.path(work, "output", "tca_weights.tsv"), "sample_id"), model$W)
  qc <- build_restart_qc_summary(tibble::tibble(metric = "gene_count", value = 2), model$W, model)
  testthat::expect_true(is.na(qc$value[qc$metric == "tca_convergence"]))
  testthat::expect_identical(qc$status[qc$metric == "tca_convergence"], "unavailable_model_restart")
})

testthat::test_that("cleanup removes the failing gene from every gene parameter only", {
  original <- cleanup_fixture()
  result <- clean_tca_model(original)
  testthat::expect_identical(result$report$gene_id, "bad")
  testthat::expect_identical(result$report$reason, "computationally_singular_variance")
  testthat::expect_equal(result$report$reciprocal_condition, 2.041062e-16,
                         tolerance = 1e-22, scale = 1)
  testthat::expect_equal(result$report$threshold, .Machine$double.eps)
  gene_fields <- c("mus_hat", "sigmas_hat", "deltas_hat", "gammas_hat",
                   "deltas_hat_pvals", "gammas_hat_pvals", "gammas_hat_pvals.joint")
  for (field in gene_fields) {
    testthat::expect_identical(result$model[[field]], original[[field]][c(1, 3), , drop = FALSE])
  }
  for (field in setdiff(names(original), gene_fields)) {
    testthat::expect_identical(result$model[[field]], original[[field]])
  }
  testthat::expect_identical(original, cleanup_fixture())
  testthat::expect_identical(result$model$gene_filter$original_gene_ids,
                             c("keep1", "bad", "keep2"))
})

testthat::test_that("cleanup uses relative variance and retains a boundary or one-gene model", {
  model <- cleanup_fixture()
  model$sigmas_hat[1, ] <- c(1e-100, 2e-100)
  model$sigmas_hat[2, ] <- c(sqrt(.Machine$double.eps), 1)
  result <- clean_tca_model(model)
  testthat::expect_equal(nrow(result$report), 0L)
  testthat::expect_identical(result$model$sigmas_hat, model$sigmas_hat)
  model$sigmas_hat[1:2, ] <- 0
  result <- clean_tca_model(model)
  testthat::expect_identical(rownames(result$model$mus_hat), "keep2")
  testthat::expect_identical(dim(result$model$deltas_hat), c(1L, 0L))
})

testthat::test_that("cleanup retains aligned nonempty covariate effects and p-values", {
  model <- cleanup_fixture()
  slopes <- matrix(seq_len(12L), 3L, dimnames = list(rownames(model$mus_hat), paste0("effect", 1:4)))
  model$gammas_hat <- slopes
  model$gammas_hat_pvals <- slopes / 20
  model$deltas_hat <- slopes[, 1, drop = FALSE]
  model$deltas_hat_pvals <- slopes[, 1, drop = FALSE] / 20
  model$gammas_hat_pvals.joint <- slopes[, 1:2, drop = FALSE] / 20
  result <- clean_tca_model(model)
  for (field in c("gammas_hat", "gammas_hat_pvals", "deltas_hat", "deltas_hat_pvals",
                  "gammas_hat_pvals.joint")) {
    testthat::expect_identical(result$model[[field]], model[[field]][c(1, 3), , drop = FALSE])
  }
  testthat::expect_error(clean_tca_model(result$model), "already contains a gene filter")
})

testthat::test_that("invalid parameters and wholly excluded models fail explicitly", {
  for (value in c(0, -1, NA_real_, Inf)) {
    model <- cleanup_fixture()
    model$tau_hat <- value
    testthat::expect_error(clean_tca_model(model), "tau_hat.*finite and positive")
  }
  for (value in c(NA_real_, Inf, -1, 1e200)) {
    model <- cleanup_fixture()
    model$sigmas_hat[1, 1] <- value
    testthat::expect_error(clean_tca_model(model), "finite|nonnegative")
  }
  model <- cleanup_fixture()
  model$sigmas_hat[,] <- 0
  testthat::expect_error(clean_tca_model(model), "No genes remain")
  model <- cleanup_fixture()
  model$mus_hat <- model$mus_hat[3:1, ]
  testthat::expect_error(clean_tca_model(model), "gene order")
  model <- cleanup_fixture()
  rownames(model$sigmas_hat)[2] <- "keep1"
  testthat::expect_error(clean_tca_model(model), "unique")
})

testthat::test_that("expression alignment removes only recorded exclusions and checks identities", {
  clean <- clean_tca_model(cleanup_fixture())$model
  X <- matrix(1:9, 3L, dimnames = list(c("keep1", "bad", "keep2"), c("s1", "s2", "s3")))
  testthat::expect_identical(align_expression_to_tca_model(X, clean), X[c(1, 3), , drop = FALSE])
  testthat::expect_error(align_expression_to_tca_model(X[3:1, ], clean), "gene order")
  wrong <- X
  rownames(wrong)[2] <- "unrecorded"
  testthat::expect_error(align_expression_to_tca_model(wrong, clean), "gene order")
  testthat::expect_error(align_expression_to_tca_model(X[, 3:1], clean), "sample order")
  clean$gene_filter$excluded_gene_ids <- character()
  testthat::expect_error(align_expression_to_tca_model(X, clean), "exclusion")
})

testthat::test_that("cleanup CLI saves a final model consumed by the real export CLI", {
  testthat::skip_if_not_installed("TCA", minimum_version = "1.2.1")
  testthat::skip_if_not_installed("ggplot2")
  withr::local_dir(pipeline_root)
  work <- tempfile("clean model ")
  dir.create(work)
  model <- cleanup_fixture()
  raw_path <- file.path(work, "raw.rds")
  saveRDS(model, raw_path)
  run_script <- function(script, arguments) {
    output <- system2(file.path(R.home("bin"), "Rscript"),
                      shQuote(c(file.path("scripts", "cell_type_specific_expression", script), arguments)),
                      stdout = TRUE, stderr = TRUE)
    status <- attr(output, "status")
    testthat::expect_true(is.null(status) || status == 0L, info = paste(output, collapse = "\n"))
    output
  }
  invisible(run_script("clean_tca_model.R", c("--model", raw_path, "--output-dir", file.path(work, "clean"))))
  clean_path <- file.path(work, "clean", "tca_model_cleaned.rds")
  testthat::expect_true(file.exists(clean_path))
  if (!file.exists(clean_path)) return(invisible(NULL))
  clean <- readRDS(clean_path)
  testthat::expect_identical(readRDS(raw_path), model)
  report <- readr::read_tsv(file.path(work, "clean", "tca_numerical_excluded_genes.tsv"),
                            show_col_types = FALSE)
  testthat::expect_identical(report$gene_id, "bad")
  coordinates <- tibble::tibble(`#chr` = "chr1", start = c(0L, 10L, 20L),
                                 end = c(1L, 11L, 21L), gene_id = rownames(model$mus_hat))
  X <- matrix(c(3, 4, 5, 10, 20, 30, 6, 7, 8), 3L, byrow = TRUE,
              dimnames = list(coordinates$gene_id, rownames(model$W)))
  bed_path <- file.path(work, "expression.bed")
  write_expression_bed(bed_path, coordinates, X)
  weights_path <- file.path(work, "weights.tsv")
  write_numeric_matrix(model$W, weights_path, "sample_id")
  testthat::expect_error(TCA::tensor(X, model, verbose = FALSE, log_file = NULL),
                         "computationally singular")
  invisible(run_script("export_tca_beds.R", c("--expression", bed_path, "--model", clean_path,
      "--weights", weights_path, "--output-dir", file.path(work, "export"))))
  exported <- read_expression_bed(file.path(work, "export", "a.bed.gz"))
  testthat::expect_identical(exported$coordinates$gene_id, c("keep1", "keep2"))
  testthat::expect_identical(exported$coordinates$start, c(0L, 20L))
  testthat::expect_equal(as.numeric(exported$cpm[1, 1]), 1 + .25 * 1.25 / 3.3125)

  # Restart uses retained model genes even if the source BED has a different row order.
  restart_bed <- file.path(work, "restart.bed")
  write_expression_bed(restart_bed, coordinates[3:1, ], X[3:1, , drop = FALSE])
  invisible(run_script("export_tca_beds.R", c("--expression", restart_bed, "--model", clean_path,
      "--weights", weights_path, "--reuse-model", "--output-dir", file.path(work, "restart"))))
  restarted <- read_expression_bed(file.path(work, "restart", "a.bed.gz"))
  testthat::expect_identical(restarted$coordinates, exported$coordinates)
  testthat::expect_equal(restarted$cpm, exported$cpm)
})

source(testthat::test_path("helper-load.R"), local = .GlobalEnv)
source(file.path(script_root, "bootstrap.R"), local = .GlobalEnv)

testthat::test_that("fit CLI writes a usable CPM model, excluded-gene report and log", {
  testthat::skip_if_not_installed("TCA")
  set.seed(20260901)
  data <- TCA::test_data(24, 20, 3, 0, 0, 0.01)
  X <- abs(data$X) + 1
  W <- data$W
  rownames(X) <- paste0("g", seq_len(nrow(X)))
  colnames(X) <- rownames(W) <- paste0("s", seq_len(ncol(X)))
  colnames(W) <- paste0("type", seq_len(ncol(W)))
  directory <- tempfile("fit cli ")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  bed <- tibble::tibble(`#chr` = "chr1", start = seq_len(nrow(X)),
                        end = seq_len(nrow(X)) + 1L, gene_id = rownames(X)) |>
    dplyr::bind_cols(tibble::as_tibble(X, .name_repair = "minimal"))
  expression <- file.path(directory, "input.bed")
  weights <- file.path(directory, "weights.tsv")
  output <- file.path(directory, "outputs")
  readr::write_tsv(bed, expression)
  W |> as.data.frame() |> tibble::rownames_to_column("sample_id") |>
    readr::write_tsv(weights)
  messages <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(
    file.path(script_root, "fit", "fit_tca.R"), "--expression", expression,
    "--weights", weights, "--output-dir", output, "--num-cores", "1", "--max-iters", "2"
  )), stdout = TRUE, stderr = TRUE)
  status <- attr(messages, "status")
  testthat::expect_true(is.null(status) || status == 0L, info = paste(messages, collapse = "\n"))
  testthat::expect_true(all(file.exists(file.path(output,
    c("tca_model.rds", "tca_excluded_genes.tsv", "tca_model.log")))))
  model <- readRDS(file.path(output, "tca_model.rds"))
  testthat::expect_equal(model$W, W, tolerance = 1e-12)
  testthat::expect_identical(rownames(model$mus_hat), rownames(X))
  testthat::expect_identical(model$expression_scale, "cpm")
  testthat::expect_false(model$tca_parallel)
  report <- readr::read_tsv(file.path(output, "tca_excluded_genes.tsv"), show_col_types = FALSE)
  testthat::expect_equal(nrow(report), 0L)
  testthat::expect_match(paste(readLines(file.path(output, "tca_model.log")), collapse = "\n"),
                        "stage_complete", fixed = TRUE)
  # Check the fitted-model handoff on the same tiny fixture. This does not
  # start HSPE, prepare-eQTL, or a WDL workflow.
  cleaned <- file.path(directory, "cleaned")
  messages <- system2(file.path(R.home("bin"), "Rscript"), shQuote(c(
    file.path(script_root, "fit", "clean_tca_model.R"),
    "--model", file.path(output, "tca_model.rds"), "--output-dir", cleaned
  )), stdout = TRUE, stderr = TRUE)
  status <- attr(messages, "status")
  testthat::expect_true(is.null(status) || status == 0L, info = paste(messages, collapse = "\n"))
  clean_model <- readRDS(file.path(cleaned, "tca_model_cleaned.rds"))
  aligned <- align_expression_to_tca_model(X, clean_model)
  testthat::expect_identical(rownames(aligned), rownames(clean_model$mus_hat))
  testthat::expect_identical(colnames(aligned), rownames(clean_model$W))
})

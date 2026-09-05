#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "bootstrap.R"))

tryCatch({
  options <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--model", type = "character", help = "Original fitted TCA model RDS."),
    optparse::make_option("--reuse-model", dest = "reuse_model", action = "store_true", default = FALSE,
                          help = "Validate a supplied CPM model; accept prior numerical cleanup."),
    optparse::make_option("--output-dir", dest = "output_dir", type = "character",
                          help = "Directory for the cleaned model and numerical exclusions.")
  )))
  for (name in c("model", "output_dir")) {
    if (is.null(options[[name]]) || !nzchar(options[[name]])) {
      stop(sprintf("Missing required option: %s", name), call. = FALSE)
    }
  }
  message(sprintf("stage=clean_tca_model event=start utc_time=%s", tca_utc_time()))
  original <- readRDS(options$model)
  result <- if (options$reuse_model) prepare_restart_model(original) else clean_tca_model(original)
  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  model_path <- file.path(options$output_dir, "tca_model_cleaned.rds")
  report_path <- file.path(options$output_dir, "tca_numerical_excluded_genes.tsv")
  saveRDS(result$model, model_path)
  write_numeric_matrix(result$model$W, file.path(options$output_dir, "tca_weights.tsv"), "sample_id")
  readr::write_tsv(result$report, report_path, na = "")
  message(sprintf(
    paste0("stage=clean_tca_model event=complete input_genes=%d retained_genes=%d ",
           "excluded_genes=%d threshold=%.17g model=%s report=%s utc_time=%s"),
    nrow(original$sigmas_hat), nrow(result$model$sigmas_hat), nrow(result$report),
    .Machine$double.eps, model_path, report_path, tca_utc_time()
  ))
}, error = function(error) {
  message(sprintf("stage=clean_tca_model status=failed utc_time=%s message=%s",
                  tca_utc_time(), conditionMessage(error)))
  quit(status = 1L)
})

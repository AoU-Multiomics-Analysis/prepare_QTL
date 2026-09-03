#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "bootstrap.R"))

tca_log_path <- NULL

run_tca_stage <- function() {
  option_list <- list(
    optparse::make_option(
      "--expression",
      type = "character",
      help = "Coordinate-preserving BED of non-negative linear CPM values."
    ),
    optparse::make_option(
      "--weights",
      type = "character",
      help = "Sample-by-group TCA weight TSV with sample_id first column."
    ),
    optparse::make_option(
      "--covariates",
      type = "character",
      default = NULL,
      help = "Optional sample-by-covariate TSV with sample_id first column."
    ),
    optparse::make_option(
      "--num-cores",
      dest = "num_cores",
      type = "integer",
      default = 1L,
      help = "Number of CPU cores for the cohort-wide TCA fit."
    ),
    optparse::make_option(
      "--parallel",
      action = "store_true",
      default = FALSE,
      help = "Enable parallel execution in TCA."
    ),
    optparse::make_option(
      "--max-iters",
      dest = "max_iters",
      type = "integer",
      default = 10L,
      help = "Maximum number of TCA optimization iterations."
    ),
    optparse::make_option(
      "--random-seed",
      dest = "random_seed",
      type = "integer",
      default = 20260901L,
      help = "Random seed for the cohort-wide TCA fit."
    ),
    optparse::make_option(
      "--output-dir",
      dest = "output_dir",
      type = "character",
      help = "Directory for the model and excluded-gene report."
    )
  )
  options <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  required_options <- c("expression", "weights", "output_dir")
  missing_options <- required_options[vapply(
    options[required_options],
    function(value) is.null(value) || !nzchar(value),
    logical(1)
  )]
  if (length(missing_options) > 0L) {
    stop(
      sprintf("Missing required options: %s", paste(missing_options, collapse = ", ")),
      call. = FALSE
    )
  }
  tca_parallel <- validate_boolean_flag(options$parallel, "parallel")

  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  output_paths <- list(
    model = file.path(options$output_dir, "tca_model.rds"),
    model_log = file.path(options$output_dir, "tca_model.log"),
    excluded_genes = file.path(options$output_dir, "tca_excluded_genes.tsv")
  )
  tca_log_path <<- output_paths$model_log
  append_tca_log(
    tca_log_path,
    sprintf(
      "stage=tca event=stage_start scale=cpm output_dir=%s",
      normalizePath(options$output_dir)
    )
  )
  message(sprintf(
    "stage=tca utc_start=%s scale=cpm",
    tca_utc_time()
  ))

  expression <- read_expression_bed(options$expression)
  X <- make_tca_expression(expression)
  W <- read_numeric_matrix(options$weights, "sample_id")
  C2 <- if (is.null(options$covariates) || !nzchar(options$covariates)) {
    NULL
  } else {
    read_numeric_matrix(options$covariates, "sample_id")
  }
  dimension_message <- sprintf(
    paste0(
      "stage=tca input_dimensions=genes:%d samples:%d groups:%d covariates:%d scale=cpm"
    ),
    nrow(X),
    ncol(X),
    ncol(W),
    if (is.null(C2)) 0L else ncol(C2)
  )
  message(dimension_message)
  append_tca_log(tca_log_path, dimension_message)
  settings_message <- sprintf(
    "stage=tca settings=num_cores:%d parallel:%s max_iters:%d random_seed:%d",
    options$num_cores,
    tolower(as.character(tca_parallel)),
    options$max_iters,
    options$random_seed
  )
  message(settings_message)
  append_tca_log(tca_log_path, settings_message)
  output_message <- sprintf(
    "stage=tca output_paths=%s",
    paste(unlist(output_paths, use.names = FALSE), collapse = ",")
  )
  message(output_message)
  append_tca_log(tca_log_path, output_message)

  result <- fit_tca_stage(
    X = X,
    W = W,
    C2 = C2,
    num_cores = options$num_cores,
    parallel = tca_parallel,
    max_iters = options$max_iters,
    random_seed = options$random_seed,
    log_file = output_paths$model_log
  )
  result$model$expression_scale <- "cpm"
  result$model$tca_parallel <- tca_parallel
  saveRDS(result$model, output_paths$model)
  readr::write_tsv(result$excluded_genes, output_paths$excluded_genes, na = "")
  complete_message <- sprintf(
    paste0(
      "stage=tca event=stage_complete output_dimensions=genes:%d samples:%d ",
      "retained_groups:%d excluded_constant_genes:%d scale=cpm"
    ),
    nrow(result$X),
    ncol(result$X),
    ncol(W),
    nrow(result$excluded_genes)
  )
  message(sprintf("%s utc_complete=%s", complete_message, tca_utc_time()))
  append_tca_log(tca_log_path, complete_message)
}

tryCatch(
  run_tca_stage(),
  error = function(error) {
    error_message <- sprintf(
      "stage=tca status=failed utc_time=%s message=%s",
      tca_utc_time(),
      conditionMessage(error)
    )
    message(error_message)
    if (!is.null(tca_log_path)) {
      append_tca_log(tca_log_path, error_message)
    }
    quit(status = 1L)
  }
)

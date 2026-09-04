#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]), fixed = TRUE)
source(file.path(dirname(normalizePath(script_path)), "bootstrap.R"))

tryCatch({
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) != 2L) stop("Usage: filter_cell_type_beds.R CONFIG_JSON OUTPUT_DIR", call. = FALSE)
  config <- jsonlite::read_json(arguments[[1L]], simplifyVector = TRUE)
  required <- c("inventory", "bed_paths", "reference_summary", "min_mean_log2_cpm1", "residual_cutoff")
  if (!all(required %in% names(config))) stop("Filter config is missing required keys", call. = FALSE)
  inventory <- readr::read_tsv(config$inventory, show_col_types = FALSE, progress = FALSE,
                               name_repair = "minimal")
  reference <- if (is.null(config$reference_summary)) NULL else
    readr::read_tsv(config$reference_summary, show_col_types = FALSE, progress = FALSE)
  chunk_size <- if (is.null(config$chunk_size)) 256L else as.integer(config$chunk_size)
  message(sprintf("stage=filter_cell_type_beds start_time=%s", tensor_utc_time()))
  filter_cell_type_beds(inventory, unlist(config$bed_paths, use.names = FALSE), arguments[[2L]],
    reference_summary = reference, min_mean_log2_cpm1 = config$min_mean_log2_cpm1,
    residual_cutoff = config$residual_cutoff, chunk_size = chunk_size)
  message(sprintf("stage=filter_cell_type_beds completion_time=%s", tensor_utc_time()))
}, error = function(error) {
  message(sprintf("stage=filter_cell_type_beds status=failed utc_time=%s message=%s",
                  tensor_utc_time(), conditionMessage(error)))
  quit(status = 1L)
})

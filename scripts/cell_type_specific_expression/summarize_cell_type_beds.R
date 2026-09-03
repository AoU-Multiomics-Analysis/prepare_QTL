#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]), fixed = TRUE)
source(file.path(dirname(normalizePath(script_path)), "bootstrap.R"))

tryCatch({
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) != 3L) {
    stop("Usage: summarize_cell_type_beds.R INVENTORY BED_PATHS OUTPUT_TSV_GZ", call. = FALSE)
  }
  message(sprintf("stage=summarize_cell_type_beds start_time=%s scale=cpm", tensor_utc_time()))
  inventory <- readr::read_tsv(
    arguments[[1]], show_col_types = FALSE, progress = FALSE,
    name_repair = "minimal"
  )
  bed_paths <- readLines(arguments[[2]], warn = FALSE)
  summarize_cell_type_beds(inventory, bed_paths, arguments[[3]])
}, error = function(error) {
  message(sprintf(
    "stage=summarize_cell_type_beds status=failed utc_time=%s message=%s",
    tensor_utc_time(), conditionMessage(error)
  ))
  quit(status = 1L)
})

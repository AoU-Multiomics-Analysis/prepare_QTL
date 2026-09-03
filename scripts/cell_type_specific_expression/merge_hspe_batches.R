#!/usr/bin/env Rscript
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", file_arg[[1]]), fixed = TRUE)
source(file.path(dirname(normalizePath(script_path)), "bootstrap.R"))

tryCatch({
  args <- commandArgs(TRUE)
  if (length(args) != 3L) stop("Usage: merge_hspe_batches.R PREPARED RESULT_PATHS OUTPUT_DIR", call. = FALSE)
  shared <- readRDS(args[1])
  paths <- readLines(args[2], warn = FALSE)
  message(sprintf("stage=merge_hspe_batches input_paths=%s batches=%d", args[1], length(paths)))
  merged <- merge_hspe_batches(purrr::map(paths, readRDS), shared)
  dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
  write_numeric_matrix(merged$proportions, file.path(args[3], "hspe_proportions.tsv"), "sample_id")
  readr::write_tsv(merged$diagnostics, file.path(args[3], "hspe_sample_diagnostics.tsv"))
  jsonlite::write_json(merged$metadata, file.path(args[3], "hspe_metadata.json"),
                       auto_unbox = TRUE, pretty = TRUE)
  message(sprintf("stage=merge_hspe_batches dimensions=samples:%d,cell_types:%d outputs=%s",
                  nrow(merged$proportions), ncol(merged$proportions), args[3]))
}, error = function(error) {
  message(sprintf("stage=merge_hspe_batches status=failed message=%s", conditionMessage(error)))
  quit(status = 1L)
})

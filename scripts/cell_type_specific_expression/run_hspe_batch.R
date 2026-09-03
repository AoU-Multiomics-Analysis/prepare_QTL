#!/usr/bin/env Rscript
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", file_arg[[1]]), fixed = TRUE)
source(file.path(dirname(normalizePath(script_path)), "bootstrap.R"))

tryCatch({
  args <- commandArgs(TRUE)
  if (length(args) != 3L) stop("Usage: run_hspe_batch.R PREPARED BATCH OUTPUT_DIR", call. = FALSE)
  shared <- readRDS(args[1])
  batch <- readRDS(args[2])
  message(sprintf("stage=hspe_batch dimensions=samples:%d,markers:%d input_paths=%s,%s",
                  nrow(batch$Y), ncol(batch$Y), args[1], args[2]))
  result <- fit_hspe_batch(batch, shared)
  dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
  saveRDS(result, file.path(args[3], "hspe_batch_result.rds"))
  message(sprintf("stage=hspe_batch outputs=%s", args[3]))
}, error = function(error) {
  message(sprintf("stage=hspe_batch status=failed message=%s", conditionMessage(error)))
  quit(status = 1L)
})

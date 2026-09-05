#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]), fixed = TRUE)
source(file.path(dirname(normalizePath(script_path)), "bootstrap.R"))

tryCatch({
  options <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--inventory", type = "character"),
    optparse::make_option("--bed-list", dest = "bed_list", type = "character"),
    optparse::make_option("--reference-summary", dest = "reference_summary", type = "character"),
    optparse::make_option("--min-mean-log2-cpm1", dest = "min_mean_log2_cpm1", type = "double", default = 0.01),
    optparse::make_option("--residual-cutoff", dest = "residual_cutoff", type = "double"),
    optparse::make_option("--chunk-size", dest = "chunk_size", type = "integer", default = 256L),
    optparse::make_option("--output-dir", dest = "output_dir", type = "character")
  )))
  for (name in c("inventory", "bed_list", "output_dir")) {
    if (is.null(options[[name]]) || !nzchar(options[[name]])) {
      stop(sprintf("Missing required option: %s", name), call. = FALSE)
    }
  }
  bed_paths <- readLines(options$bed_list, warn = FALSE)
  if (length(bed_paths) == 0L || any(!nzchar(bed_paths))) {
    stop("BED list must contain one non-empty localized path per line", call. = FALSE)
  }
  inventory <- readr::read_tsv(options$inventory, show_col_types = FALSE, progress = FALSE,
                               name_repair = "minimal")
  reference <- if (is.null(options$reference_summary)) NULL else
    readr::read_tsv(options$reference_summary, show_col_types = FALSE, progress = FALSE)
  message(sprintf("stage=filter_cell_type_beds start_time=%s", tensor_utc_time()))
  filter_cell_type_beds(inventory, bed_paths, options$output_dir,
    reference_summary = reference, min_mean_log2_cpm1 = options$min_mean_log2_cpm1,
    residual_cutoff = options$residual_cutoff, chunk_size = options$chunk_size)
  message(sprintf("stage=filter_cell_type_beds completion_time=%s", tensor_utc_time()))
}, error = function(error) {
  message(sprintf("stage=filter_cell_type_beds status=failed utc_time=%s message=%s",
                  tensor_utc_time(), conditionMessage(error)))
  quit(status = 1L)
})

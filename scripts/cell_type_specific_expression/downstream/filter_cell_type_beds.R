#!/usr/bin/env Rscript

# BED inputs use linear CPM.
# Expression thresholds use the mean of log2(CPM + 1) across samples.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]), fixed = TRUE)
source(file.path(dirname(dirname(normalizePath(script_path))), "bootstrap.R"))
source(file.path(dirname(normalizePath(script_path)), "filter_bed_reports.R"))

tryCatch({
  options <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--inventory", type = "character"),
    optparse::make_option("--bed-list", dest = "bed_list", type = "character"),
    optparse::make_option("--single-bed", dest = "single_bed", type = "character"),
    optparse::make_option("--reference-summary", dest = "reference_summary", type = "character"),
    optparse::make_option("--min-mean-log2-cpm1", dest = "min_mean_log2_cpm1", type = "double", default = 0.01),
    optparse::make_option("--residual-cutoff", dest = "residual_cutoff", type = "double"),
    optparse::make_option("--chunk-size", dest = "chunk_size", type = "integer", default = 256L),
    optparse::make_option("--output-dir", dest = "output_dir", type = "character")
  )))
  for (name in c("inventory", "output_dir")) {
    if (is.null(options[[name]]) || !nzchar(options[[name]])) {
      stop(sprintf("Missing required option: %s", name), call. = FALSE)
    }
  }
  if (is.null(options$bed_list) == is.null(options$single_bed)) {
    stop("Supply exactly one of --bed-list or --single-bed", call. = FALSE)
  }
  require_filter_files(c(options$inventory, options$bed_list, options$single_bed,
                         options$reference_summary))
  bed_paths <- if (is.null(options$single_bed)) readLines(options$bed_list, warn = FALSE) else options$single_bed
  if (length(bed_paths) == 0L || any(!nzchar(bed_paths))) {
    stop("BED list must contain one non-empty localized path per line", call. = FALSE)
  }
  inventory <- readr::read_tsv(options$inventory, show_col_types = FALSE, progress = FALSE,
                               name_repair = "minimal")
  require_filter_files(bed_paths)
  if (!is.null(options$single_bed)) inventory <- select_bed_inventory(inventory, options$single_bed)
  reference <- if (is.null(options$reference_summary)) NULL else
    readr::read_tsv(options$reference_summary, show_col_types = FALSE, progress = FALSE)
  message(sprintf("stage=filter_cell_type_beds start_time=%s", tensor_utc_time()))
  result <- filter_cell_type_beds(inventory, bed_paths, options$output_dir,
    reference_summary = reference, min_mean_log2_cpm1 = options$min_mean_log2_cpm1,
    residual_cutoff = options$residual_cutoff, chunk_size = options$chunk_size,
    make_plots = is.null(options$single_bed))
  writeLines(result$samples, file.path(options$output_dir, "sample_ids.txt"))
  message(sprintf("stage=filter_cell_type_beds completion_time=%s", tensor_utc_time()))
}, error = function(error) {
  message(sprintf("stage=filter_cell_type_beds status=failed utc_time=%s message=%s",
                  tensor_utc_time(), conditionMessage(error)))
  quit(status = 1L)
})

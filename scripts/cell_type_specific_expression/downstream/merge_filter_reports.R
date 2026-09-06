#!/usr/bin/env Rscript
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]), fixed = TRUE)
source(file.path(dirname(dirname(normalizePath(script_path))), "bootstrap.R"))
source(file.path(dirname(normalizePath(script_path)), "filter_bed_reports.R"))

tryCatch({
  options <- optparse::parse_args(optparse::OptionParser(option_list = list(
    optparse::make_option("--inventory", type = "character"),
    optparse::make_option("--inventories", type = "character"),
    optparse::make_option("--comparisons", type = "character"),
    optparse::make_option("--metrics", type = "character"),
    optparse::make_option("--samples", type = "character"),
    optparse::make_option("--logs", type = "character"),
    optparse::make_option("--post-residual", dest = "post_residual", action = "store_true", default = FALSE),
    optparse::make_option("--output-dir", dest = "output_dir", type = "character")
  )))
  for (name in c("inventory", "inventories", "comparisons", "metrics", "samples", "logs", "output_dir")) {
    if (is.null(options[[name]]) || !nzchar(options[[name]])) stop("Missing required option: ", name)
  }
  require_filter_files(unlist(options[c("inventory", "inventories", "comparisons", "metrics", "samples", "logs")]))
  lists <- purrr::map(options[c("inventories", "comparisons", "metrics", "samples", "logs")],
    ~ readLines(.x, warn = FALSE))
  purrr::walk(lists, require_filter_files)
  message(sprintf("stage=merge_filter_reports start_time=%s", tensor_utc_time()))
  original <- readr::read_tsv(options$inventory, show_col_types = FALSE, progress = FALSE)
  merge_filter_reports(original, lists$inventories, lists$comparisons, lists$metrics,
    lists$samples, options$output_dir, options$post_residual)
  purrr::walk(lists$logs, ~ cat(readLines(.x, warn = FALSE), sep = "\n"))
  message(sprintf("stage=merge_filter_reports completion_time=%s", tensor_utc_time()))
}, error = function(error) {
  message(sprintf("stage=merge_filter_reports status=failed utc_time=%s message=%s",
    tensor_utc_time(), conditionMessage(error)))
  quit(status = 1L)
})

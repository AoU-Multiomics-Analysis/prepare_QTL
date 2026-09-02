#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "bootstrap.R"))

scatter_log_path <- NULL

append_scatter_log <- function(log_file, message_text) {
  if (is.null(log_file)) {
    return(invisible(NULL))
  }
  cat(
    sprintf("utc_time=%s %s\n", tensor_utc_time(), message_text),
    file = log_file,
    append = TRUE
  )
  invisible(log_file)
}

write_scatter_contract <- function(contract, output_dir) {
  output_files <- list(
    cell_types = file.path(output_dir, "cell_types.txt"),
    cell_type_slugs = file.path(output_dir, "cell_type_slugs.txt"),
    expression_beds = file.path(output_dir, "expression_beds.txt"),
    output_prefixes = file.path(output_dir, "output_prefixes.txt")
  )
  writeLines(contract$cell_type, output_files$cell_types)
  writeLines(contract$cell_type_slug, output_files$cell_type_slugs)
  writeLines(contract$expression_bed, output_files$expression_beds)
  writeLines(contract$output_prefix, output_files$output_prefixes)
  unname(unlist(output_files, use.names = FALSE))
}

run_prepare_scatter_inputs <- function() {
  option_list <- list(
    optparse::make_option(
      "--inventory",
      type = "character",
      help = "Cell-type BED inventory TSV."
    ),
    optparse::make_option(
      "--bed-paths",
      dest = "bed_paths",
      type = "character",
      help = "Newline-delimited localized cell-type BED paths."
    ),
    optparse::make_option(
      "--output-prefix",
      dest = "output_prefix",
      type = "character",
      default = NULL,
      help = "Shared safe output token for direct command-line use."
    ),
    optparse::make_option(
      "--output-prefix-file",
      dest = "output_prefix_file",
      type = "character",
      default = NULL,
      help = "File that contains one shared safe output token."
    ),
    optparse::make_option(
      "--output-dir",
      dest = "output_dir",
      type = "character",
      help = "Directory for aligned scatter metadata files."
    ),
    optparse::make_option(
      "--log-file",
      dest = "log_file",
      type = "character",
      default = NULL,
      help = "Optional scatter validation log path."
    )
  )
  options <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  required_options <- c("inventory", "bed_paths", "output_dir")
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
  prefix_sources <- c(
    !is.null(options[["output_prefix"]]),
    !is.null(options[["output_prefix_file"]])
  )
  if (sum(prefix_sources) != 1L) {
    stop(
      "Supply exactly one of --output-prefix or --output-prefix-file",
      call. = FALSE
    )
  }
  output_prefix <- if (!is.null(options[["output_prefix_file"]])) {
    prefix_text <- readr::read_file(options[["output_prefix_file"]])
    sub("\n\\z", "", prefix_text, perl = TRUE)
  } else {
    options[["output_prefix"]]
  }

  dir.create(options$output_dir, recursive = TRUE, showWarnings = FALSE)
  scatter_log_path <<- if (is.null(options$log_file)) {
    file.path(options$output_dir, "prepare_scatter_inputs.log")
  } else {
    options$log_file
  }
  inventory <- readr::read_tsv(
    options$inventory,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )
  bed_paths <- readLines(options$bed_paths, warn = FALSE)
  start_message <- sprintf("stage=prepare_scatter_inputs start_time=%s", tensor_utc_time())
  dimensions_message <- sprintf(
    "stage=prepare_scatter_inputs input_dimensions=inventory_rows:%d,bed_paths:%d",
    nrow(inventory),
    length(bed_paths)
  )
  message(start_message)
  message(dimensions_message)
  append_scatter_log(scatter_log_path, start_message)
  append_scatter_log(scatter_log_path, dimensions_message)

  contract <- prepare_scatter_contract(
    inventory = inventory,
    bed_paths = bed_paths,
    output_prefix = output_prefix
  )
  output_files <- write_scatter_contract(contract, options$output_dir)
  completion_message <- sprintf(
    paste0(
      "stage=prepare_scatter_inputs validated_cell_count=%d ",
      "output_paths=%s completion_time=%s"
    ),
    nrow(contract),
    paste(output_files, collapse = ","),
    tensor_utc_time()
  )
  message(completion_message)
  append_scatter_log(scatter_log_path, completion_message)
}

tryCatch(
  run_prepare_scatter_inputs(),
  error = function(error) {
    error_message <- sprintf(
      "stage=prepare_scatter_inputs status=failed utc_time=%s message=%s",
      tensor_utc_time(),
      conditionMessage(error)
    )
    message(error_message)
    if (!is.null(scatter_log_path)) {
      append_scatter_log(scatter_log_path, error_message)
    }
    quit(status = 1L)
  }
)

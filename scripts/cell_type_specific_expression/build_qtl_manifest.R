#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "bootstrap.R"))

run_qtl_manifest_stage <- function() {
  option_list <- list(
    optparse::make_option("--cell-types", dest = "cell_types", type = "character"),
    optparse::make_option(
      "--cell-type-slugs",
      dest = "cell_type_slugs",
      type = "character"
    ),
    optparse::make_option("--int-beds", dest = "int_beds", type = "character"),
    optparse::make_option("--scaled-beds", dest = "scaled_beds", type = "character"),
    optparse::make_option("--int-pcs", dest = "int_pcs", type = "character"),
    optparse::make_option(
      "--int-pcs-all",
      dest = "int_pcs_all",
      type = "character"
    ),
    optparse::make_option(
      "--scaled-pcs",
      dest = "scaled_pcs",
      type = "character"
    ),
    optparse::make_option(
      "--scaled-pcs-all",
      dest = "scaled_pcs_all",
      type = "character"
    ),
    optparse::make_option(
      "--int-covariates",
      dest = "int_covariates",
      type = "character"
    ),
    optparse::make_option(
      "--scaled-covariates",
      dest = "scaled_covariates",
      type = "character"
    ),
    optparse::make_option(
      "--int-outliers",
      dest = "int_outliers",
      type = "character"
    ),
    optparse::make_option(
      "--scaled-outliers",
      dest = "scaled_outliers",
      type = "character"
    ),
    optparse::make_option("--output", type = "character")
  )
  options <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
  required_options <- c(
    "cell_types", "cell_type_slugs", "int_beds", "scaled_beds", "int_pcs",
    "int_pcs_all", "scaled_pcs", "scaled_pcs_all", "int_covariates",
    "scaled_covariates", "int_outliers", "scaled_outliers", "output"
  )
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

  input_option_names <- setdiff(required_options, "output")
  inputs <- purrr::map(options[input_option_names], readLines, warn = FALSE)
  message(sprintf(
    "stage=build_qtl_manifest start_time=%s",
    tensor_utc_time()
  ))
  message(sprintf(
    "stage=build_qtl_manifest input_dimensions=cell_types:%d file_categories:%d",
    length(inputs$cell_types),
    length(input_option_names) - 2L
  ))
  message(sprintf(
    "stage=build_qtl_manifest input_paths=%s",
    paste(unlist(options[input_option_names], use.names = FALSE), collapse = ",")
  ))

  manifest <- do.call(build_cell_type_qtl_manifest, inputs)
  dir.create(dirname(options$output), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(manifest, options$output)
  message(sprintf(
    paste0(
      "stage=build_qtl_manifest validated_cell_count=%d manifest_path=%s ",
      "completion_time=%s"
    ),
    nrow(manifest),
    options$output,
    tensor_utc_time()
  ))
}

tryCatch(
  run_qtl_manifest_stage(),
  error = function(error) {
    message(sprintf(
      "stage=build_qtl_manifest status=failed utc_time=%s message=%s",
      tensor_utc_time(), conditionMessage(error)
    ))
    quit(status = 1L)
  }
)

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
    optparse::make_option("--source-beds", dest = "source_beds", type = "character"),
    optparse::make_option("--source-inventory", dest = "source_inventory", type = "character"),
    optparse::make_option("--filtered-beds", dest = "filtered_beds", type = "character"),
    optparse::make_option("--filtered-inventory", dest = "filtered_inventory", type = "character"),
    optparse::make_option("--negative-summary", dest = "negative_summary", type = "character"),
    optparse::make_option("--gene-comparison", dest = "gene_comparison", type = "character"),
    optparse::make_option("--filter-metrics", dest = "filter_metrics", type = "character"),
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
  # Each CLI input is a JSON array serialized inside the WDL task. Keep paths
  # as metadata; never open the upstream cloud files from this script.
  inputs <- purrr::map(options[input_option_names], jsonlite::read_json, simplifyVector = TRUE)
  filter_option_names <- c(
    "source_beds", "source_inventory", "filtered_beds", "filtered_inventory",
    "negative_summary", "gene_comparison", "filter_metrics"
  )
  provided_filter_options <- !vapply(options[filter_option_names], is.null, logical(1))
  if (any(provided_filter_options) && !all(provided_filter_options)) {
    stop("Filter manifest options must be supplied together", call. = FALSE)
  }
  if (all(provided_filter_options)) {
    source_inventory <- readr::read_tsv(
      options$source_inventory, show_col_types = FALSE, progress = FALSE
    )
    filtered_inventory <- readr::read_tsv(
      options$filtered_inventory, show_col_types = FALSE, progress = FALSE
    )
    inputs$source_beds <- jsonlite::read_json(options$source_beds, simplifyVector = TRUE)
    inputs$source_bed_slugs <- source_inventory$slug[match(
      basename(inputs$source_beds), basename(source_inventory$path)
    )]
    inputs$filtered_beds <- jsonlite::read_json(options$filtered_beds, simplifyVector = TRUE)
    inputs$filtered_bed_slugs <- filtered_inventory$slug[match(
      basename(inputs$filtered_beds), basename(filtered_inventory$path)
    )]
    inputs$negative_summary <- jsonlite::read_json(options$negative_summary, simplifyVector = TRUE)
    inputs$gene_comparison <- jsonlite::read_json(options$gene_comparison, simplifyVector = TRUE)
    inputs$filter_metrics <- jsonlite::read_json(options$filter_metrics, simplifyVector = TRUE)
  }
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

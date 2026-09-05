#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]))
source(file.path(dirname(script_path), "bootstrap.R"))

manifest_log_path <- NULL

run_build_manifest <- function() {
  message(sprintf("stage=manifest utc_start=%s", tensor_utc_time()))
  option_list <- list(
    optparse::make_option(
      "--outputs",
      type = "character",
      help = "TSV output inventory from direct BED export."
    ),
    optparse::make_option(
      "--export-qc-summary",
      dest = "export_qc_summary",
      type = "character",
      help = "Direct BED export QC summary TSV."
    ),
    optparse::make_option(
      "--original-proportions",
      dest = "original_proportions",
      type = "character",
      help = "Sample-by-LM22 proportion TSV."
    ),
    optparse::make_option(
      "--combined-proportions",
      dest = "combined_proportions",
      type = "character",
      help = "Sample-by-combined-group proportion TSV."
    ),
    optparse::make_option(
      "--tca-weights",
      dest = "tca_weights",
      type = "character",
      help = "Sample-by-retained-group TCA weight TSV."
    ),
    optparse::make_option(
      "--filter-report",
      dest = "filter_report",
      type = "character",
      help = "Cell-group filter and zero-adjustment report TSV."
    ),
    optparse::make_option(
      "--model",
      type = "character",
      help = "Fitted TCA model RDS."
    ),
    optparse::make_option(
      "--model-log",
      dest = "model_log",
      type = "character",
      help = "TCA model log with internal iteration records."
    ),
    optparse::make_option(
      "--hspe-metadata",
      dest = "hspe_metadata",
      type = "character",
      default = NULL,
      help = "Optional hspe metadata JSON with LM22 QC."
    ),
    optparse::make_option(
      "--tca-version",
      dest = "tca_version",
      type = "character",
      default = "1.2.1",
      help = "TCA package version."
    ),
    optparse::make_option(
      "--proportion-mode",
      dest = "proportion_mode",
      type = "character",
      help = "Proportion source: hspe or precomputed."
    ),
    optparse::make_option(
      "--log2-pseudocount",
      dest = "log2_pseudocount",
      type = "double",
      help = "Pseudocount used before the log2 transform."
    ),
    optparse::make_option(
      "--min-lm22-overlap",
      dest = "min_lm22_overlap",
      type = "double",
      help = "Minimum required LM22 gene overlap."
    ),
    optparse::make_option(
      "--hspe-marker-fraction",
      dest = "hspe_marker_fraction",
      type = "double",
      help = "Fraction of markers selected for hspe."
    ),
    optparse::make_option(
      "--hspe-marker-method",
      dest = "hspe_marker_method",
      type = "character",
      help = "Marker selection method used by hspe."
    ),
    optparse::make_option(
      "--hspe-quantile-normalize",
      dest = "hspe_quantile_normalize",
      type = "character",
      help = "Whether hspe quantile normalization was enabled."
    ),
    optparse::make_option(
      "--group-mean-threshold",
      dest = "group_mean_threshold",
      type = "double",
      help = "Mean-proportion threshold used to retain major groups."
    ),
    optparse::make_option(
      "--zero-floor",
      dest = "zero_floor",
      type = "double",
      help = "Replacement value used for exact zero proportions."
    ),
    optparse::make_option(
      "--tca-max-iters",
      dest = "tca_max_iters",
      type = "integer",
      help = "Maximum TCA iteration count."
    ),
    optparse::make_option(
      "--tca-parallel",
      dest = "tca_parallel",
      type = "character",
      help = "Whether TCA parallel execution was enabled."
    ),
    optparse::make_option(
      "--gene-types",
      dest = "gene_types",
      type = "character",
      help = "Comma-separated GTF gene types retained for deconvolution."
    ),
    optparse::make_option(
      "--random-seed",
      dest = "random_seed",
      type = "integer",
      help = "Random seed used for TCA."
    ),
    optparse::make_option(
      "--scale",
      type = "character",
      help = "Expression scale used by the model."
    ),
    optparse::make_option(
      "--effective-parameters-output",
      dest = "effective_parameters_output",
      type = "character",
      help = "Output JSON path for the effective parameters."
    ),
    optparse::make_option(
      "--container-image",
      dest = "container_image",
      type = "character",
      help = "Container image name and immutable digest."
    ),
    optparse::make_option(
      "--output",
      type = "character",
      help = "Output manifest JSON path."
    ),
    optparse::make_option(
      "--qc-output",
      dest = "qc_output",
      type = "character",
      help = "Final pipeline QC summary TSV path."
    ),
    optparse::make_option(
      "--log-file",
      dest = "log_file",
      type = "character",
      default = NULL,
      help = "Manifest log path."
    )
  )
  options <- optparse::parse_args(optparse::OptionParser(
    option_list = option_list
  ))
  required_options <- c(
    "outputs", "export_qc_summary",
    "tca_weights", "model",
    "proportion_mode", "log2_pseudocount", "min_lm22_overlap",
    "hspe_marker_fraction", "hspe_marker_method",
    "hspe_quantile_normalize", "group_mean_threshold", "zero_floor",
    "tca_max_iters", "tca_parallel", "gene_types", "random_seed", "scale",
    "effective_parameters_output",
    "container_image", "output", "qc_output"
  )
  restart <- identical(options$proportion_mode, "precomputed_model")
  if (!restart) {
    required_options <- c(required_options, "original_proportions", "combined_proportions",
                          "filter_report", "model_log")
  } else {
    required_options <- setdiff(required_options, c(
      "log2_pseudocount", "min_lm22_overlap", "hspe_marker_fraction", "hspe_marker_method",
      "hspe_quantile_normalize", "group_mean_threshold", "zero_floor", "tca_max_iters",
      "gene_types", "random_seed"
    ))
  }
  missing_options <- required_options[vapply(
    options[required_options],
    function(value) is.null(value) || !nzchar(value),
    logical(1)
  )]
  if (length(missing_options) > 0L) {
    stop(
      sprintf(
        "Missing required options: %s",
        paste(missing_options, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  purrr::walk(
    unique(c(
      dirname(options$output),
      dirname(options$qc_output),
      dirname(options$effective_parameters_output)
    )),
    ~ dir.create(.x, recursive = TRUE, showWarnings = FALSE)
  )
  manifest_log_path <<- if (is.null(options$log_file)) {
    paste0(options$output, ".log")
  } else {
    options$log_file
  }
  outputs <- readr::read_tsv(
    options$outputs,
    col_types = readr::cols(
      logical_name = readr::col_character(),
      path = readr::col_character(),
      n_genes = readr::col_integer(),
      n_samples = readr::col_integer(),
      scale = readr::col_character(),
      cell_group = readr::col_character()
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
  export_qc_summary <- readr::read_tsv(
    options$export_qc_summary,
    show_col_types = FALSE,
    progress = FALSE
  )
  original_proportions <- if (restart) NULL else read_numeric_matrix(
    options$original_proportions,
    "sample_id"
  )
  combined_proportions <- if (restart) NULL else read_numeric_matrix(
    options$combined_proportions,
    "sample_id"
  )
  tca_weights <- read_numeric_matrix(options$tca_weights, "sample_id")
  filter_report <- if (restart) NULL else readr::read_tsv(
    options$filter_report,
    col_types = readr::cols(
      cell_group = readr::col_character(),
      cohort_mean = readr::col_double(),
      threshold = readr::col_double(),
      retained = readr::col_logical(),
      filter_reason = readr::col_character(),
      zero_count_before = readr::col_integer(),
      zero_floor = readr::col_double()
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
  tca_model <- readRDS(options$model)
  tca_log_lines <- if (restart) character() else readLines(options$model_log, warn = FALSE)
  hspe_metadata <- if (is.null(options$hspe_metadata) ||
      !nzchar(options$hspe_metadata)) {
    NULL
  } else {
    jsonlite::read_json(options$hspe_metadata, simplifyVector = FALSE)
  }
  quantile_normalize <- if (restart) "false" else tolower(options$hspe_quantile_normalize)
  if (!quantile_normalize %in% c("true", "false")) {
    stop(
      "hspe_quantile_normalize must be true or false",
      call. = FALSE
    )
  }
  tca_parallel <- tolower(options$tca_parallel)
  if (!tca_parallel %in% c("true", "false")) {
    stop("tca_parallel must be true or false", call. = FALSE)
  }
  gene_types <- if (restart) character() else parse_gene_types(options$gene_types)
  parameters <- if (restart) {
    list(proportion_mode = "precomputed_model", scale = "cpm",
         tca_parallel = identical(tca_parallel, "true"))
  } else list(
    proportion_mode = options$proportion_mode,
    log2_pseudocount = validate_log2_pseudocount(
      options$log2_pseudocount
    ),
    min_lm22_overlap = options$min_lm22_overlap,
    hspe_marker_fraction = options$hspe_marker_fraction,
    hspe_marker_method = options$hspe_marker_method,
    hspe_quantile_normalize = identical(quantile_normalize, "true"),
    group_mean_threshold = options$group_mean_threshold,
    zero_floor = options$zero_floor,
    tca_max_iters = options$tca_max_iters,
    tca_parallel = identical(tca_parallel, "true"),
    gene_type = gene_types,
    random_seed = options$random_seed,
    scale = options$scale
  )
  dimensions_message <- sprintf(
    paste0(
      "stage=manifest input_dimensions=outputs:%d samples:%d ",
      "lm22_types:%d retained_groups:%d log2_pseudocount:%g gene_types:%s"
    ),
    nrow(outputs),
    nrow(tca_weights),
    if (restart) 0L else ncol(original_proportions),
    ncol(tca_weights),
    if (restart) NA_real_ else parameters$log2_pseudocount,
    paste(gene_types, collapse = ",")
  )
  paths_message <- sprintf(
    "stage=manifest input_paths=%s output_paths=%s",
    options$outputs,
    paste(c(options$output, options$qc_output), collapse = ",")
  )
  message(dimensions_message)
  message(paths_message)
  append_tensor_log(manifest_log_path, dimensions_message)
  append_tensor_log(manifest_log_path, paths_message)

  manifest <- build_output_manifest(
    outputs = outputs,
    tca_version = options$tca_version,
    parameters = parameters,
    container_image = options$container_image
  )
  qc_summary <- if (restart) {
    build_restart_qc_summary(export_qc_summary, tca_weights, tca_model)
  } else build_pipeline_qc_summary(
    export_summary = export_qc_summary,
    original_proportions = original_proportions,
    combined_proportions = combined_proportions,
    tca_weights = tca_weights,
    filter_report = filter_report,
    tca_model = tca_model,
    tca_log_lines = tca_log_lines,
    hspe_metadata = hspe_metadata
  )
  write_output_manifest(options$output, manifest)
  jsonlite::write_json(
    parameters,
    options$effective_parameters_output,
    auto_unbox = TRUE,
    pretty = TRUE,
    digits = NA
  )
  readr::write_tsv(qc_summary, options$qc_output, na = "")
  complete_message <- sprintf(
    paste0(
      "stage=manifest event=stage_complete outputs=%d qc_metrics=%d ",
      "output_paths=%s"
    ),
    length(manifest$outputs),
    nrow(qc_summary),
    paste(
      c(
        options$output,
        options$qc_output,
        options$effective_parameters_output
      ),
      collapse = ","
    )
  )
  message(sprintf("%s utc_complete=%s", complete_message, tensor_utc_time()))
  append_tensor_log(manifest_log_path, complete_message)
}

tryCatch(
  run_build_manifest(),
  error = function(error) {
    error_message <- sprintf(
      "stage=manifest status=failed utc_time=%s message=%s",
      tensor_utc_time(),
      conditionMessage(error)
    )
    message(error_message)
    if (!is.null(manifest_log_path)) {
      append_tensor_log(manifest_log_path, error_message)
    }
    quit(status = 1L)
  }
)

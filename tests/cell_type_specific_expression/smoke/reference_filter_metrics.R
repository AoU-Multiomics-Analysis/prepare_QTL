read_reference_filter_metrics <- function(path) {
  readr::read_tsv(
    path,
    col_types = readr::cols(
      cell_type = readr::col_character(),
      slug = readr::col_character(),
      comparison_status = readr::col_character(),
      deconvolution_n_samples = readr::col_integer(),
      reference_n_samples = readr::col_integer(),
      n_original = readr::col_integer(),
      n_negative_excluded = readr::col_integer(),
      n_reference_excluded = readr::col_integer(),
      n_residual_excluded = readr::col_integer(),
      n_retained = readr::col_integer(),
      metric_set = readr::col_character(),
      n_genes = readr::col_integer(),
      pearson_r = readr::col_double(),
      spearman_rho = readr::col_double(),
      r_squared = readr::col_double(),
      intercept = readr::col_double(),
      slope = readr::col_double(),
      .default = readr::col_character()
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
}

validate_reference_filter_metrics <- function(metrics, expected_cell_types,
                                              mapped_cell_types,
                                              reference_enabled) {
  required_columns <- c(
    "cell_type", "slug", "comparison_status", "deconvolution_n_samples",
    "reference_n_samples", "n_original", "n_negative_excluded",
    "n_reference_excluded", "n_residual_excluded", "n_retained",
    "metric_set", "n_genes", "pearson_r", "spearman_rho", "r_squared",
    "intercept", "slope"
  )
  if (!inherits(metrics, "data.frame") ||
      !all(required_columns %in% names(metrics))) {
    stop("The reference filter metrics do not contain all documented columns", call. = FALSE)
  }

  baseline <- dplyr::filter(metrics, .data$metric_set == "baseline")
  if (!identical(baseline$cell_type, expected_cell_types)) {
    stop("The reference filter metrics must contain one ordered baseline row per cell type",
         call. = FALSE)
  }
  count_columns <- c(
    "deconvolution_n_samples", "n_original", "n_negative_excluded",
    "n_reference_excluded", "n_residual_excluded", "n_retained", "n_genes"
  )
  metric_columns <- c(
    count_columns, "reference_n_samples", "pearson_r", "spearman_rho",
    "r_squared", "intercept", "slope"
  )
  if (any(!vapply(baseline[metric_columns], is.numeric, logical(1))) ||
      anyNA(baseline$slug) || any(!nzchar(baseline$slug)) ||
      anyDuplicated(baseline$slug) > 0L) {
    stop("The reference filter metrics contain invalid identifiers or metric values",
         call. = FALSE)
  }
  counts <- as.matrix(dplyr::select(baseline, dplyr::all_of(count_columns)))
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts)) ||
      any(baseline$deconvolution_n_samples == 0L) ||
      any(baseline$n_original == 0L) || any(baseline$n_retained == 0L) ||
      any(baseline$n_original != baseline$n_negative_excluded +
            baseline$n_reference_excluded + baseline$n_residual_excluded +
            baseline$n_retained)) {
    stop("The reference filter metrics contain invalid or inconsistent counts", call. = FALSE)
  }

  mapped <- baseline$cell_type %in% mapped_cell_types
  if (isTRUE(reference_enabled)) {
    allowed_mapped_status <- c("available", "insufficient_genes_or_variation")
    available <- mapped & baseline$comparison_status == "available"
    regression_columns <- c(
      "pearson_r", "spearman_rho", "r_squared", "intercept", "slope"
    )
    if (any(!baseline$comparison_status[mapped] %in% allowed_mapped_status) ||
        any(!is.finite(baseline$reference_n_samples[mapped])) ||
        any(baseline$reference_n_samples[mapped] <= 0) ||
        any(baseline$n_genes[available] < 3L) ||
        any(!is.finite(as.matrix(baseline[available, regression_columns, drop = FALSE]))) ||
        any(baseline$comparison_status[!mapped] != "no_reference_cell_type") ||
        any(!is.na(baseline$reference_n_samples[!mapped]))) {
      stop("The reference-enabled metrics have an invalid reference-comparison status",
           call. = FALSE)
    }
  } else if (any(baseline$comparison_status != "reference_not_provided") ||
             any(!is.na(baseline$reference_n_samples)) ||
             any(baseline$n_reference_excluded != 0L) ||
             any(baseline$n_residual_excluded != 0L)) {
    stop("The negative-only metrics must report that a reference was not provided",
         call. = FALSE)
  }
  invisible(metrics)
}

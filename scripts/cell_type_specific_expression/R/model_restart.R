# Restart only from pipeline models with an explicit linear-CPM scale.
prepare_restart_model <- function(model) {
  if (!identical(model$expression_scale, "cpm")) {
    stop("Restart requires a model recorded in linear CPM space", call. = FALSE)
  }
  validate_tca_gene_parameters(model)
  W <- model$W
  if (!is.numeric(W) || any(!is.finite(W)) || any(W < 0) ||
      nrow(W) == 0L || any(abs(rowSums(W) - 1) > 1e-6)) {
    stop("Model weights must be finite nonnegative proportions with row sums of one", call. = FALSE)
  }
  validate_matrix_identifiers(W, "Model sample", "Model cell type")
  for (name in c("C1", "C2")) {
    C <- model[[name]]
    if (!is.null(C) && (!is.matrix(C) || !is.numeric(C) ||
        any(!is.finite(C)) || nrow(C) != nrow(W) ||
        (ncol(C) > 0L && !identical(rownames(C), rownames(W))))) {
      stop(sprintf("Model %s must contain finite covariates in model sample order", name), call. = FALSE)
    }
  }
  if (is.null(model$gene_filter)) {
    result <- clean_tca_model(model)
    result$model$gene_filter$report <- result$report
    return(result)
  }
  filter <- model$gene_filter
  genes <- rownames(model$mus_hat)
  if (!identical(filter$method, "diagonal_variance_reciprocal_condition") ||
      !identical(filter$threshold, .Machine$double.eps) ||
      !is.character(filter$original_gene_ids) || anyNA(filter$original_gene_ids) ||
      anyDuplicated(filter$original_gene_ids) ||
      !is.character(filter$excluded_gene_ids) || anyNA(filter$excluded_gene_ids) ||
      anyDuplicated(filter$excluded_gene_ids) ||
      !all(filter$excluded_gene_ids %in% filter$original_gene_ids) ||
      !identical(setdiff(filter$original_gene_ids, filter$excluded_gene_ids), genes)) {
    stop("Cleaned model gene exclusion record is invalid", call. = FALSE)
  }
  # Check the retained parameters again without altering the existing audit record.
  candidate <- model
  candidate$gene_filter <- NULL
  checked <- clean_tca_model(candidate)
  if (nrow(checked$report) > 0L) {
    stop("Cleaned model still contains numerically excluded genes", call. = FALSE)
  }
  report <- filter$report
  if (is.null(report)) {
    # Older cleaned models did not store the removed parameter values.
    report <- tibble::tibble(
      gene_id = filter$excluded_gene_ids,
      reason = rep("previously_excluded_variance", length(filter$excluded_gene_ids)),
      min_variance = NA_real_, max_variance = NA_real_,
      reciprocal_condition = NA_real_, threshold = .Machine$double.eps
    )
  }
  list(model = model, report = report)
}

align_restart_expression <- function(X, model) {
  if (!identical(colnames(X), rownames(model$W))) {
    stop("Model sample order must match expression sample order exactly", call. = FALSE)
  }
  genes <- rownames(model$mus_hat)
  missing <- setdiff(genes, rownames(X))
  if (length(missing) > 0L) {
    stop(sprintf("BED is missing model genes: %s", paste(head(missing, 10L), collapse = ", ")),
         call. = FALSE)
  }
  X[genes, , drop = FALSE]
}

build_restart_qc_summary <- function(export_summary, weights, model) {
  dplyr::bind_rows(
    tibble::as_tibble(export_summary) |>
      dplyr::mutate(value = as.numeric(.data$value), status = "observed"),
    tibble::tibble(
      metric = c("adjusted_weight_max_row_sum_error", "tca_tau_hat"),
      value = c(max(abs(rowSums(weights) - 1)), model$tau_hat),
      status = "from_supplied_model"
    ),
    tibble::tibble(
      metric = c("input_proportion_max_row_sum_error", "combined_proportion_max_row_sum_error",
                 "normalization_adjustment_max_abs", "zero_values_adjusted",
                 "tca_internal_iterations", "tca_max_internal_iterations", "tca_convergence",
                 "lm22_gene_count", "lm22_cell_type_count", "lm22_value_min",
                 "lm22_value_max", "lm22_value_validation"),
      value = NA_real_, status = "unavailable_model_restart"
    )
  )
}

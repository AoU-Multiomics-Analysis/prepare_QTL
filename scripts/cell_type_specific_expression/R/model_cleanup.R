# TCA 1.2.1 stores these fields by gene. Do not infer fields by row count:
# a sample matrix can have the same row count as the gene matrices.
tca_gene_parameter_fields <- function() {
  c("mus_hat", "sigmas_hat", "deltas_hat", "gammas_hat",
    "deltas_hat_pvals", "gammas_hat_pvals", "gammas_hat_pvals.joint")
}

validate_tca_gene_parameters <- function(model) {
  required <- c("mus_hat", "sigmas_hat", "deltas_hat", "gammas_hat")
  if (!is.list(model) || !all(required %in% names(model))) {
    stop("TCA model is missing required gene parameter matrices", call. = FALSE)
  }
  sigmas <- model$sigmas_hat
  if (!is.matrix(sigmas) || !is.numeric(sigmas) ||
      nrow(sigmas) == 0L || ncol(sigmas) < 2L) {
    stop("TCA sigmas_hat must be a non-empty gene-by-cell-type matrix", call. = FALSE)
  }
  validate_matrix_identifiers(sigmas, "Model gene", "Model cell type")
  fields <- intersect(tca_gene_parameter_fields(), names(model))
  for (field in fields) {
    value <- model[[field]]
    if (!is.matrix(value) || !is.numeric(value) ||
        !identical(rownames(value), rownames(sigmas))) {
      stop(sprintf("TCA %s gene order must match sigmas_hat", field), call. = FALSE)
    }
  }
  if (!identical(dimnames(model$mus_hat), dimnames(sigmas)) ||
      !is.matrix(model$W) || !identical(colnames(model$W), colnames(sigmas))) {
    stop("TCA model cell-type order must match across means, sigmas, and weights", call. = FALSE)
  }
  if (any(!is.finite(sigmas)) || any(sigmas < 0)) {
    stop("TCA fitted standard deviations must be finite and nonnegative", call. = FALSE)
  }
  for (field in c("mus_hat", "deltas_hat", "gammas_hat")) {
    if (any(!is.finite(model[[field]]))) {
      stop(sprintf("TCA %s values must be finite", field), call. = FALSE)
    }
  }
  if (!is.numeric(model$tau_hat) || length(model$tau_hat) != 1L ||
      !is.finite(model$tau_hat) || model$tau_hat <= 0) {
    stop("TCA tau_hat must be finite and positive", call. = FALSE)
  }
  invisible(fields)
}

clean_tca_model <- function(model) {
  fields <- validate_tca_gene_parameters(model)
  if (!is.null(model$gene_filter)) {
    stop("Model already contains a gene filter; supply the original fitted model", call. = FALSE)
  }
  variances <- model$sigmas_hat^2
  if (any(!is.finite(variances))) {
    stop("TCA fitted variances must be finite after squaring", call. = FALSE)
  }
  minimum <- apply(variances, 1L, min)
  maximum <- apply(variances, 1L, max)
  reciprocal_condition <- ifelse(maximum == 0, 0, minimum / maximum)
  threshold <- .Machine$double.eps
  keep <- reciprocal_condition >= threshold
  gene_ids <- rownames(variances)
  report <- tibble::tibble(
    gene_id = gene_ids,
    reason = "computationally_singular_variance",
    min_variance = unname(minimum),
    max_variance = unname(maximum),
    reciprocal_condition = unname(reciprocal_condition),
    threshold = threshold
  ) |>
    dplyr::filter(!keep)
  if (!any(keep)) {
    stop("No genes remain after TCA numerical filtering", call. = FALSE)
  }
  cleaned <- model
  for (field in fields) {
    cleaned[[field]] <- model[[field]][keep, , drop = FALSE]
  }
  cleaned$gene_filter <- list(
    method = "diagonal_variance_reciprocal_condition",
    threshold = threshold,
    original_gene_ids = gene_ids,
    excluded_gene_ids = gene_ids[!keep]
  )
  list(model = cleaned, report = report)
}

align_expression_to_tca_model <- function(X, model) {
  if (!identical(colnames(X), rownames(model$W))) {
    stop("Model sample order must match expression sample order exactly", call. = FALSE)
  }
  model_genes <- rownames(model$mus_hat)
  if (is.null(model_genes) || anyNA(model_genes) || any(!nzchar(model_genes)) ||
      anyDuplicated(model_genes) > 0L) {
    stop("Model gene identifiers must be unique and non-empty", call. = FALSE)
  }
  filter <- model$gene_filter
  if (is.null(filter)) {
    # Legacy complete models are supported; never silently intersect genes.
    if (!identical(rownames(X), model_genes)) {
      stop("Expression gene order must match the fitted model", call. = FALSE)
    }
    return(X)
  }
  validate_tca_gene_parameters(model)
  if (!identical(filter$method, "diagonal_variance_reciprocal_condition") ||
      !identical(filter$threshold, .Machine$double.eps) ||
      !is.character(filter$original_gene_ids) || anyNA(filter$original_gene_ids) ||
      anyDuplicated(filter$original_gene_ids) > 0L ||
      !identical(rownames(X), filter$original_gene_ids)) {
    stop("Expression gene order must match the original genes recorded in the cleaned model", call. = FALSE)
  }
  excluded <- filter$excluded_gene_ids
  if (!is.character(excluded) || anyNA(excluded) || anyDuplicated(excluded) > 0L ||
      !all(excluded %in% filter$original_gene_ids) ||
      !identical(filter$original_gene_ids[!filter$original_gene_ids %in% excluded], model_genes)) {
    stop("Cleaned model gene exclusion record does not match retained genes", call. = FALSE)
  }
  X[model_genes, , drop = FALSE]
}

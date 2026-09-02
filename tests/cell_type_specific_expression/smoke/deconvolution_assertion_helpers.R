workflow_input_proportion_row_sum_tolerance <- 1e-6

require_row_sums_within_tolerance <- function(values, tolerance, label) {
  if (!is.matrix(values) || !is.numeric(values) || nrow(values) == 0L ||
      ncol(values) == 0L || any(!is.finite(values))) {
    stop(sprintf("The %s values must be a finite numeric matrix", label), call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("The row-sum tolerance must be non-negative", call. = FALSE)
  }
  maximum_error <- max(abs(rowSums(values) - 1))
  if (maximum_error > tolerance + .Machine$double.eps) {
    stop(
      sprintf(
        "The %s rows must sum to one within %.3g; maximum error is %.17g",
        label,
        tolerance,
        maximum_error
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

canonical_lm22_group_map <- function() {
  list(
    "B cells" = c("B cells naive", "B cells memory", "Plasma cells"),
    "CD4 T cells" = c(
      "T cells CD4 naive", "T cells CD4 memory resting",
      "T cells CD4 memory activated", "T cells follicular helper",
      "T cells regulatory (Tregs)"
    ),
    "CD8 T cells" = "T cells CD8",
    "Gamma-delta T cells" = "T cells gamma delta",
    "NK cells" = c("NK cells resting", "NK cells activated"),
    "Monocyte/myeloid" = c(
      "Monocytes", "Macrophages M0", "Macrophages M1", "Macrophages M2"
    ),
    "Neutrophils" = "Neutrophils",
    "Eosinophils" = "Eosinophils",
    "Dendritic cells" = c("Dendritic cells resting", "Dendritic cells activated"),
    "Mast cells" = c("Mast cells resting", "Mast cells activated")
  )
}

derive_expected_proportion_outputs <- function(
    proportions,
    mean_threshold,
    zero_floor) {
  group_map <- canonical_lm22_group_map()
  required_cell_types <- unlist(group_map, use.names = FALSE)
  if (!is.matrix(proportions) || !is.numeric(proportions) ||
      is.null(rownames(proportions)) || is.null(colnames(proportions)) ||
      !setequal(colnames(proportions), required_cell_types)) {
    stop("The authoritative proportions must contain the canonical LM22 types", call. = FALSE)
  }
  if (!is.numeric(mean_threshold) || length(mean_threshold) != 1L ||
      !is.finite(mean_threshold) || mean_threshold < 0 ||
      !is.numeric(zero_floor) || length(zero_floor) != 1L ||
      !is.finite(zero_floor) || zero_floor <= 0) {
    stop("The threshold and zero floor must be valid", call. = FALSE)
  }

  combined <- vapply(
    group_map,
    function(cell_types) {
      rowSums(proportions[, cell_types, drop = FALSE])
    },
    numeric(nrow(proportions))
  )
  if (nrow(proportions) == 1L) {
    combined <- matrix(
      combined,
      nrow = 1L,
      dimnames = list(rownames(proportions), names(group_map))
    )
  } else {
    rownames(combined) <- rownames(proportions)
  }
  cohort_means <- colMeans(combined)
  retained <- cohort_means >= mean_threshold
  adjusted <- combined[, retained, drop = FALSE]
  adjusted[adjusted == 0] <- zero_floor
  tca_weights <- sweep(adjusted, 1L, rowSums(adjusted), "/")

  list(
    combined = combined,
    cohort_means = cohort_means,
    retained = retained,
    retained_groups = colnames(combined)[retained],
    zero_counts = colSums(combined == 0),
    tca_weights = tca_weights
  )
}

require_matrix_equal <- function(observed, expected, tolerance, label) {
  if (!is.matrix(observed) || !is.numeric(observed) ||
      !is.matrix(expected) || !is.numeric(expected)) {
    stop(sprintf("The %s values must be numeric matrices", label), call. = FALSE)
  }
  if (!identical(rownames(observed), rownames(expected)) ||
      !identical(colnames(observed), colnames(expected))) {
    stop(sprintf("The %s row or column order is incorrect", label), call. = FALSE)
  }
  if (!is.numeric(tolerance) || length(tolerance) != 1L ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("The matrix comparison tolerance must be non-negative", call. = FALSE)
  }
  maximum_difference <- max(abs(observed - expected))
  if (!is.finite(maximum_difference) || maximum_difference > tolerance) {
    stop(
      sprintf(
        "The %s values exceed tolerance %.3g; maximum difference is %.17g",
        label,
        tolerance,
        maximum_difference
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

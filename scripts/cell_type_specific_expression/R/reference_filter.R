reference_supported_cell_types <- unique(unname(reference_population_map))

summarize_bed_values <- function(values) {
  if (!is.matrix(values) || !is.numeric(values) || any(!is.finite(values))) {
    stop("BED expression values must be a finite numeric matrix", call. = FALSE)
  }
  negative <- values < 0
  negative_count <- rowSums(negative)
  logged <- suppressWarnings(log2(values + 1))
  valid_log <- rowSums(negative) == 0L
  tibble::tibble(
    negative_count = as.integer(negative_count),
    negative_percentage = 100 * negative_count / ncol(values),
    minimum_cpm = apply(values, 1L, min),
    mean_negative_cpm = vapply(seq_len(nrow(values)), function(i) {
      if (negative_count[[i]] == 0L) NA_real_ else mean(values[i, negative[i, ]])
    }, numeric(1)),
    n_samples = ncol(values), nonnegative = negative_count == 0L,
    mean_log2_cpm1 = ifelse(valid_log, rowMeans(logged), NA_real_),
    median_log2_cpm1 = vapply(seq_len(nrow(values)), function(i) {
      if (valid_log[[i]]) stats::median(logged[i, ]) else NA_real_
    }, numeric(1))
  )
}

apply_reference_rules <- function(bed_summary, cell_type, reference_summary,
                                  threshold = 0.01) {
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
    stop("Expression threshold must be one finite number", call. = FALSE)
  }
  result <- bed_summary
  result$reference_gene_matched <- NA
  result$reference_mean_log2_cpm1 <- NA_real_
  result$low_deconvolution_expression <- FALSE
  result$low_reference_expression <- FALSE
  result$comparison_status <- "reference_not_provided"
  result$retained <- result$nonnegative
  if (is.null(reference_summary)) return(result)
  if (!(cell_type %in% reference_supported_cell_types)) {
    result$comparison_status <- "no_reference_cell_type"
    return(result)
  }
  ref <- dplyr::filter(reference_summary, .data$cell_type == .env$cell_type)
  if (nrow(ref) == 0L) {
    result$comparison_status <- "reference_cell_type_unavailable"
    return(result)
  }
  bed_keys <- reference_gene_key(result$gene_id, unique = TRUE)
  ref_keys <- reference_gene_key(ref$gene_id, unique = TRUE)
  index <- match(bed_keys, ref_keys)
  matched <- !is.na(index)
  result$reference_gene_matched <- matched
  result$reference_mean_log2_cpm1[matched] <- ref$mean_log2_cpm1[index[matched]]
  result$comparison_status <- ifelse(matched, "compared", "reference_gene_unmatched")
  result$low_deconvolution_expression <- result$nonnegative &
    result$mean_log2_cpm1 <= threshold
  result$low_reference_expression <- matched &
    result$reference_mean_log2_cpm1 <= threshold
  result$retained <- result$nonnegative & matched &
    result$mean_log2_cpm1 > threshold & result$reference_mean_log2_cpm1 > threshold
  result
}

fit_reference_regression <- function(comparison) {
  needed <- c("reference_mean_log2_cpm1", "deconvolution_mean_log2_cpm1")
  if (nrow(comparison) < 3L) stop("Regression requires at least three shared genes", call. = FALSE)
  if (any(vapply(comparison[needed], function(x) !all(is.finite(x)), logical(1)))) {
    stop("Regression values must be finite", call. = FALSE)
  }
  if (any(vapply(comparison[needed], function(x) stats::var(x) == 0, logical(1)))) {
    stop("Regression requires variation in both profiles", call. = FALSE)
  }
  fit <- stats::lm(deconvolution_mean_log2_cpm1 ~ reference_mean_log2_cpm1,
                   data = comparison)
  fit_summary <- suppressWarnings(summary(fit))
  standardized <- if (!is.finite(fit_summary$sigma) || fit_summary$sigma <= sqrt(.Machine$double.eps)) {
    rep(NA_real_, nrow(comparison))
  } else suppressWarnings(stats::rstandard(fit))
  genes <- comparison |>
    dplyr::mutate(fitted_value = unname(stats::fitted(fit)),
                  residual = unname(stats::residuals(fit)),
                  standardized_residual = standardized)
  metrics <- tibble::tibble(
    n_genes = nrow(comparison),
    pearson_r = stats::cor(comparison[[needed[[1L]]]], comparison[[needed[[2L]]]], method = "pearson"),
    spearman_rho = stats::cor(comparison[[needed[[1L]]]], comparison[[needed[[2L]]]], method = "spearman"),
    r_squared = unname(fit_summary$r.squared), intercept = unname(stats::coef(fit)[[1L]]),
    slope = unname(stats::coef(fit)[[2L]]))
  list(genes = genes, metrics = metrics, fit = fit)
}

apply_residual_cutoff <- function(genes, cutoff) {
  if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff) || cutoff <= 0) {
    stop("Residual cutoff must be one positive finite number", call. = FALSE)
  }
  if (anyNA(genes$standardized_residual)) {
    stop("Standardized residuals are undefined; residual filtering cannot be applied", call. = FALSE)
  }
  dplyr::mutate(genes, residual_excluded = abs(.data$standardized_residual) > cutoff)
}

# First-update screening for the pinned TCA 1.2.1 least-squares estimator.
# Fit-stage only: do not change shared TCA or export implementations.
tca_variance_flags <- function(gene_id, estimates, lower_bound) {
  tibble::tibble(
    gene_id = gene_id, component = names(estimates), estimate = unname(estimates),
    lower_bound = lower_bound,
    reason = dplyr::case_when(
      !is.finite(estimates) ~ "nonfinite_variance",
      estimates < 0 ~ "negative_variance",
      estimates < lower_bound ~ "below_variance_bound",
      TRUE ~ NA_character_
    )
  ) |>
    dplyr::filter(!is.na(.data$reason)) |>
    dplyr::mutate(action = dplyr::if_else(.data$reason == "below_variance_bound", "warn", "exclude"))
}

screen_tca_variances <- function(X, W, C2 = NULL, log_file) {
  validate_tca_version()
  validate_tca_inputs(X, W, C2)
  cfg <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)
  n <- nrow(W)
  k <- ncol(W)
  if (is.null(C2)) C2 <- matrix(numeric(), n, 0L)
  W_norms <- rowSums(W^2)^0.5
  # Same first-update designs, bounds and column normalization as
  # TCA:::tca.fit_means_vars with vars.mle=FALSE and constrain_mu=TRUE.
  design <- cbind(W / t(pracma::repmat(W_norms, k, 1)),
                  if (ncol(C2)) C2 / t(pracma::repmat(W_norms, ncol(C2), 1)) else C2)
  variance_design <- cbind(W^2, rep(1, n))
  norms <- colSums(variance_design^2)^0.5
  variance_design <- variance_design / pracma::repmat(norms, n, 1)
  if (any(!is.finite(variance_design))) stop("Non-finite preflight variance design", call. = FALSE)
  append_tca_log(log_file, sprintf("stage=tca_preflight event=start genes=%d samples=%d groups=%d", nrow(X), n, k))
  solve_gene <- function(gene, phase, ...) {
    tryCatch(pracma::lsqlincon(...), error = function(error) {
      stop(sprintf("TCA preflight gene=%s phase=%s solver error: %s", gene, phase, conditionMessage(error)), call. = FALSE)
    })
  }
  mean_lower <- c(rep(min(X) + cfg$mu_epsilon, k), rep(-cfg$lsqlincon_inf, ncol(C2)))
  mean_upper <- c(rep(max(X) - cfg$mu_epsilon, k), rep(cfg$lsqlincon_inf, ncol(C2)))
  coefficients <- t(vapply(seq_len(nrow(X)), function(j) {
    gene <- rownames(X)[j]
    mean_coefficients <- solve_gene(gene, "mean", design, X[j, ] / W_norms,
      lb = mean_lower, ub = mean_upper)
    if (any(!is.finite(mean_coefficients))) {
      stop(sprintf("TCA preflight gene=%s phase=mean solver returned non-finite coefficients", gene), call. = FALSE)
    }
    if (j %% 1000L == 0L) message(sprintf("stage=tca_preflight event=progress phase=means genes_checked=%d total_genes=%d", j, nrow(X)))
    mean_coefficients
  }, numeric(ncol(design))))
  # Preserve upstream's full matrix products and separate C2 addition. A
  # per-gene dot product can change last-bit rounding at these extreme scales.
  U <- (tcrossprod(W, coefficients[, seq_len(k), drop = FALSE]) +
    tcrossprod(C2, coefficients[, k + seq_len(ncol(C2)), drop = FALSE]) - t(X))^2
  reports <- vector("list", nrow(X))
  for (j in seq_len(nrow(X))) {
    gene <- rownames(X)[j]
    if (any(!is.finite(U[, j]))) stop(sprintf("TCA preflight gene=%s phase=variance non-finite squared residuals", gene), call. = FALSE)
    estimates <- solve_gene(gene, "variance", variance_design, U[, j],
      lb = rep(cfg$min_sd, k + 1L) * norms) / norms
    names(estimates) <- c(colnames(W), "tau_squared_contribution")
    reports[[j]] <- tca_variance_flags(gene, estimates, cfg$min_sd)
    if (j %% 1000L == 0L) {
      message(sprintf("stage=tca_preflight event=progress phase=variances genes_checked=%d total_genes=%d", j, nrow(X)))
    }
  }
  report <- dplyr::bind_rows(reports)
  excluded <- report |>
    dplyr::filter(.data$action == "exclude") |>
    dplyr::distinct(.data$gene_id) |>
    dplyr::mutate(reason = "invalid_first_update_variance")
  # The model retains every flagged coefficient, including warning-only rows.
  purrr::pwalk(report, function(gene_id, component, estimate, lower_bound, reason, action) {
    append_tca_log(log_file, sprintf(
      "stage=tca_preflight event=flag gene_id=%s component=%s estimate=%.17g lower_bound=%.17g reason=%s action=%s",
      gene_id, component, estimate, lower_bound, reason, action))
  })
  keep <- !rownames(X) %in% excluded$gene_id
  append_tca_log(log_file, sprintf("stage=tca_preflight event=complete excluded_genes=%d warning_coefficients=%d retained_genes=%d",
    nrow(excluded), sum(report$action == "warn"), sum(keep)))
  if (!any(keep)) stop("No genes remain after TCA variance preflight", call. = FALSE)
  list(X = X[keep, , drop = FALSE], excluded = excluded, report = report)
}

attach_tca_preflight_exclusions <- function(result, original) {
  preflight <- original$tca_preflight
  if (is.null(preflight) || !is.null(original$gene_filter)) return(result)
  genes <- preflight$original_gene_ids
  removed <- preflight$excluded_gene_ids
  if (!identical(preflight$method, "first_iteration_variance") ||
      !is.character(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes) ||
      !is.character(removed) || anyNA(removed) || any(!nzchar(removed)) || anyDuplicated(removed) ||
      !all(removed %in% genes) ||
      !identical(genes[!genes %in% removed], rownames(original$mus_hat))) {
    stop("Invalid TCA preflight gene exclusion record", call. = FALSE)
  }
  filter <- result$model$gene_filter
  combined <- unique(c(removed, filter$excluded_gene_ids))
  filter$original_gene_ids <- genes
  filter$excluded_gene_ids <- genes[genes %in% combined]
  if (!identical(genes[!genes %in% combined], rownames(result$model$mus_hat))) {
    stop("TCA preflight and cleanup gene exclusions do not match retained genes", call. = FALSE)
  }
  prefit_report <- tibble::tibble(gene_id = removed, reason = rep("invalid_first_update_variance", length(removed)))
  result$report <- dplyr::bind_rows(prefit_report, result$report)
  filter$report <- result$report
  result$model$gene_filter <- filter
  result
}

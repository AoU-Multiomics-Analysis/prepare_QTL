hspe_sample_seeds <- function(sample_ids, random_seed) {
  if (length(random_seed) != 1L || !is.numeric(random_seed) ||
      !is.finite(random_seed) || random_seed < 1 ||
      random_seed > .Machine$integer.max || random_seed != trunc(random_seed)) {
    stop("random_seed must be one positive R integer", call. = FALSE)
  }
  purrr::map_int(sample_ids, function(id) {
    key <- paste0(as.integer(random_seed), ":", enc2utf8(id))
    hash <- digest::digest(key, algo = "sha256", serialize = FALSE)
    value <- as.double(strtoi(substr(hash, 1, 7), base = 16)) * 16 +
      strtoi(substr(hash, 8, 8), base = 16)
    as.integer(value %% 2147483646 + 1)
  })
}

prepare_hspe_batches <- function(inputs, batch_size = 100L,
                                 marker_fraction = .1, random_seed = 20260901L) {
  if (length(batch_size) != 1L || !is.numeric(batch_size) ||
      !is.finite(batch_size) || batch_size < 1 ||
      batch_size > .Machine$integer.max || batch_size != trunc(batch_size)) {
    stop("batch_size must be one positive integer", call. = FALSE)
  }
  if (length(marker_fraction) != 1L || !is.finite(marker_fraction) ||
      marker_fraction <= 0 || marker_fraction > 1) {
    stop("marker_fraction must be in (0, 1]", call. = FALSE)
  }
  version <- validate_hspe_version()
  validate_bulk_log(t(inputs$Y))
  assert_identical_ids(colnames(inputs$Y), colnames(inputs$references), "HSPE gene")
  cell_types <- rownames(inputs$references)
  pure <- stats::setNames(as.list(seq_along(cell_types)), cell_types)
  ranked <- hspe::find_markers(Y = inputs$references, pure_samples = pure,
                              marker_method = "ratio")$L
  if (any(lengths(ranked) == 0L)) {
    stop("Every cell type must have at least one selected marker", call. = FALSE)
  }
  # Match HSPE 0.1: fractions below one select a fraction; one selects one marker.
  counts <- if (marker_fraction < 1) {
    pmax(1L, as.integer(floor(marker_fraction * lengths(ranked))))
  } else {
    rep(1L, length(ranked))
  }
  selected <- purrr::map2(ranked, counts, ~ unname(.x[seq_len(.y)]))
  gene_indices <- unique(unlist(selected, use.names = FALSE))
  markers <- purrr::map(selected, ~ match(.x, gene_indices))
  marker_table <- purrr::imap_dfr(selected, function(indices, cell_type) {
    tibble::tibble(cell_type = cell_type, marker_rank = seq_along(indices),
                   gene_symbol = colnames(inputs$references)[indices])
  })
  ids <- rownames(inputs$Y)
  samples <- tibble::tibble(sample_id = ids,
                            random_seed = hspe_sample_seeds(ids, random_seed))
  indices <- split(seq_along(ids), ceiling(seq_along(ids) / batch_size))
  batches <- purrr::map(indices, function(rows) {
    list(Y = inputs$Y[rows, gene_indices, drop = FALSE])
  })
  metadata <- list(
    hspe_version = version, optimizer = "DEoptimR", random_seed = random_seed,
    seed_strategy = "sha256(base_seed:UTF8_sample_id), first 8 hex digits modulo 2147483646 plus 1",
    marker_method = "ratio", marker_fraction = marker_fraction,
    marker_counts = as.list(stats::setNames(counts, cell_types)),
    overlap_count = inputs$overlap_count, overlap_fraction = inputs$overlap_fraction,
    quantile_normalize = inputs$quantile_normalize, sample_count = length(ids),
    batch_size = batch_size, batch_count = length(batches),
    selected_marker_count = length(gene_indices),
    reference_dimensions = list(cell_types = length(cell_types), genes = ncol(inputs$references)),
    lm22_qc = inputs$lm22_qc
  )
  shared <- list(references = inputs$references[, gene_indices, drop = FALSE],
                 markers = markers, marker_table = marker_table, samples = samples,
                 metadata = metadata)
  list(shared = shared, batches = unname(batches))
}

fit_hspe_batch <- function(batch, shared) {
  version <- validate_hspe_version()
  validate_bulk_log(t(batch$Y))
  assert_identical_ids(colnames(shared$references), colnames(batch$Y), "HSPE marker")
  ids <- rownames(batch$Y)
  indices <- match(ids, shared$samples$sample_id)
  if (anyNA(indices)) {
    stop("Batch contains an unexpected sample", call. = FALSE)
  }
  fits <- purrr::map(seq_along(ids), function(i) {
    seed <- shared$samples$random_seed[indices[i]]
    message(sprintf("stage=hspe_batch sample_index=%d sample_count=%d random_seed=%d start_time=%s",
                    i, length(ids), seed, format(Sys.time(), tz = "UTC", usetz = TRUE)))
    fit <- hspe::hspe(Y = batch$Y[i, , drop = FALSE], references = shared$references,
                      markers = shared$markers, n_markers = lengths(shared$markers),
                      seed = seed, sto = TRUE)
    values <- as.numeric(fit$estimates)
    if (length(values) != nrow(shared$references) || any(!is.finite(values)) ||
        any(values < 0) || abs(sum(values) - 1) > 1e-8) {
      stop("HSPE returned invalid proportions", call. = FALSE)
    }
    optimizer <- fit$diag$opt[1, "opt"][[1]]$opt_out
    if (is.null(optimizer) || length(optimizer$convergence) != 1L) {
      stop("HSPE did not return optimizer diagnostics", call. = FALSE)
    }
    diagnostics <- tibble::tibble(sample_id = ids[i], random_seed = seed,
      iterations = optimizer$iter, convergence = optimizer$convergence,
      loss = optimizer$value)
    message(sprintf("stage=hspe_batch sample_index=%d iterations=%d convergence=%d completion_time=%s",
                    i, optimizer$iter, optimizer$convergence,
                    format(Sys.time(), tz = "UTC", usetz = TRUE)))
    list(values = values, diagnostics = diagnostics)
  })
  proportions <- do.call(rbind, purrr::map(fits, "values"))
  dimnames(proportions) <- list(ids, rownames(shared$references))
  list(proportions = proportions, diagnostics = purrr::map_dfr(fits, "diagnostics"),
       hspe_version = version, optimizer_version = as.character(utils::packageVersion("DEoptimR")))
}

merge_hspe_batches <- function(results, shared) {
  if (length(results) == 0L) stop("No batch results supplied", call. = FALSE)
  purrr::walk(results, function(result) {
    p <- result$proportions
    assert_identical_ids(rownames(shared$references), colnames(p), "HSPE cell-type column")
    if (!is.matrix(p) || !is.numeric(p) || any(!is.finite(p))) {
      stop("Batch proportions must be finite numeric matrices", call. = FALSE)
    }
    if (any(p < 0) || any(abs(rowSums(p) - 1) > 1e-8)) {
      stop("Batch proportions must be nonnegative and sum to one", call. = FALSE)
    }
    assert_identical_ids(rownames(p), result$diagnostics$sample_id, "HSPE diagnostics sample")
    if (!identical(result$hspe_version, shared$metadata$hspe_version)) {
      stop("HSPE versions differ between preparation and fitting", call. = FALSE)
    }
  })
  proportions <- do.call(rbind, purrr::map(results, "proportions"))
  ids <- rownames(proportions)
  if (anyDuplicated(ids)) stop("Batch results contain duplicate samples", call. = FALSE)
  if (!setequal(ids, shared$samples$sample_id)) {
    stop("Batch results have missing or unexpected samples", call. = FALSE)
  }
  diagnostics <- purrr::map_dfr(results, "diagnostics")
  order <- match(shared$samples$sample_id, ids)
  proportions <- proportions[order, , drop = FALSE]
  diagnostics <- diagnostics[order, ]
  if (!identical(as.integer(diagnostics$random_seed), shared$samples$random_seed)) {
    stop("Batch sample seeds differ from preparation", call. = FALSE)
  }
  versions <- unique(purrr::map_chr(results, "optimizer_version"))
  if (length(versions) != 1L) stop("Batch optimizer versions differ", call. = FALSE)
  metadata <- shared$metadata
  metadata$batch_count <- length(results)
  metadata$optimizer_version <- versions[[1]]
  metadata$nonconverged_samples <- sum(diagnostics$convergence != 0)
  list(proportions = proportions, diagnostics = diagnostics, metadata = metadata)
}

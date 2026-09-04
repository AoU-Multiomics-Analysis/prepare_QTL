read_bed_header <- function(path) {
  names(readr::read_tsv(path, n_max = 0, name_repair = "minimal", show_col_types = FALSE,
                        col_types = readr::cols(.default = readr::col_character())))
}

scan_bed_summary <- function(path, cell_type, expected_genes, expected_samples, chunk_size) {
  header <- read_bed_header(path)
  if (length(header) <= 4L || !identical(header[1:4], expression_bed_columns()) ||
      anyNA(header) || any(!nzchar(header)) || anyDuplicated(header) > 0L) {
    stop(sprintf("Invalid BED header for '%s'", cell_type), call. = FALSE)
  }
  samples <- header[-(1:4)]
  if (length(samples) != expected_samples) stop("BED sample count does not match inventory", call. = FALSE)
  seen <- character()
  pieces <- list()
  callback <- readr::SideEffectChunkCallback$new(function(chunk, position) {
    if (nrow(readr::problems(chunk)) > 0L) stop("BED parse error", call. = FALSE)
    coordinates <- dplyr::select(chunk, dplyr::all_of(expression_bed_columns()))
    if (anyNA(coordinates) || any(!nzchar(coordinates$gene_id)) ||
        any(!nzchar(coordinates[["#chr"]])) || any(coordinates$start < 0L) ||
        any(coordinates$end <= coordinates$start)) stop("BED coordinates and IDs must be valid", call. = FALSE)
    if (anyDuplicated(coordinates$gene_id) || any(coordinates$gene_id %in% seen)) {
      stop("BED gene IDs must be unique", call. = FALSE)
    }
    seen <<- c(seen, coordinates$gene_id)
    values <- as.matrix(dplyr::select(chunk, dplyr::all_of(samples)))
    statistics <- summarize_bed_values(values)
    pieces[[length(pieces) + 1L]] <<- dplyr::bind_cols(coordinates, statistics)
  })
  readr::read_tsv_chunked(path, callback = callback, chunk_size = chunk_size,
    col_types = readr::cols(`#chr` = readr::col_character(), start = readr::col_integer(),
      end = readr::col_integer(), gene_id = readr::col_character(), .default = readr::col_double()),
    progress = FALSE)
  result <- dplyr::bind_rows(pieces)
  if (nrow(result) != expected_genes) stop("BED gene count does not match inventory", call. = FALSE)
  list(summary = result, samples = samples)
}

write_filtered_bed <- function(input_path, output_path, retained_ids, samples, chunk_size) {
  temporary_path <- tempfile("filtered-bed-", tmpdir = dirname(output_path), fileext = ".bed")
  on.exit(unlink(temporary_path), add = TRUE)
  state <- new.env(parent = emptyenv())
  state$wrote <- FALSE
  callback <- readr::SideEffectChunkCallback$new(function(chunk, position) {
    kept <- dplyr::filter(chunk, .data$gene_id %in% retained_ids)
    if (nrow(kept) > 0L) {
      chunk_connection <- file(temporary_path, if (state$wrote) "ab" else "wb")
      if (!state$wrote) writeLines(paste(names(kept), collapse = "\t"), chunk_connection)
      lines <- apply(kept, 1L, function(row) paste(row, collapse = "\t"))
      writeLines(lines, chunk_connection)
      close(chunk_connection)
      state$wrote <- TRUE
    }
  })
  readr::read_tsv_chunked(input_path, callback = callback, chunk_size = chunk_size,
    col_types = readr::cols(`#chr` = readr::col_character(), start = readr::col_integer(),
      end = readr::col_integer(), gene_id = readr::col_character(), .default = readr::col_double()),
    progress = FALSE)
  if (!state$wrote) stop("Filtering left an empty BED", call. = FALSE)
  input_connection <- file(temporary_path, "rt")
  output_connection <- gzfile(output_path, "wt")
  repeat {
    lines <- readLines(input_connection, n = 1024L, warn = FALSE)
    if (length(lines) == 0L) break
    writeLines(lines, output_connection)
  }
  close(input_connection)
  close(output_connection)
  invisible(output_path)
}

filter_cell_type_beds <- function(inventory, bed_paths, output_dir, reference_summary = NULL,
                                  min_mean_log2_cpm1 = 0.01, residual_cutoff = NULL,
                                  chunk_size = 256L) {
  validate_scatter_inventory(inventory)
  chunk_size <- validate_tensor_positive_integer(chunk_size, "chunk_size")
  bed_paths <- validate_scatter_bed_paths(bed_paths, inventory$path)
  if (!all(file.exists(bed_paths))) stop("Every cell-type BED must exist", call. = FALSE)
  actual_hashes <- purrr::map_chr(bed_paths, ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE))
  if (!identical(actual_hashes, inventory$sha256)) stop("BED checksum does not match inventory", call. = FALSE)
  if (!is.null(residual_cutoff) && is.null(reference_summary)) {
    stop("A residual cutoff requires a reference", call. = FALSE)
  }
  if (!is.null(reference_summary)) {
    required_reference <- c("gene_id", "cell_type", "n_samples", "mean_log2_cpm1", "median_log2_cpm1")
    if (!inherits(reference_summary, "data.frame") ||
        !all(required_reference %in% names(reference_summary)) ||
        anyNA(reference_summary$gene_id) || anyNA(reference_summary$cell_type) ||
        any(!is.finite(reference_summary$mean_log2_cpm1))) {
      stop("Reference summary is invalid or incomplete", call. = FALSE)
    }
    duplicate_pairs <- paste(reference_summary$gene_id, reference_summary$cell_type, sep = "\r")
    if (anyDuplicated(duplicate_pairs) > 0L) stop("Reference gene and cell-type pairs must be unique", call. = FALSE)
  }
  dir.create(file.path(output_dir, "beds"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  comparisons <- list(); metrics <- list(); output_paths <- character(nrow(inventory))
  cohort_samples <- NULL
  for (i in seq_len(nrow(inventory))) {
    cell_type <- inventory$cell_group[[i]]; slug <- inventory$slug[[i]]
    message(sprintf("stage=filter_cell_type_beds cell_type=%s start_time=%s", cell_type, tensor_utc_time()))
    scanned <- scan_bed_summary(bed_paths[[i]], cell_type, inventory$n_genes[[i]],
                                inventory$n_samples[[i]], chunk_size)
    if (is.null(cohort_samples)) cohort_samples <- scanned$samples
    assert_identical_ids(cohort_samples, scanned$samples, "BED sample")
    decision <- apply_reference_rules(scanned$summary, cell_type, reference_summary,
                                      min_mean_log2_cpm1)
    decision$cell_type <- cell_type
    decision$residual_excluded <- FALSE
    decision$fitted_value <- NA_real_; decision$residual <- NA_real_
    decision$standardized_residual <- NA_real_
    fit_rows <- which(decision$retained & decision$comparison_status == "compared")
    status <- if (is.null(reference_summary)) "reference_not_provided" else
      if (!(cell_type %in% reference_supported_cell_types)) "no_reference_cell_type" else "available"
    baseline <- tibble::tibble(n_genes = length(fit_rows), pearson_r = NA_real_, spearman_rho = NA_real_,
      r_squared = NA_real_, intercept = NA_real_, slope = NA_real_)
    if (length(fit_rows) >= 3L && stats::var(decision$reference_mean_log2_cpm1[fit_rows]) > 0 &&
        stats::var(decision$mean_log2_cpm1[fit_rows]) > 0) {
      fit_input <- tibble::tibble(gene_id = decision$gene_id[fit_rows],
        reference_mean_log2_cpm1 = decision$reference_mean_log2_cpm1[fit_rows],
        deconvolution_mean_log2_cpm1 = decision$mean_log2_cpm1[fit_rows])
      fitted <- fit_reference_regression(fit_input)
      baseline <- fitted$metrics
      decision$fitted_value[fit_rows] <- fitted$genes$fitted_value
      decision$residual[fit_rows] <- fitted$genes$residual
      decision$standardized_residual[fit_rows] <- fitted$genes$standardized_residual
      if (!is.null(residual_cutoff)) {
        residual_result <- apply_residual_cutoff(fitted$genes, residual_cutoff)
        excluded <- residual_result$gene_id[residual_result$residual_excluded]
        decision$residual_excluded <- decision$gene_id %in% excluded
        decision$retained <- decision$retained & !decision$residual_excluded
      }
    } else if (!is.null(residual_cutoff) && status == "available") {
      stop(sprintf("Cell type '%s' has unavailable regression residuals", cell_type), call. = FALSE)
    } else if (status == "available") status <- "insufficient_genes_or_variation"
    retained_ids <- decision$gene_id[decision$retained]
    output_paths[[i]] <- file.path(output_dir, "beds", sprintf("%s.filtered.bed.gz", slug))
    write_filtered_bed(bed_paths[[i]], output_paths[[i]], retained_ids, scanned$samples, chunk_size)
    comparisons[[i]] <- decision
    metric_rows <- dplyr::mutate(baseline, cell_type = cell_type, slug = slug,
      comparison_status = status, n_original = nrow(decision),
      n_negative_excluded = sum(!decision$nonnegative),
      n_reference_excluded = sum(decision$nonnegative & !decision$retained & !decision$residual_excluded),
      n_residual_excluded = sum(decision$residual_excluded), n_retained = sum(decision$retained),
      metric_set = "baseline", .before = 1)
    if (!is.null(residual_cutoff) && status == "available") {
      retained_rows <- which(decision$retained)
      retained_input <- tibble::tibble(gene_id = decision$gene_id[retained_rows],
        reference_mean_log2_cpm1 = decision$reference_mean_log2_cpm1[retained_rows],
        deconvolution_mean_log2_cpm1 = decision$mean_log2_cpm1[retained_rows])
      retained_metrics <- if (nrow(retained_input) >= 3L &&
          stats::var(retained_input$reference_mean_log2_cpm1) > 0 &&
          stats::var(retained_input$deconvolution_mean_log2_cpm1) > 0) {
        fit_reference_regression(retained_input)$metrics
      } else tibble::tibble(n_genes = nrow(retained_input), pearson_r = NA_real_, spearman_rho = NA_real_,
        r_squared = NA_real_, intercept = NA_real_, slope = NA_real_)
      metric_rows <- dplyr::bind_rows(metric_rows,
        dplyr::mutate(retained_metrics, cell_type = cell_type, slug = slug,
          comparison_status = if (nrow(retained_input) >= 3L) "available" else "insufficient_retained_genes",
          n_original = nrow(decision), n_negative_excluded = sum(!decision$nonnegative),
          n_reference_excluded = sum(decision$nonnegative & !decision$retained & !decision$residual_excluded),
          n_residual_excluded = sum(decision$residual_excluded), n_retained = sum(decision$retained),
          metric_set = "retained", .before = 1))
    }
    metrics[[i]] <- metric_rows
    message(sprintf("stage=filter_cell_type_beds cell_type=%s retained=%d completion_time=%s",
                    cell_type, length(retained_ids), tensor_utc_time()))
  }
  all_comparisons <- dplyr::bind_rows(comparisons)
  negative <- dplyr::select(all_comparisons, "cell_type", "gene_id",
    "negative_count", "negative_percentage", "minimum_cpm",
    "mean_negative_cpm", "n_samples")
  readr::write_tsv(negative, file.path(output_dir, "negative_summary.tsv.gz"), na = "NA")
  readr::write_tsv(all_comparisons, file.path(output_dir, "gene_comparison.tsv.gz"), na = "NA")
  readr::write_tsv(dplyr::bind_rows(metrics), file.path(output_dir, "filter_metrics.tsv"), na = "NA")
  save_negative_plots(negative, file.path(output_dir, "plots"))
  purrr::walk2(seq_len(nrow(inventory)), inventory$slug, function(i, slug) {
    save_reference_plots(comparisons[[i]], inventory$cell_group[[i]], slug, file.path(output_dir, "plots"))
    if (!is.null(residual_cutoff)) {
      retained_plot <- dplyr::filter(comparisons[[i]], .data$retained)
      if (nrow(retained_plot) > 0L) {
        save_reference_plots(retained_plot, inventory$cell_group[[i]], slug,
                             file.path(output_dir, "plots"), post = TRUE)
      }
    }
  })
  output_basenames <- basename(output_paths)
  filtered_inventory <- inventory |>
    dplyr::mutate(path = output_basenames,
      sha256 = purrr::map_chr(output_paths, ~ digest::digest(file = .x, algo = "sha256", serialize = FALSE)),
      n_genes = purrr::map_int(comparisons, ~ sum(.x$retained)))
  readr::write_tsv(filtered_inventory, file.path(output_dir, "filtered_inventory.tsv"))
  writeLines(file.path("beds", output_basenames), file.path(output_dir, "filtered_beds.txt"))
  invisible(list(inventory = filtered_inventory, comparisons = all_comparisons, metrics = dplyr::bind_rows(metrics)))
}

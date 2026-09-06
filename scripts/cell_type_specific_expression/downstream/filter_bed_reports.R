require_filter_files <- function(paths) {
  if (length(paths) == 0L || anyNA(paths) || any(!nzchar(paths))) {
    stop("Required local input paths are empty", call. = FALSE)
  }
  if (any(grepl("^[a-zA-Z][a-zA-Z0-9+.-]*://", paths))) {
    stop("Localization error: required input is still a cloud URI", call. = FALSE)
  }
  if (any(file.access(paths, 4L) != 0L)) stop("Required local input is not readable", call. = FALSE)
  invisible(paths)
}

select_bed_inventory <- function(inventory, bed) {
  validate_scatter_inventory(inventory)
  selected <- dplyr::filter(inventory, .data$path == basename(bed))
  if (length(bed) != 1L || nrow(selected) != 1L) {
    stop("BED must match exactly one inventory row", call. = FALSE)
  }
  selected
}

merge_filter_reports <- function(original, inventories, comparisons, metrics, samples,
                                 output_dir, post_residual = FALSE) {
  validate_scatter_inventory(original)
  groups <- list(inventories, comparisons, metrics, samples)
  if (any(lengths(groups) != nrow(original))) {
    stop("Report shard count must match the original inventory", call. = FALSE)
  }
  require_filter_files(unlist(groups, use.names = FALSE))
  read_table <- function(path) readr::read_tsv(path, show_col_types = FALSE, progress = FALSE,
    name_repair = "minimal", col_types = readr::cols(.default = readr::col_guess(),
      gene_id = readr::col_character(), cell_type = readr::col_character()))
  # Inventory has no gene_id/cell_type fields; use its existing schema.
  parts <- purrr::map(inventories, ~ readr::read_tsv(.x, show_col_types = FALSE, progress = FALSE))
  if (any(purrr::map_int(parts, nrow) != 1L)) stop("Each shard must contain one cell type", call. = FALSE)
  combined <- dplyr::bind_rows(parts)
  validate_scatter_inventory(combined)
  if (!setequal(combined$slug, original$slug)) stop("Filtered cell types do not match inventory", call. = FALSE)
  order <- match(original$slug, combined$slug)
  combined <- combined[order, ]
  if (!identical(combined$cell_group, original$cell_group) ||
      any(combined$n_samples != original$n_samples)) stop("Shard labels or sample counts changed", call. = FALSE)
  sample_ids <- purrr::map(samples[order], ~ readLines(.x, warn = FALSE))
  purrr::walk2(sample_ids, combined$n_samples, function(ids, n) {
    if (length(ids) != n || any(!nzchar(ids)) || anyDuplicated(ids)) {
      stop("Invalid shard sample IDs", call. = FALSE)
    }
    assert_identical_ids(sample_ids[[1L]], ids, "BED sample")
  })
  comparison_parts <- purrr::map(comparisons[order], read_table)
  metric_parts <- purrr::map(metrics[order], ~ readr::read_tsv(.x, show_col_types = FALSE, progress = FALSE))
  purrr::walk(seq_len(nrow(combined)), function(i) {
    x <- comparison_parts[[i]]
    m <- metric_parts[[i]]
    if (nrow(x) != original$n_genes[[i]] || anyDuplicated(x$gene_id) ||
        anyNA(x$cell_type) || !all(x$cell_type == combined$cell_group[[i]]) ||
        nrow(m) < 1L || !all(m$cell_type == combined$cell_group[[i]]) ||
        sum(x$retained) != combined$n_genes[[i]]) {
      stop("Shard reports do not match their inventory", call. = FALSE)
    }
  })
  comparison <- dplyr::bind_rows(comparison_parts)
  negative <- dplyr::select(comparison, "cell_type", "gene_id", "negative_count",
    "negative_percentage", "minimum_cpm", "mean_negative_cpm", "n_samples")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(combined, file.path(output_dir, "filtered_inventory.tsv"))
  readr::write_tsv(comparison, file.path(output_dir, "gene_comparison.tsv.gz"), na = "NA")
  readr::write_tsv(negative, file.path(output_dir, "negative_summary.tsv.gz"), na = "NA")
  readr::write_tsv(dplyr::bind_rows(metric_parts), file.path(output_dir, "filter_metrics.tsv"), na = "NA")
  plot_dir <- file.path(output_dir, "plots")
  save_negative_plots(negative, plot_dir)
  purrr::walk(seq_len(nrow(combined)), function(i) {
    save_reference_plots(comparison_parts[[i]], combined$cell_group[[i]], combined$slug[[i]], plot_dir)
    if (post_residual) {
      save_reference_plots(dplyr::filter(comparison_parts[[i]], .data$retained),
        combined$cell_group[[i]], combined$slug[[i]], plot_dir, post = TRUE)
    }
  })
  invisible(combined)
}

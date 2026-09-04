save_negative_plots <- function(negative_summary, plot_dir) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  overview <- negative_summary |>
    dplyr::group_by(.data$cell_type) |>
    dplyr::summarise(percent_genes_with_negatives = 100 * mean(.data$negative_count > 0), .groups = "drop") |>
    ggplot2::ggplot(ggplot2::aes(x = .data$cell_type, y = .data$percent_genes_with_negatives)) +
    ggplot2::geom_col(fill = "#386cb0") +
    ggplot2::labs(x = NULL, y = "Genes with negative estimates (%)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(file.path(plot_dir, "negative_overview.pdf"), overview,
                  width = 7, height = 4.5, units = "in")
  ranked <- negative_summary |>
    dplyr::group_by(.data$gene_id) |>
    dplyr::summarise(rank_fraction = max(.data$negative_percentage), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$rank_fraction), .data$gene_id) |>
    dplyr::slice_head(n = 100L)
  heat_data <- dplyr::semi_join(negative_summary, ranked, by = "gene_id") |>
    dplyr::mutate(gene_id = factor(.data$gene_id, levels = rev(ranked$gene_id)))
  heat <- ggplot2::ggplot(heat_data, ggplot2::aes(x = .data$cell_type, y = .data$gene_id,
                                                  fill = .data$negative_percentage)) +
    ggplot2::geom_tile() + ggplot2::scale_fill_viridis_c(name = "Negative (%)") +
    ggplot2::labs(x = NULL, y = "Gene ID") + ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(plot_dir, "negative_fraction_heatmap.pdf"), heat,
                  width = 7, height = max(4, min(14, nrow(ranked) * 0.13)), units = "in", limitsize = FALSE)
}

save_reference_plots <- function(comparison, cell_type, slug, plot_dir, post = FALSE) {
  eligible <- dplyr::filter(comparison, .data$comparison_status == "compared",
                            .data$nonnegative, !.data$low_deconvolution_expression,
                            !.data$low_reference_expression,
                            is.finite(.data$reference_mean_log2_cpm1))
  if (nrow(eligible) < 1L) return(invisible(character()))
  suffix <- if (post) "post_residual" else "baseline"
  scatter <- ggplot2::ggplot(eligible, ggplot2::aes(
    x = .data$reference_mean_log2_cpm1, y = .data$mean_log2_cpm1)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = 2) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$residual_excluded), alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black") +
    ggplot2::labs(x = "Reference mean log2(CPM + 1)",
                  y = "Deconvolution mean log2(CPM + 1)", colour = "Residual excluded") +
    ggplot2::theme_minimal()
  scatter_path <- file.path(plot_dir, sprintf("%s.%s.scatter.pdf", slug, suffix))
  ggplot2::ggsave(scatter_path, scatter, width = 5.5, height = 5, units = "in")
  paths <- scatter_path
  if ("standardized_residual" %in% names(eligible) && any(is.finite(eligible$standardized_residual))) {
    residual <- ggplot2::ggplot(eligible, ggplot2::aes(
      x = .data$reference_mean_log2_cpm1, y = .data$standardized_residual)) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60") + ggplot2::geom_point(alpha = 0.7) +
      ggplot2::labs(x = "Reference mean log2(CPM + 1)", y = "Internally standardized residual") +
      ggplot2::theme_minimal()
    residual_path <- file.path(plot_dir, sprintf("%s.%s.residual.pdf", slug, suffix))
    ggplot2::ggsave(residual_path, residual, width = 5.5, height = 5, units = "in")
    paths <- c(paths, residual_path)
  }
  invisible(paths)
}

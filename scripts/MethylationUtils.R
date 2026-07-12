# Shared helpers for pb-CpG methylation filtering and cohort merging.

suppressPackageStartupMessages(library(data.table))

load_methylation_data <- function(file_path, filter_chroms = "X|Y|M|_", fence_k = 3, decomp_threads = 1L) {
    # pb-CpG-tools BED columns are `#chrom begin end mod_score type cov
    # est_mod_count est_unmod_count discretized_mod_score`; only the first six are
    # used downstream, so select them and drop the model-mode extras. The file is
    # bgzipped (BGZF), so decompress its blocks in parallel with bgzip and let fread
    # parse the stream. bgzip -d also handles plain gzip (single-threaded), so a
    # non-BGZF input still reads correctly. The leading `##key=value` metadata lines
    # mean skip = "#chrom" must remain so fread locates the real column header.
    keep_cols <- c("#chrom", "begin", "end", "mod_score", "type", "cov")
    loaded_data <- fread(
        cmd = sprintf("bgzip -c -d -@ %d %s", decomp_threads, shQuote(file_path)),
        skip = "#chrom",
        select = keep_cols
    )
    n_input_rows <- nrow(loaded_data)
    required_columns <- c("#chrom", "begin", "end", "mod_score", "type", "cov")
    missing_columns <- setdiff(required_columns, names(loaded_data))
    if (length(missing_columns) > 0) {
        stop("Missing required column(s) in ", file_path, ": ",
             paste(missing_columns, collapse = ", "),
             ". Expected pb-CpG-tools columns #chrom, begin, end, mod_score, type, and cov.")
    }

    loaded_data[, cov := as.numeric(cov)]
    if (all(is.na(loaded_data$cov))) stop("Column 'cov' is not numeric in ", file_path)
    if (!is.null(filter_chroms) && nzchar(filter_chroms)) {
        loaded_data <- loaded_data[!grepl(filter_chroms, `#chrom`)]
    }
    if (nrow(loaded_data) == 0) stop("No rows remain after chromosome filtering in ", file_path)

    call_types <- unique(as.character(loaded_data$type[!is.na(loaded_data$type)]))
    if (length(call_types) != 1) {
        stop("Expected one pb-CpG-tools 'type' per input file, but found: ",
             paste(call_types, collapse = ", "), ". Use one .combined.bed.gz file per sample.")
    }
    median_cov <- median(loaded_data$cov, na.rm = TRUE)
    message("Median cov for ", basename(file_path), ": ", round(median_cov, 3))
    if (!is.finite(median_cov) || median_cov <= 0) {
        stop("Median coverage must be positive after chromosome filtering in ", file_path)
    }

    logc <- log10(loaded_data$cov[!is.na(loaded_data$cov) & loaded_data$cov > 0])
    if (length(logc) == 0) {
        extreme_cut <- Inf
        message("No positive coverage values in ", basename(file_path), "; setting extreme-coverage cutoff to Inf")
    } else {
        qs <- quantile(logc, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
        extreme_cut <- 10^(qs[2] + fence_k * (qs[2] - qs[1]))
        message("Extreme-coverage cutoff for ", basename(file_path), ": ",
                round(extreme_cut), "x (implied CN ~", round(2 * extreme_cut / median_cov), ")")
    }

    loaded_data[, `:=`(
        implied_cn = 2 * cov / median_cov,
        extreme_cov_flag = fifelse(!is.na(cov) & cov >= extreme_cut, "extreme", "ok"),
        site_key = paste(`#chrom`, begin, end, sep = "*")
    )]
    setattr(loaded_data, "median_cov", median_cov)
    setattr(loaded_data, "extreme_cut", extreme_cut)
    setattr(loaded_data, "n_input_rows", n_input_rows)
    setattr(loaded_data, "call_type", call_types[[1]])
    loaded_data
}

read_manifest <- function(manifest_path) {
    manifest <- fread(manifest_path)
    missing_columns <- setdiff(c("sample_id", "file_path"), names(manifest))
    if (length(missing_columns) > 0) {
        stop("Input manifest must contain columns 'sample_id' and 'file_path'. Missing: ",
             paste(missing_columns, collapse = ", "))
    }
    manifest <- manifest[, .(sample_id = as.character(sample_id), file_path = as.character(file_path))]
    if (anyNA(manifest$sample_id) || any(!nzchar(manifest$sample_id))) stop("Input manifest contains an empty sample_id")
    if (anyDuplicated(manifest$sample_id)) stop("Each sample_id must occur exactly once in the input manifest")
    if (anyNA(manifest$file_path) || any(!nzchar(manifest$file_path))) stop("Input manifest contains an empty file_path")

    manifest_dir <- dirname(normalizePath(manifest_path))
    is_absolute <- grepl("^(/|~)", manifest$file_path)
    manifest[!is_absolute, file_path := file.path(manifest_dir, file_path)]
    manifest[, file_path := path.expand(file_path)]
    missing_files <- manifest[!file.exists(file_path), file_path]
    if (length(missing_files) > 0) stop("Input BED file(s) do not exist: ", paste(missing_files, collapse = ", "))
    manifest
}

read_file_list <- function(list_path, label) {
    if (!file.exists(list_path)) stop(label, " file list does not exist: ", list_path)
    paths <- scan(list_path, what = character(), quiet = TRUE)
    if (length(paths) == 0) stop(label, " file list is empty: ", list_path)
    list_dir <- dirname(normalizePath(list_path))
    is_absolute <- grepl("^(/|~)", paths)
    paths[!is_absolute] <- file.path(list_dir, paths[!is_absolute])
    paths <- path.expand(paths)
    missing_paths <- paths[!file.exists(paths)]
    if (length(missing_paths) > 0) stop(label, " file(s) do not exist: ", paste(missing_paths, collapse = ", "))
    paths
}

read_call_tables <- function(call_paths, label, additional_required_columns = character()) {
    required_columns <- c("sample_id", "#chrom", "begin", "end", "cov", "site_key", additional_required_columns)
    call_tables <- vector("list", length(call_paths))
    reference_columns <- NULL
    for (i in seq_along(call_paths)) {
        call_table <- fread(call_paths[i])
        missing_columns <- setdiff(required_columns, names(call_table))
        if (length(missing_columns) > 0) {
            stop(label, " file ", call_paths[i], " is missing: ", paste(missing_columns, collapse = ", "))
        }
        call_table[, `:=`(sample_id = as.character(sample_id), cov = as.numeric(cov))]
        if (is.null(reference_columns)) {
            reference_columns <- names(call_table)
        } else if (!identical(reference_columns, names(call_table))) {
            stop(label, " files do not use the same schema")
        }
        call_tables[[i]] <- call_table
    }
    all_calls <- rbindlist(call_tables, use.names = TRUE)
    duplicated_calls <- all_calls[, .N, by = .(sample_id, site_key)][N > 1]
    if (nrow(duplicated_calls) > 0) {
        stop("Found duplicated sample/site calls across ", label, " shards. Each sample must be present in exactly one shard and each site must occur once per sample.")
    }
    list(all_calls = all_calls, reference_columns = reference_columns)
}

read_sample_qc <- function(qc_paths, total_samples) {
    sample_qc <- rbindlist(lapply(qc_paths, fread), use.names = TRUE)
    if (!("sample_id" %in% names(sample_qc))) stop("Filtered sample-QC files must contain a sample_id column")
    if (anyDuplicated(sample_qc$sample_id)) stop("A sample appears in more than one filtered sample-QC file")
    if (nrow(sample_qc) != total_samples) {
        stop("Filtered sample-QC files contain ", nrow(sample_qc), " samples, but --TotalSamples is ",
             total_samples, ". These must match so the fraction threshold has the correct denominator.")
    }
    sample_qc
}

safe_mean <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else mean(x) }
safe_sd <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2) NA_real_ else sd(x) }
safe_cv <- function(x) { average <- safe_mean(x); deviation <- safe_sd(x); if (is.na(average) || is.na(deviation) || average == 0) NA_real_ else deviation / average }
safe_median <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else median(x) }
safe_min <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else min(x) }
safe_max <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else max(x) }
safe_mad <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else median(abs(x - median(x))) }
safe_spearman <- function(x, y) {
    valid <- is.finite(x) & is.finite(y)
    if (sum(valid) < 3) return(NA_real_)
    x <- x[valid]
    y <- y[valid]
    if (uniqueN(x) < 2 || uniqueN(y) < 2) return(NA_real_)
    suppressWarnings(cor(x, y, method = "spearman"))
}

inverse_normal_transform <- function(x) {
    transformed <- rep(NA_real_, length(x))
    valid <- is.finite(x)
    if (any(valid)) {
        ranks <- rank(x[valid], ties.method = "average")
        transformed[valid] <- qnorm((ranks - 0.5) / sum(valid))
    }
    transformed
}

write_filter_plots <- function(site_metadata, output_prefix) {
    failure_levels <- c("Insufficient minimum coverage", "Extreme coverage exclusion", "Low methylation MAD", "Pass all cohort filters")
    filter_summary <- merge(
        data.table(failure_reason = failure_levels), site_metadata[, .N, by = failure_reason],
        by = "failure_reason", all.x = TRUE, sort = FALSE
    )
    filter_summary[is.na(N), N := 0L]
    filter_summary[, failure_reason := factor(failure_reason, levels = failure_levels)]
    summary_output <- paste0(output_prefix, ".methylation.filter_summary.tsv")
    count_plot_output <- paste0(output_prefix, ".methylation.filter_counts.png")
    upset_plot_output <- paste0(output_prefix, ".methylation.filter_upset.png")
    fwrite(filter_summary, summary_output, sep = "\t")

    count_plot <- ggplot2::ggplot(filter_summary, ggplot2::aes(x = failure_reason, y = N, fill = failure_reason)) +
        ggplot2::geom_col(show.legend = FALSE) +
        ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.3, size = 3) +
        ggplot2::labs(title = "Cohort methylation site filtering", x = NULL, y = "Number of sites") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
    ggplot2::ggsave(count_plot_output, count_plot, width = 10, height = 6, dpi = 200, bg = "white")

    if (!requireNamespace("ggupset", quietly = TRUE)) {
        stop("Package 'ggupset' is required to write the methylation filter UpSet plot.")
    }

    upset_labels <- c("At least one missing/low-coverage call", "At least one extreme-coverage call", "Fails cohort sample-presence filter", "Fails methylation MAD filter")
    pattern_code <- as.integer(site_metadata$has_missing_or_low_coverage) +
        2L * as.integer(site_metadata$has_extreme_coverage_loss) +
        4L * as.integer(!site_metadata$pass_sample_presence_filter) +
        8L * as.integer(!site_metadata$pass_methylation_mad_filter)
    intersection_counts <- data.table(pattern_code = pattern_code)[, .N, by = pattern_code]
    setorder(intersection_counts, -N, pattern_code)
    intersection_counts[, filter_failures := lapply(
        pattern_code,
        function(code) upset_labels[bitwAnd(code, as.integer(2^(seq_along(upset_labels) - 1))) > 0]
    )]

    upset_plot <- ggplot2::ggplot(intersection_counts, ggplot2::aes(x = filter_failures, y = N)) +
        ggplot2::geom_col(fill = "#2C7FB8", width = 0.75) +
        ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.35, size = 3) +
        ggupset::scale_x_upset(sets = upset_labels, intersections = intersection_counts$filter_failures) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
        ggplot2::labs(title = "Overlap of cohort-level filter failures", x = NULL, y = "Number of sites") +
        ggplot2::theme_minimal(base_size = 12)
    ggplot2::ggsave(upset_plot_output, upset_plot, width = 12, height = 8, dpi = 200, bg = "white")
    message("Wrote filter summary: ", summary_output)
    message("Wrote filter-count plot: ", count_plot_output)
    message("Wrote filter UpSet plot: ", upset_plot_output)
}

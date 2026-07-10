#!/usr/bin/env Rscript

# Merge per-sample PacBio 5mC BED calls after coverage-based QC.
#
# The input manifest is a tab-delimited file with these columns:
#   sample_id    file_path
# Paths may be absolute or relative to the manifest's directory.
#
# Example:
# Rscript scripts/MergeMethylationCalls.R \
#   --InputManifest methylation_manifest.tsv \
#   --OutputPrefix results/cohort \
#   --MinCoverage 10 \
#   --MinSampleFraction 0.95 \
#   --ValueColumn mod_score
#
# The script writes:
#   <prefix>.methylation.filtered.long.tsv.gz  calls that passed all QC
#   <prefix>.methylation.site_qc.tsv.gz        per-site cohort QC summary
#   <prefix>.methylation.site_metadata.tsv.gz  all-site coverage/methylation metrics
#   <prefix>.methylation.sample_qc.tsv         per-sample QC summary
#   <prefix>.methylation.filter_summary.tsv     sequential cohort-filter counts
#   <prefix>.methylation.filter_counts.png      sequential cohort-filter bar chart
#   <prefix>.methylation.filter_upset.png       QC-condition intersection chart
#   <prefix>.methylation.raw.bed.gz             raw beta-value phenotype BED
#   <prefix>.methylation.INT.bed.gz             inverse-normal phenotype BED
#
# For sharded execution, run once per shard with --PerSampleOnly. Then call
# the script a final time with --FilteredCallList and --TotalSamples to apply
# the cohort-level site filter across every sample, not within each shard.

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

load_methylation_data <- function(file_path,
                                  filter_chroms = "X|Y|M|_",
                                  fence_k = 3) {
    # pb-CpG-tools prepends ## metadata lines before the #chrom header. Using
    # the header marker also supports already-cleaned BED files.
    loaded_data <- data.table::fread(file_path, skip = "#chrom")
    n_input_rows <- nrow(loaded_data)
    required_columns <- c("#chrom", "begin", "end", "mod_score", "type", "cov")
    missing_columns <- setdiff(required_columns, names(loaded_data))
    if (length(missing_columns) > 0) {
        stop(
            "Missing required column(s) in ", file_path, ": ",
            paste(missing_columns, collapse = ", "),
            ". Expected pb-CpG-tools columns #chrom, begin, end, mod_score, type, and cov."
        )
    }

    loaded_data[, cov := as.numeric(cov)]
    if (all(is.na(loaded_data$cov))) {
        stop("Column 'cov' is not numeric in ", file_path)
    }

    # Set --FilterChroms '' to retain every contig.
    if (!is.null(filter_chroms) && nzchar(filter_chroms)) {
        loaded_data <- loaded_data[!grepl(filter_chroms, `#chrom`)]
    }
    if (nrow(loaded_data) == 0) {
        stop("No rows remain after chromosome filtering in ", file_path)
    }

    call_types <- unique(as.character(loaded_data$type[!is.na(loaded_data$type)]))
    if (length(call_types) != 1) {
        stop(
            "Expected one pb-CpG-tools 'type' per input file, but found: ",
            paste(call_types, collapse = ", "), ". For a standard meQTL, use one .combined.bed.gz file per sample."
        )
    }

    median_cov <- median(loaded_data$cov, na.rm = TRUE)
    message("Median cov for ", basename(file_path), ": ", round(median_cov, 3))
    if (!is.finite(median_cov) || median_cov <= 0) {
        stop("Median coverage must be positive after chromosome filtering in ", file_path)
    }

    # Tukey's far-out fence on log10 coverage.  Zero-coverage rows are not
    # included in the log transform and will fail any positive MinCoverage.
    logc <- log10(loaded_data$cov[!is.na(loaded_data$cov) & loaded_data$cov > 0])
    if (length(logc) == 0) {
        extreme_cut <- Inf
        message("No positive coverage values in ", basename(file_path),
                "; setting extreme-coverage cutoff to Inf")
    } else {
        qs <- quantile(logc, c(0.25, 0.75), na.rm = TRUE, names = FALSE)
        extreme_cut <- 10^(qs[2] + fence_k * (qs[2] - qs[1]))
        message(
            "Extreme-coverage cutoff for ", basename(file_path), ": ",
            round(extreme_cut), "x (implied CN ~",
            round(2 * extreme_cut / median_cov), ")"
        )
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
    manifest <- data.table::fread(manifest_path)
    required_columns <- c("sample_id", "file_path")
    missing_columns <- setdiff(required_columns, names(manifest))
    if (length(missing_columns) > 0) {
        stop(
            "Input manifest must contain columns 'sample_id' and 'file_path'. Missing: ",
            paste(missing_columns, collapse = ", ")
        )
    }

    manifest <- manifest[, .(sample_id = as.character(sample_id), file_path = as.character(file_path))]
    if (anyNA(manifest$sample_id) || any(!nzchar(manifest$sample_id))) {
        stop("Input manifest contains an empty sample_id")
    }
    if (anyDuplicated(manifest$sample_id)) {
        stop("Each sample_id must occur exactly once in the input manifest")
    }
    if (anyNA(manifest$file_path) || any(!nzchar(manifest$file_path))) {
        stop("Input manifest contains an empty file_path")
    }

    manifest_dir <- dirname(normalizePath(manifest_path))
    is_absolute <- grepl("^(/|~)", manifest$file_path)
    manifest[!is_absolute, file_path := file.path(manifest_dir, file_path)]
    manifest[, file_path := path.expand(file_path)]
    missing_files <- manifest[!file.exists(file_path), file_path]
    if (length(missing_files) > 0) {
        stop("Input BED file(s) do not exist: ", paste(missing_files, collapse = ", "))
    }
    manifest
}

read_file_list <- function(list_path, label) {
    if (!file.exists(list_path)) {
        stop(label, " file list does not exist: ", list_path)
    }
    paths <- scan(list_path, what = character(), quiet = TRUE)
    if (length(paths) == 0) {
        stop(label, " file list is empty: ", list_path)
    }
    list_dir <- dirname(normalizePath(list_path))
    is_absolute <- grepl("^(/|~)", paths)
    paths[!is_absolute] <- file.path(list_dir, paths[!is_absolute])
    paths <- path.expand(paths)
    missing_paths <- paths[!file.exists(paths)]
    if (length(missing_paths) > 0) {
        stop(label, " file(s) do not exist: ", paste(missing_paths, collapse = ", "))
    }
    paths
}

read_call_tables <- function(call_paths, label, additional_required_columns = character()) {
    required_columns <- c("sample_id", "#chrom", "begin", "end", "cov", "site_key")
    required_columns <- c(required_columns, additional_required_columns)
    call_tables <- vector("list", length(call_paths))
    reference_columns <- NULL

    for (i in seq_along(call_paths)) {
        call_table <- data.table::fread(call_paths[i])
        missing_columns <- setdiff(required_columns, names(call_table))
        if (length(missing_columns) > 0) {
            stop(
                label, " file ", call_paths[i], " is missing: ",
                paste(missing_columns, collapse = ", ")
            )
        }
        call_table[, `:=`(
            sample_id = as.character(sample_id),
            cov = as.numeric(cov)
        )]
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
        stop(
            "Found duplicated sample/site calls across ", label, " shards. ",
            "Each sample must be present in exactly one shard and each site must occur once per sample."
        )
    }
    list(all_calls = all_calls, reference_columns = reference_columns)
}

read_sample_qc <- function(qc_paths, total_samples) {
    qc_tables <- lapply(qc_paths, data.table::fread)
    sample_qc <- rbindlist(qc_tables, use.names = TRUE)
    if (!("sample_id" %in% names(sample_qc))) {
        stop("Filtered sample-QC files must contain a sample_id column")
    }
    if (anyDuplicated(sample_qc$sample_id)) {
        stop("A sample appears in more than one filtered sample-QC file")
    }
    if (nrow(sample_qc) != total_samples) {
        stop(
            "Filtered sample-QC files contain ", nrow(sample_qc), " samples, but --TotalSamples is ",
            total_samples, ". These must match so the fraction threshold has the correct denominator."
        )
    }
    sample_qc
}

safe_mean <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else mean(x)
}

safe_sd <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 2) NA_real_ else sd(x)
}

safe_cv <- function(x) {
    average <- safe_mean(x)
    deviation <- safe_sd(x)
    if (is.na(average) || is.na(deviation) || average == 0) NA_real_ else deviation / average
}

safe_median <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else median(x)
}

safe_min <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else min(x)
}

safe_max <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else max(x)
}

safe_mad <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) NA_real_ else median(abs(x - median(x)))
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
    failure_levels <- c(
        "Insufficient minimum coverage",
        "Extreme coverage exclusion",
        "Low methylation MAD",
        "Pass all cohort filters"
    )
    filter_summary <- site_metadata[, .N, by = failure_reason]
    filter_summary <- merge(
        data.table(failure_reason = failure_levels),
        filter_summary,
        by = "failure_reason",
        all.x = TRUE,
        sort = FALSE
    )
    filter_summary[is.na(N), N := 0L]
    filter_summary[, failure_reason := factor(failure_reason, levels = failure_levels)]

    summary_output <- paste0(output_prefix, ".methylation.filter_summary.tsv")
    count_plot_output <- paste0(output_prefix, ".methylation.filter_counts.png")
    upset_plot_output <- paste0(output_prefix, ".methylation.filter_upset.png")
    fwrite(filter_summary, summary_output, sep = "\t")

    count_plot <- ggplot2::ggplot(
        filter_summary,
        ggplot2::aes(x = failure_reason, y = N, fill = failure_reason)
    ) +
        ggplot2::geom_col(show.legend = FALSE) +
        ggplot2::geom_text(ggplot2::aes(label = N), vjust = -0.3, size = 3) +
        ggplot2::labs(
            title = "Cohort methylation site filtering",
            x = NULL,
            y = "Number of sites"
        ) +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
    ggplot2::ggsave(
        count_plot_output,
        count_plot,
        width = 10,
        height = 6,
        dpi = 200,
        bg = "white"
    )

    upset_labels <- c(
        "At least one missing/low-coverage call",
        "At least one extreme-coverage call",
        "Fails cohort sample-presence filter",
        "Fails methylation MAD filter"
    )
    pattern_code <-
        as.integer(site_metadata$has_missing_or_low_coverage) +
        2L * as.integer(site_metadata$has_extreme_coverage_loss) +
        4L * as.integer(!site_metadata$pass_sample_presence_filter) +
        8L * as.integer(!site_metadata$pass_methylation_mad_filter)
    intersection_counts <- data.table(pattern_code = pattern_code)[, .N, by = pattern_code]
    setorder(intersection_counts, -N, pattern_code)
    intersection_counts[, intersection := factor(pattern_code, levels = pattern_code)]

    # Draw this simple UpSet-style plot with base graphics. This avoids a
    # layout-package dependency and keeps the script usable in lightweight R
    # environments while preserving the bar/matrix representation.
    grDevices::png(upset_plot_output, width = 2400, height = 1600, res = 200)
    graphics::layout(matrix(c(1, 2), ncol = 1), heights = c(2, 1.25))
    graphics::par(mar = c(2.5, 4.5, 3, 1))
    bar_positions <- graphics::barplot(
        intersection_counts$N,
        names.arg = rep("", nrow(intersection_counts)),
        col = "#2C7FB8",
        border = NA,
        ylab = "Number of sites",
        main = "Overlap of cohort-level filter failures"
    )
    graphics::text(
        bar_positions,
        intersection_counts$N,
        labels = intersection_counts$N,
        pos = 3,
        cex = 0.8
    )

    graphics::par(mar = c(4.5, 19, 0.5, 1))
    n_intersections <- nrow(intersection_counts)
    graphics::plot(
        NA,
        xlim = c(0.5, n_intersections + 0.5),
        ylim = c(0.5, length(upset_labels) + 0.5),
        xaxt = "n",
        yaxt = "n",
        xlab = "Filter-failure intersections",
        ylab = ""
    )
    graphics::axis(
        2,
        at = rev(seq_along(upset_labels)),
        labels = upset_labels,
        las = 1,
        tick = FALSE,
        cex.axis = 0.75
    )
    for (i in seq_len(n_intersections)) {
        y_positions <- rev(seq_along(upset_labels))
        graphics::points(rep(i, length(y_positions)), y_positions, pch = 16, col = "grey85", cex = 1.4)
        included_indices <- which(bitwAnd(
            intersection_counts$pattern_code[i],
            as.integer(2^(seq_along(upset_labels) - 1))
        ) > 0)
        if (length(included_indices) > 1) {
            included_y <- rev(seq_along(upset_labels))[included_indices]
            graphics::segments(i, min(included_y), i, max(included_y), col = "#2C7FB8", lwd = 2)
        }
        if (length(included_indices) > 0) {
            graphics::points(
                rep(i, length(included_indices)),
                rev(seq_along(upset_labels))[included_indices],
                pch = 16,
                col = "#2C7FB8",
                cex = 1.4
            )
        }
    }
    grDevices::dev.off()

    message("Wrote filter summary: ", summary_output)
    message("Wrote filter-count plot: ", count_plot_output)
    message("Wrote filter UpSet plot: ", upset_plot_output)
}

option_list <- list(
    make_option("--InputManifest", type = "character", default = NULL,
                help = "TSV with sample_id and file_path columns (normal/per-shard mode)"),
    make_option("--FilteredCallList", type = "character", default = NULL,
                help = "One filtered-call file path per line (final sharded-merge mode)"),
    make_option("--AllCallList", type = "character", default = NULL,
                help = "One per-sample-QC call file path per line for all-site metadata in final sharded-merge mode"),
    make_option("--FilteredSampleQcList", type = "character", default = NULL,
                help = "One per-shard sample-QC file path per line (required in final sharded-merge mode)"),
    make_option("--TotalSamples", type = "integer", default = 0,
                help = "Total input sample count for final sharded-merge mode [default: %default]"),
    make_option("--PerSampleOnly", action = "store_true", default = FALSE,
                help = "Apply only per-sample QC and write a shard intermediate"),
    make_option("--OutputPrefix", type = "character", default = NULL,
                help = "Prefix for output files [required]"),
    make_option("--MinCoverage", type = "double", default = 10,
                help = "Minimum per-call coverage to retain [default: %default]"),
    make_option("--MinSampleFraction", type = "double", default = 0.95,
                help = "Minimum fraction of all samples passing per-site QC [default: %default]"),
    make_option("--MinSamples", type = "integer", default = 0,
                help = "Additional minimum number of samples passing per-site QC [default: %default]"),
    make_option("--MinMethylationMAD", type = "double", default = 0.003,
                help = "Minimum cohort methylation MAD among per-sample-QC-passing calls [default: %default]"),
    make_option("--FilterChroms", type = "character", default = "X|Y|M|_",
                help = "Regex for chromosomes/contigs to remove; use '' to keep all [default: %default]"),
    make_option("--FenceK", type = "double", default = 3,
                help = "Tukey log10-coverage far-out fence multiplier [default: %default]"),
    make_option("--ValueColumn", type = "character", default = "mod_score",
                help = "pb-CpG-tools methylation column; mod_score is the recommended model-pileup value [default: %default]"),
    make_option("--ValueMultiplier", type = "double", default = 0.01,
                help = "Multiplier applied to ValueColumn before QTL output; 0.01 converts pb-CpG mod_score percent to beta values [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
has_manifest <- !is.null(opt$InputManifest)
has_filtered_call_list <- !is.null(opt$FilteredCallList)
if (is.null(opt$OutputPrefix) || (has_manifest == has_filtered_call_list)) {
    stop(
        "--OutputPrefix and exactly one of --InputManifest or --FilteredCallList are required. ",
        "Run with --help for usage."
    )
}
if (opt$PerSampleOnly && !has_manifest) {
    stop("--PerSampleOnly can only be used with --InputManifest")
}
if (has_filtered_call_list && opt$TotalSamples <= 0) {
    stop("--TotalSamples must be positive when using --FilteredCallList")
}
if (has_filtered_call_list && is.null(opt$AllCallList)) {
    stop("--AllCallList is required with --FilteredCallList to calculate all-site metadata")
}
if (has_filtered_call_list && is.null(opt$FilteredSampleQcList)) {
    stop(
        "--FilteredSampleQcList is required with --FilteredCallList to construct ",
        "a complete cohort-level QTL phenotype matrix"
    )
}
if (has_manifest && opt$TotalSamples != 0) {
    warning("--TotalSamples is ignored when --InputManifest is supplied")
}
if (has_manifest && !is.null(opt$AllCallList)) {
    warning("--AllCallList is ignored when --InputManifest is supplied")
}
if (!is.finite(opt$MinCoverage) || opt$MinCoverage < 0) {
    stop("--MinCoverage must be a non-negative number")
}
if (!is.finite(opt$MinSampleFraction) || opt$MinSampleFraction <= 0 || opt$MinSampleFraction > 1) {
    stop("--MinSampleFraction must be in (0, 1]")
}
if (is.na(opt$MinSamples) || opt$MinSamples < 0) {
    stop("--MinSamples must be a non-negative integer")
}
if (!is.finite(opt$MinMethylationMAD) || opt$MinMethylationMAD < 0) {
    stop("--MinMethylationMAD must be a non-negative number")
}
if (!is.finite(opt$FenceK) || opt$FenceK < 0) {
    stop("--FenceK must be a non-negative number")
}
if (!is.finite(opt$ValueMultiplier) || opt$ValueMultiplier <= 0) {
    stop("--ValueMultiplier must be a positive number")
}

sample_qc <- NULL
if (has_manifest) {
    manifest <- read_manifest(opt$InputManifest)
    n_samples <- nrow(manifest)
    cohort_sample_ids <- manifest$sample_id
    message("Processing ", n_samples, " samples")

    filtered_calls <- vector("list", n_samples)
    site_metric_calls <- vector("list", n_samples)
    sample_qc_tables <- vector("list", n_samples)
    reference_columns <- NULL

    for (i in seq_len(n_samples)) {
        sample_id <- manifest$sample_id[i]
        file_path <- manifest$file_path[i]
        message("[", i, "/", n_samples, "] Loading ", sample_id, ": ", file_path)

        methylation_data <- load_methylation_data(
            file_path = file_path,
            filter_chroms = opt$FilterChroms,
            fence_k = opt$FenceK
        )
        input_columns <- names(methylation_data)
        if (is.null(reference_columns)) {
            reference_columns <- copy(input_columns)
        } else if (!identical(reference_columns, input_columns)) {
            stop(
                "BED columns in sample ", sample_id,
                " do not match the first input file. All input BED files must use the same schema."
            )
        }
        if ("sample_id" %in% input_columns) {
            stop("Input BED files must not already contain a 'sample_id' column")
        }

        duplicate_sites <- methylation_data[, .N, by = site_key][N > 1]
        if (nrow(duplicate_sites) > 0) {
            stop(
                "Found ", nrow(duplicate_sites), " duplicated #chrom/begin/end site(s) in ",
                sample_id, ". Aggregate duplicate calls before merging so a site is counted once per sample."
            )
        }

        coverage_pass <- !is.na(methylation_data$cov) & methylation_data$cov >= opt$MinCoverage
        extreme_pass <- methylation_data$extreme_cov_flag == "ok"
        n_below_min_coverage <- sum(!coverage_pass)
        n_extreme_coverage <- sum(!extreme_pass)
        n_extreme_after_min_coverage <- sum(coverage_pass & !extreme_pass)
        current_sample_id <- sample_id
        methylation_data[, `:=`(
            sample_id = current_sample_id,
            meets_min_coverage = coverage_pass,
            per_sample_qc_pass = coverage_pass & extreme_pass
        )]
        retained <- methylation_data[per_sample_qc_pass == TRUE]

        sample_qc_tables[[i]] <- data.table(
            sample_id = sample_id,
            file_path = file_path,
            n_input_rows = attr(methylation_data, "n_input_rows"),
            n_rows_after_chrom_filter = nrow(methylation_data),
            n_removed_by_chrom_filter = attr(methylation_data, "n_input_rows") - nrow(methylation_data),
            pb_cpg_type = attr(methylation_data, "call_type"),
            median_cov = attr(methylation_data, "median_cov"),
            extreme_coverage_cutoff = attr(methylation_data, "extreme_cut"),
            n_below_min_coverage = n_below_min_coverage,
            n_extreme_coverage = n_extreme_coverage,
            n_extreme_coverage_after_min_coverage = n_extreme_after_min_coverage,
            n_passing_per_sample_qc = nrow(retained)
        )
        filtered_calls[[i]] <- retained
        site_metric_calls[[i]] <- methylation_data
        message(
            "  Input sites: ", attr(methylation_data, "n_input_rows"),
            "; removed by chromosome filter: ",
            attr(methylation_data, "n_input_rows") - nrow(methylation_data),
            "; evaluated for coverage: ", nrow(methylation_data)
        )
        message(
            "  Per-sample thresholds: ", n_below_min_coverage,
            " fail MinCoverage (<", opt$MinCoverage, "); ",
            n_extreme_after_min_coverage,
            " fail extreme coverage after MinCoverage; ",
            nrow(retained), " pass both thresholds"
        )
    }
    all_calls <- rbindlist(filtered_calls, use.names = TRUE)
    all_site_calls <- rbindlist(site_metric_calls, use.names = TRUE)
    sample_qc <- rbindlist(sample_qc_tables, use.names = TRUE)
} else {
    filtered_call_paths <- read_file_list(opt$FilteredCallList, "Filtered-call")
    filtered_call_data <- read_call_tables(filtered_call_paths, "Filtered-call")
    all_calls <- filtered_call_data$all_calls
    reference_columns <- filtered_call_data$reference_columns
    all_call_paths <- read_file_list(opt$AllCallList, "All-call")
    all_call_data <- read_call_tables(
        all_call_paths,
        "All-call",
        c("meets_min_coverage", "per_sample_qc_pass")
    )
    all_site_calls <- all_call_data$all_calls
    all_site_calls[, `:=`(
        meets_min_coverage = as.logical(meets_min_coverage),
        per_sample_qc_pass = as.logical(per_sample_qc_pass)
    )]
    if (anyNA(all_site_calls$meets_min_coverage) || anyNA(all_site_calls$per_sample_qc_pass)) {
        stop("All-call files contain non-logical meets_min_coverage or per_sample_qc_pass values")
    }
    n_samples <- opt$TotalSamples

    observed_samples <- unique(all_calls$sample_id)
    if (length(observed_samples) > n_samples) {
        stop(
            "Filtered-call files contain ", length(observed_samples),
            " samples, which exceeds --TotalSamples (", n_samples, ")"
        )
    }
    if (!is.null(opt$FilteredSampleQcList)) {
        filtered_sample_qc_paths <- read_file_list(opt$FilteredSampleQcList, "Filtered sample-QC")
        sample_qc <- read_sample_qc(filtered_sample_qc_paths, n_samples)
        missing_qc_samples <- setdiff(observed_samples, sample_qc$sample_id)
        if (length(missing_qc_samples) > 0) {
            stop(
                "Filtered-call files contain sample(s) not present in the filtered sample-QC files: ",
                paste(missing_qc_samples, collapse = ", ")
            )
        }
    }
    cohort_sample_ids <- sample_qc$sample_id
    message("Reading per-sample-QC-passing calls and all-site metadata from ",
            length(filtered_call_paths), " shard(s) for ", n_samples, " total samples")
}

output_dir <- dirname(opt$OutputPrefix)
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}
sample_qc_output <- paste0(opt$OutputPrefix, ".methylation.sample_qc.tsv")

if (opt$PerSampleOnly) {
    per_sample_output <- paste0(opt$OutputPrefix, ".methylation.per_sample_filtered.long.tsv.gz")
    all_call_output <- paste0(opt$OutputPrefix, ".methylation.per_sample_qc.long.tsv.gz")
    fwrite(all_calls, per_sample_output, sep = "\t", na = "NA")
    fwrite(all_site_calls, all_call_output, sep = "\t", na = "NA")
    fwrite(sample_qc, sample_qc_output, sep = "\t", na = "NA")
    message("Wrote per-sample-QC-passing shard calls: ", per_sample_output)
    message("Wrote all per-sample-QC shard calls for site metadata: ", all_call_output)
    message("Wrote sample QC: ", sample_qc_output)
    quit(save = "no", status = 0)
}

required_samples <- max(ceiling(n_samples * opt$MinSampleFraction), opt$MinSamples)
message(
    "A site must pass per-sample QC in at least ", required_samples, " of ", n_samples,
    " samples (fraction threshold = ", opt$MinSampleFraction,
    "; count threshold = ", opt$MinSamples, ")"
)

if (!is.null(opt$ValueColumn)) {
    if (!(opt$ValueColumn %in% names(all_site_calls))) {
        stop(
            "--ValueColumn '", opt$ValueColumn,
            "' is not present in the input BED files. Available columns: ",
            paste(reference_columns, collapse = ", ")
        )
    }
    methylation_values <- suppressWarnings(as.numeric(all_site_calls[[opt$ValueColumn]]))
    if (all(is.na(methylation_values)) && any(!is.na(all_site_calls[[opt$ValueColumn]]))) {
        stop("--ValueColumn '", opt$ValueColumn, "' must be numeric")
    }
    all_site_calls[, methylation_value_for_metrics := methylation_values * opt$ValueMultiplier]
} else {
    all_site_calls[, methylation_value_for_metrics := NA_real_]
}

# All calls remaining after chromosome filtering are represented here, including
# sites that fail MinCoverage, extreme-coverage QC, or cohort-level QC.
site_metadata <- all_site_calls[, {
    per_sample_pass <- per_sample_qc_pass == TRUE
    list(
        n_samples_observed = uniqueN(sample_id),
        fraction_samples_observed = uniqueN(sample_id) / n_samples,
        mean_cov_all_calls = safe_mean(cov),
        sd_cov_all_calls = safe_sd(cov),
        cv_cov_all_calls = safe_cv(cov),
        mean_methylation_all_calls = safe_mean(methylation_value_for_metrics),
        sd_methylation_all_calls = safe_sd(methylation_value_for_metrics),
        cv_methylation_all_calls = safe_cv(methylation_value_for_metrics),
        n_samples_min_coverage = sum(meets_min_coverage == TRUE),
        fraction_samples_min_coverage = sum(meets_min_coverage == TRUE) / n_samples,
        n_samples_passing_per_sample_qc = sum(per_sample_pass),
        fraction_samples_passing_per_sample_qc = sum(per_sample_pass) / n_samples,
        mean_cov_passing_per_sample_qc = safe_mean(cov[per_sample_pass]),
        sd_cov_passing_per_sample_qc = safe_sd(cov[per_sample_pass]),
        cv_cov_passing_per_sample_qc = safe_cv(cov[per_sample_pass]),
        median_cov_passing_per_sample_qc = safe_median(cov[per_sample_pass]),
        min_cov_passing_per_sample_qc = safe_min(cov[per_sample_pass]),
        max_cov_passing_per_sample_qc = safe_max(cov[per_sample_pass]),
        mean_methylation_passing_per_sample_qc = safe_mean(methylation_value_for_metrics[per_sample_pass]),
        sd_methylation_passing_per_sample_qc = safe_sd(methylation_value_for_metrics[per_sample_pass]),
        cv_methylation_passing_per_sample_qc = safe_cv(methylation_value_for_metrics[per_sample_pass]),
        methylation_mad_passing_per_sample_qc = safe_mad(methylation_value_for_metrics[per_sample_pass])
    )
}, by = .(`#chrom`, begin, end, site_key)]
site_metadata[, `:=`(
    n_samples_required = required_samples,
    pass_minimum_coverage_filter = n_samples_min_coverage >= required_samples,
    pass_sample_presence_filter = n_samples_passing_per_sample_qc >= required_samples,
    pass_methylation_mad_filter = !is.na(methylation_mad_passing_per_sample_qc) &
        methylation_mad_passing_per_sample_qc >= opt$MinMethylationMAD
)]
site_metadata[, `:=`(
    has_missing_or_low_coverage = n_samples_min_coverage < n_samples,
    has_extreme_coverage_loss = n_samples_passing_per_sample_qc < n_samples_min_coverage,
    keep_site = pass_sample_presence_filter & pass_methylation_mad_filter
)]
site_metadata[, failure_reason := fcase(
    !pass_minimum_coverage_filter, "Insufficient minimum coverage",
    !pass_sample_presence_filter, "Extreme coverage exclusion",
    !pass_methylation_mad_filter, "Low methylation MAD",
    default = "Pass all cohort filters"
)]
setorder(site_metadata, `#chrom`, begin, end)

n_sites_total <- nrow(site_metadata)
n_sites_with_min_coverage <- site_metadata[n_samples_min_coverage > 0, .N]
n_sites_with_per_sample_qc <- site_metadata[n_samples_passing_per_sample_qc > 0, .N]
n_sites_failing_sample_presence <- site_metadata[pass_sample_presence_filter == FALSE, .N]
n_sites_passing_sample_presence <- site_metadata[pass_sample_presence_filter == TRUE, .N]
n_sites_failing_mad <- site_metadata[
    pass_sample_presence_filter == TRUE & pass_methylation_mad_filter == FALSE,
    .N
]
n_sites_passing_cohort_qc <- site_metadata[keep_site == TRUE, .N]
message(
    "Cohort site summary: ", n_sites_total,
    " sites observed after chromosome filtering; ",
    n_sites_with_min_coverage, " have >=1 sample meeting MinCoverage; ",
    n_sites_with_per_sample_qc, " have >=1 sample passing per-sample QC"
)
message(
    "Cohort sample-presence threshold: ", n_sites_failing_sample_presence,
    " fail the required ", required_samples, "/", n_samples,
    " sample threshold; ", n_sites_passing_sample_presence, " sites pass"
)
message(
    "Cohort methylation MAD filter: ", n_sites_failing_mad,
    " sample-presence-passing site(s) fail MAD < ", opt$MinMethylationMAD,
    "; ", n_sites_passing_cohort_qc, " sites pass all cohort filters"
)

kept_site_keys <- site_metadata[keep_site == TRUE, site_key]
merged_calls <- all_calls[site_key %chin% kept_site_keys]
setcolorder(merged_calls, c("sample_id", setdiff(names(merged_calls), "sample_id")))

long_output <- paste0(opt$OutputPrefix, ".methylation.filtered.long.tsv.gz")
site_qc_output <- paste0(opt$OutputPrefix, ".methylation.site_qc.tsv.gz")
site_metadata_output <- paste0(opt$OutputPrefix, ".methylation.site_metadata.tsv.gz")
raw_bed_output <- paste0(opt$OutputPrefix, ".methylation.raw.bed.gz")
int_bed_output <- paste0(opt$OutputPrefix, ".methylation.INT.bed.gz")

site_metadata[, n_samples_imputed_in_qtl_bed := 0L]

if (!is.null(opt$ValueColumn)) {
    merged_calls[, methylation_value_for_qtl :=
        suppressWarnings(as.numeric(get(opt$ValueColumn))) * opt$ValueMultiplier]
    matrix_formula <- as.formula("`#chrom` + begin + end + site_key ~ sample_id")
    raw_methylation_bed <- dcast(
        merged_calls,
        formula = matrix_formula,
        value.var = "methylation_value_for_qtl"
    )
    setorder(raw_methylation_bed, `#chrom`, begin, end)
    setnames(
        raw_methylation_bed,
        c("#chrom", "begin", "end", "site_key"),
        c("#chr", "start", "end", "phenotype_id")
    )
    phenotype_columns <- c("#chr", "start", "end", "phenotype_id")
    missing_sample_columns <- setdiff(cohort_sample_ids, names(raw_methylation_bed))
    if (length(missing_sample_columns) > 0) {
        raw_methylation_bed[, (missing_sample_columns) := NA_real_]
    }
    setcolorder(raw_methylation_bed, c(phenotype_columns, cohort_sample_ids))

    sample_columns <- cohort_sample_ids
    raw_values <- as.matrix(raw_methylation_bed[, ..sample_columns])
    n_samples_imputed <- rowSums(is.na(raw_values))
    if (nrow(raw_methylation_bed) > 0 && length(sample_columns) > 0) {
        for (row_index in which(n_samples_imputed > 0)) {
            feature_mean <- mean(raw_values[row_index, ], na.rm = TRUE)
            if (!is.finite(feature_mean)) {
                stop(
                    "Cannot impute a retained QTL feature with no observed methylation values: ",
                    raw_methylation_bed$phenotype_id[[row_index]]
                )
            }
            raw_values[row_index, is.na(raw_values[row_index, ])] <- feature_mean
        }
        for (column_index in seq_along(sample_columns)) {
            set(
                raw_methylation_bed,
                j = sample_columns[[column_index]],
                value = raw_values[, column_index]
            )
        }
    }

    imputation_summary <- data.table(
        site_key = raw_methylation_bed$phenotype_id,
        n_samples_imputed_in_qtl_bed = as.integer(n_samples_imputed)
    )
    site_metadata[
        imputation_summary,
        on = .(site_key),
        n_samples_imputed_in_qtl_bed := i.n_samples_imputed_in_qtl_bed
    ]
    message(
        "Cohort mean imputation: ", sum(n_samples_imputed), " sample/site value(s) imputed ",
        "across ", sum(n_samples_imputed > 0), " retained QTL feature(s)"
    )

    int_methylation_bed <- copy(raw_methylation_bed)
    if (nrow(int_methylation_bed) > 0 && length(sample_columns) > 0) {
        int_values <- t(vapply(
            seq_len(nrow(raw_methylation_bed)),
            function(row_index) inverse_normal_transform(as.numeric(raw_values[row_index, ])),
            FUN.VALUE = numeric(length(sample_columns))
        ))
        for (column_index in seq_along(sample_columns)) {
            set(
                int_methylation_bed,
                j = sample_columns[[column_index]],
                value = int_values[, column_index]
            )
        }
    }
}

# Preserve the compact QC output while the metadata output contains the
# all-call and passing-call coverage/methylation summaries and imputation count.
site_qc <- site_metadata[, .(
    `#chrom`, begin, end, site_key,
    n_samples_passing = n_samples_passing_per_sample_qc,
    fraction_samples_passing = fraction_samples_passing_per_sample_qc,
    median_cov_passing = median_cov_passing_per_sample_qc,
    min_cov_passing = min_cov_passing_per_sample_qc,
    max_cov_passing = max_cov_passing_per_sample_qc,
    n_samples_required,
    keep_site
)]
setorder(site_qc, `#chrom`, begin, end)

fwrite(merged_calls, long_output, sep = "\t", na = "NA")
fwrite(site_qc, site_qc_output, sep = "\t", na = "NA")
fwrite(site_metadata, site_metadata_output, sep = "\t", na = "NA")
if (!is.null(sample_qc)) {
    fwrite(sample_qc, sample_qc_output, sep = "\t", na = "NA")
}
if (!is.null(opt$ValueColumn)) {
    fwrite(raw_methylation_bed, raw_bed_output, sep = "\t", na = "NA")
    fwrite(int_methylation_bed, int_bed_output, sep = "\t", na = "NA")
}
write_filter_plots(site_metadata, opt$OutputPrefix)
message("Kept ", length(kept_site_keys), " / ", nrow(site_qc), " sites after cohort-level QC")
message("Wrote filtered long calls: ", long_output)
message("Wrote site QC: ", site_qc_output)
message("Wrote all-site metadata: ", site_metadata_output)
if (!is.null(sample_qc)) {
    message("Wrote sample QC: ", sample_qc_output)
}
if (!is.null(opt$ValueColumn)) {
    message("Wrote TensorQTL-compatible raw beta-value BED: ", raw_bed_output)
    message("Wrote TensorQTL-compatible inverse-normal BED: ", int_bed_output)
}

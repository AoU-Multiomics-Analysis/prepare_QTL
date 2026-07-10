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
#   --MinSampleFraction 0.8 \
#   --ValueColumn methylation_fraction
#
# The script writes:
#   <prefix>.methylation.filtered.long.tsv.gz  calls that passed all QC
#   <prefix>.methylation.site_qc.tsv.gz        per-site cohort QC summary
#   <prefix>.methylation.site_metadata.tsv.gz  all-site coverage/methylation metrics
#   <prefix>.methylation.sample_qc.tsv         per-sample QC summary
#   <prefix>.methylation.matrix.bed.gz         optional site-by-sample matrix
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
    loaded_data <- data.table::fread(file_path)
    n_input_rows <- nrow(loaded_data)
    required_columns <- c("#chrom", "begin", "end", "cov")
    missing_columns <- setdiff(required_columns, names(loaded_data))
    if (length(missing_columns) > 0) {
        stop(
            "Missing required column(s) in ", file_path, ": ",
            paste(missing_columns, collapse = ", "),
            ". Expected #chrom, begin, end, and cov."
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

option_list <- list(
    make_option("--InputManifest", type = "character", default = NULL,
                help = "TSV with sample_id and file_path columns (normal/per-shard mode)"),
    make_option("--FilteredCallList", type = "character", default = NULL,
                help = "One filtered-call file path per line (final sharded-merge mode)"),
    make_option("--AllCallList", type = "character", default = NULL,
                help = "One per-sample-QC call file path per line for all-site metadata in final sharded-merge mode"),
    make_option("--FilteredSampleQcList", type = "character", default = NULL,
                help = "Optional one per-shard sample-QC file path per line (final sharded-merge mode)"),
    make_option("--TotalSamples", type = "integer", default = 0,
                help = "Total input sample count for final sharded-merge mode [default: %default]"),
    make_option("--PerSampleOnly", action = "store_true", default = FALSE,
                help = "Apply only per-sample QC and write a shard intermediate"),
    make_option("--OutputPrefix", type = "character", default = NULL,
                help = "Prefix for output files [required]"),
    make_option("--MinCoverage", type = "double", default = 10,
                help = "Minimum per-call coverage to retain [default: %default]"),
    make_option("--MinSampleFraction", type = "double", default = 0.8,
                help = "Minimum fraction of all samples passing per-site QC [default: %default]"),
    make_option("--MinSamples", type = "integer", default = 0,
                help = "Additional minimum number of samples passing per-site QC [default: %default]"),
    make_option("--FilterChroms", type = "character", default = "X|Y|M|_",
                help = "Regex for chromosomes/contigs to remove; use '' to keep all [default: %default]"),
    make_option("--FenceK", type = "double", default = 3,
                help = "Tukey log10-coverage far-out fence multiplier [default: %default]"),
    make_option("--ValueColumn", type = "character", default = NULL,
                help = "Optional call column to pivot into the site-by-sample BED matrix")
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
if (!is.finite(opt$FenceK) || opt$FenceK < 0) {
    stop("--FenceK must be a non-negative number")
}

sample_qc <- NULL
if (has_manifest) {
    manifest <- read_manifest(opt$InputManifest)
    n_samples <- nrow(manifest)
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
            median_cov = attr(methylation_data, "median_cov"),
            extreme_coverage_cutoff = attr(methylation_data, "extreme_cut"),
            n_below_min_coverage = sum(!coverage_pass),
            n_extreme_coverage = sum(!extreme_pass),
            n_passing_per_sample_qc = nrow(retained)
        )
        filtered_calls[[i]] <- retained
        site_metric_calls[[i]] <- methylation_data
        message("  Retained ", nrow(retained), " / ", nrow(methylation_data),
                " calls after per-sample coverage QC")
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
    all_site_calls[, methylation_value_for_metrics := methylation_values]
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
        cv_methylation_passing_per_sample_qc = safe_cv(methylation_value_for_metrics[per_sample_pass])
    )
}, by = .(`#chrom`, begin, end, site_key)]
site_metadata[, `:=`(
    n_samples_required = required_samples,
    keep_site = n_samples_passing_per_sample_qc >= required_samples
)]
setorder(site_metadata, `#chrom`, begin, end)

# Preserve the compact QC output while the metadata output below contains the
# all-call and passing-call coverage/methylation summaries.
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

kept_site_keys <- site_metadata[keep_site == TRUE, site_key]
merged_calls <- all_calls[site_key %chin% kept_site_keys]
setcolorder(merged_calls, c("sample_id", setdiff(names(merged_calls), "sample_id")))

long_output <- paste0(opt$OutputPrefix, ".methylation.filtered.long.tsv.gz")
site_qc_output <- paste0(opt$OutputPrefix, ".methylation.site_qc.tsv.gz")
site_metadata_output <- paste0(opt$OutputPrefix, ".methylation.site_metadata.tsv.gz")

fwrite(merged_calls, long_output, sep = "\t", na = "NA")
fwrite(site_qc, site_qc_output, sep = "\t", na = "NA")
fwrite(site_metadata, site_metadata_output, sep = "\t", na = "NA")
if (!is.null(sample_qc)) {
    fwrite(sample_qc, sample_qc_output, sep = "\t", na = "NA")
}
message("Kept ", length(kept_site_keys), " / ", nrow(site_qc), " sites after cohort-level QC")
message("Wrote filtered long calls: ", long_output)
message("Wrote site QC: ", site_qc_output)
message("Wrote all-site metadata: ", site_metadata_output)
if (!is.null(sample_qc)) {
    message("Wrote sample QC: ", sample_qc_output)
}

if (!is.null(opt$ValueColumn)) {
    matrix_output <- paste0(opt$OutputPrefix, ".methylation.matrix.bed.gz")
    matrix_formula <- as.formula("`#chrom` + begin + end + site_key ~ sample_id")
    methylation_matrix <- dcast(
        merged_calls,
        formula = matrix_formula,
        value.var = opt$ValueColumn
    )
    setorder(methylation_matrix, `#chrom`, begin, end)
    fwrite(methylation_matrix, matrix_output, sep = "\t", na = "NA")
    message("Wrote site-by-sample methylation matrix: ", matrix_output)
}

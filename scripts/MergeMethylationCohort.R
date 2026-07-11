#!/usr/bin/env Rscript

# Reduce shard-level pb-CpG call tables to cohort-wide QC outputs and QTL BEDs.

suppressPackageStartupMessages(library(optparse))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_file)), "MethylationUtils.R"))

option_list <- list(
    make_option("--AllCallList", type = "character", help = "One per-shard all-call file path per line [required]"),
    make_option("--CohortSamples", type = "character", help = "TSV containing one sample_id column for the complete cohort [required]"),
    make_option("--TotalSamples", type = "integer", help = "Total input sample count [required]"),
    make_option("--OutputPrefix", type = "character", help = "Prefix for output files [required]"),
    make_option("--Chromosome", type = "character", default = "",
                help = "Optional chromosome name to merge from the input call files [default: all chromosomes]"),
    make_option("--MinSampleFraction", type = "double", default = 0.95,
                help = "Minimum fraction of all samples passing per-site QC [default: %default]"),
    make_option("--MinSamples", type = "integer", default = 0,
                help = "Additional minimum number of samples passing per-site QC [default: %default]"),
    make_option("--MinMethylationMAD", type = "double", default = 0.003,
                help = "Minimum cohort methylation MAD among per-sample-QC-passing calls [default: %default]"),
    make_option("--ValueColumn", type = "character", default = "mod_score",
                help = "pb-CpG-tools methylation column [default: %default]"),
    make_option("--ValueMultiplier", type = "double", default = 0.01,
                help = "Multiplier applied to ValueColumn [default: %default]"),
    make_option("--SkipFilterPlots", action = "store_true", default = FALSE,
                help = "Do not write filter summary/count/UpSet outputs")
)
opt <- parse_args(OptionParser(option_list = option_list))
required_options <- c("AllCallList", "CohortSamples", "OutputPrefix")
if (any(vapply(required_options, function(name) is.null(opt[[name]]), logical(1)))) {
    stop("--AllCallList, --CohortSamples, and --OutputPrefix are required")
}
if (is.null(opt$TotalSamples) || opt$TotalSamples <= 0) stop("--TotalSamples must be positive")
if (!is.finite(opt$MinSampleFraction) || opt$MinSampleFraction <= 0 || opt$MinSampleFraction > 1) stop("--MinSampleFraction must be in (0, 1]")
if (is.na(opt$MinSamples) || opt$MinSamples < 0) stop("--MinSamples must be a non-negative integer")
if (!is.finite(opt$MinMethylationMAD) || opt$MinMethylationMAD < 0) stop("--MinMethylationMAD must be a non-negative number")
if (!is.finite(opt$ValueMultiplier) || opt$ValueMultiplier <= 0) stop("--ValueMultiplier must be a positive number")

all_call_paths <- read_file_list(opt$AllCallList, "All-call")
all_call_data <- read_call_tables(all_call_paths, "All-call", c("meets_min_coverage", "per_sample_qc_pass"))
all_site_calls <- all_call_data$all_calls
reference_columns <- all_call_data$reference_columns
if (nzchar(opt$Chromosome)) {
    all_site_calls <- all_site_calls[`#chrom` == opt$Chromosome]
}
all_site_calls[, `:=`(
    meets_min_coverage = as.logical(meets_min_coverage),
    per_sample_qc_pass = as.logical(per_sample_qc_pass)
)]
if (anyNA(all_site_calls$meets_min_coverage) || anyNA(all_site_calls$per_sample_qc_pass)) {
    stop("All-call files contain non-logical meets_min_coverage or per_sample_qc_pass values")
}
all_calls <- all_site_calls[per_sample_qc_pass == TRUE]
n_samples <- opt$TotalSamples
cohort_samples <- fread(opt$CohortSamples, colClasses = "character")
if (!identical(names(cohort_samples), "sample_id")) stop("--CohortSamples must contain exactly one column named sample_id")
if (anyNA(cohort_samples$sample_id) || any(!nzchar(cohort_samples$sample_id))) stop("--CohortSamples contains an empty sample_id")
if (anyDuplicated(cohort_samples$sample_id)) stop("--CohortSamples contains duplicate sample_id values")
if (nrow(cohort_samples) != n_samples) {
    stop("--CohortSamples contains ", nrow(cohort_samples), " samples, but --TotalSamples is ", n_samples)
}
observed_samples <- unique(all_calls$sample_id)
if (length(observed_samples) > n_samples) stop("Filtered-call files contain more samples than --TotalSamples")
unknown_samples <- setdiff(observed_samples, cohort_samples$sample_id)
if (length(unknown_samples) > 0) {
    stop("All-call files contain sample(s) not present in --CohortSamples: ",
         paste(unknown_samples, collapse = ", "))
}
cohort_sample_ids <- cohort_samples$sample_id
chromosome_label <- if (nzchar(opt$Chromosome)) paste0(" for ", opt$Chromosome) else ""
message("Reading per-sample-QC-passing calls and all-site metadata from ",
        length(all_call_paths), " shard(s)", chromosome_label, " for ", n_samples, " total samples")

output_dir <- dirname(opt$OutputPrefix)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
required_samples <- max(ceiling(n_samples * opt$MinSampleFraction), opt$MinSamples)
message("A site must pass per-sample QC in at least ", required_samples, " of ", n_samples,
        " samples (fraction threshold = ", opt$MinSampleFraction, "; count threshold = ", opt$MinSamples, ")")

if (!(opt$ValueColumn %in% names(all_site_calls))) {
    stop("--ValueColumn '", opt$ValueColumn, "' is not present in the input BED files. Available columns: ",
         paste(reference_columns, collapse = ", "))
}
methylation_values <- suppressWarnings(as.numeric(all_site_calls[[opt$ValueColumn]]))
if (all(is.na(methylation_values)) && any(!is.na(all_site_calls[[opt$ValueColumn]]))) {
    stop("--ValueColumn '", opt$ValueColumn, "' must be numeric")
}
all_site_calls[, methylation_value_for_metrics := methylation_values * opt$ValueMultiplier]

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
    pass_methylation_mad_filter = !is.na(methylation_mad_passing_per_sample_qc) & methylation_mad_passing_per_sample_qc >= opt$MinMethylationMAD
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
n_sites_failing_mad <- site_metadata[pass_sample_presence_filter == TRUE & pass_methylation_mad_filter == FALSE, .N]
n_sites_passing_cohort_qc <- site_metadata[keep_site == TRUE, .N]
message("Cohort site summary: ", n_sites_total, " sites observed after chromosome filtering; ",
        n_sites_with_min_coverage, " have >=1 sample meeting MinCoverage; ",
        n_sites_with_per_sample_qc, " have >=1 sample passing per-sample QC")
message("Cohort sample-presence threshold: ", n_sites_failing_sample_presence,
        " fail the required ", required_samples, "/", n_samples, " sample threshold; ",
        n_sites_passing_sample_presence, " sites pass")
message("Cohort methylation MAD filter: ", n_sites_failing_mad,
        " sample-presence-passing site(s) fail MAD < ", opt$MinMethylationMAD,
        "; ", n_sites_passing_cohort_qc, " sites pass all cohort filters")

kept_site_keys <- site_metadata[keep_site == TRUE, site_key]
merged_calls <- all_calls[site_key %chin% kept_site_keys]
setcolorder(merged_calls, c("sample_id", setdiff(names(merged_calls), "sample_id")))
long_output <- paste0(opt$OutputPrefix, ".methylation.filtered.long.tsv.gz")
site_qc_output <- paste0(opt$OutputPrefix, ".methylation.site_qc.tsv.gz")
site_metadata_output <- paste0(opt$OutputPrefix, ".methylation.site_metadata.tsv.gz")
raw_bed_output <- paste0(opt$OutputPrefix, ".methylation.raw.bed.gz")
int_bed_output <- paste0(opt$OutputPrefix, ".methylation.INT.bed.gz")

site_metadata[, n_samples_imputed_in_qtl_bed := 0L]
merged_calls[, methylation_value_for_qtl := suppressWarnings(as.numeric(get(opt$ValueColumn))) * opt$ValueMultiplier]
matrix_formula <- as.formula("`#chrom` + begin + end + site_key ~ sample_id")
phenotype_columns <- c("#chr", "start", "end", "phenotype_id")
if (nrow(merged_calls) == 0) {
    raw_methylation_bed <- data.table(
        `#chr` = character(),
        start = integer(),
        end = integer(),
        phenotype_id = character()
    )
} else {
    raw_methylation_bed <- dcast(merged_calls, formula = matrix_formula, value.var = "methylation_value_for_qtl")
    setorder(raw_methylation_bed, `#chrom`, begin, end)
    setnames(raw_methylation_bed, c("#chrom", "begin", "end", "site_key"), phenotype_columns)
}
missing_sample_columns <- setdiff(cohort_sample_ids, names(raw_methylation_bed))
if (length(missing_sample_columns) > 0) raw_methylation_bed[, (missing_sample_columns) := NA_real_]
setcolorder(raw_methylation_bed, c(phenotype_columns, cohort_sample_ids))

sample_columns <- cohort_sample_ids
raw_values <- as.matrix(raw_methylation_bed[, ..sample_columns])
n_samples_imputed <- rowSums(is.na(raw_values))
if (nrow(raw_methylation_bed) > 0 && length(sample_columns) > 0) {
    for (row_index in which(n_samples_imputed > 0)) {
        feature_mean <- mean(raw_values[row_index, ], na.rm = TRUE)
        if (!is.finite(feature_mean)) stop("Cannot impute a retained QTL feature with no observed methylation values: ", raw_methylation_bed$phenotype_id[[row_index]])
        raw_values[row_index, is.na(raw_values[row_index, ])] <- feature_mean
    }
    for (column_index in seq_along(sample_columns)) {
        set(raw_methylation_bed, j = sample_columns[[column_index]], value = raw_values[, column_index])
    }
}
imputation_summary <- data.table(site_key = raw_methylation_bed$phenotype_id, n_samples_imputed_in_qtl_bed = as.integer(n_samples_imputed))
site_metadata[imputation_summary, on = .(site_key), n_samples_imputed_in_qtl_bed := i.n_samples_imputed_in_qtl_bed]
message("Cohort mean imputation: ", sum(n_samples_imputed), " sample/site value(s) imputed across ",
        sum(n_samples_imputed > 0), " retained QTL feature(s)")

int_methylation_bed <- copy(raw_methylation_bed)
if (nrow(int_methylation_bed) > 0 && length(sample_columns) > 0) {
    int_values <- t(vapply(seq_len(nrow(raw_methylation_bed)),
                           function(row_index) inverse_normal_transform(as.numeric(raw_values[row_index, ])),
                           FUN.VALUE = numeric(length(sample_columns))))
    for (column_index in seq_along(sample_columns)) {
        set(int_methylation_bed, j = sample_columns[[column_index]], value = int_values[, column_index])
    }
}

site_qc <- site_metadata[, .(
    `#chrom`, begin, end, site_key,
    n_samples_passing = n_samples_passing_per_sample_qc,
    fraction_samples_passing = fraction_samples_passing_per_sample_qc,
    median_cov_passing = median_cov_passing_per_sample_qc,
    min_cov_passing = min_cov_passing_per_sample_qc,
    max_cov_passing = max_cov_passing_per_sample_qc,
    n_samples_required, keep_site
)]
setorder(site_qc, `#chrom`, begin, end)
fwrite(merged_calls, long_output, sep = "\t", na = "NA")
fwrite(site_qc, site_qc_output, sep = "\t", na = "NA")
fwrite(site_metadata, site_metadata_output, sep = "\t", na = "NA")
fwrite(raw_methylation_bed, raw_bed_output, sep = "\t", na = "NA")
fwrite(int_methylation_bed, int_bed_output, sep = "\t", na = "NA")
if (!opt$SkipFilterPlots) write_filter_plots(site_metadata, opt$OutputPrefix)
message("Kept ", length(kept_site_keys), " / ", nrow(site_qc), " sites after cohort-level QC")
message("Wrote filtered long calls: ", long_output)
message("Wrote site QC: ", site_qc_output)
message("Wrote all-site metadata: ", site_metadata_output)
message("Wrote TensorQTL-compatible raw beta-value BED: ", raw_bed_output)
message("Wrote TensorQTL-compatible inverse-normal BED: ", int_bed_output)

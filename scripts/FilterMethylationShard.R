#!/usr/bin/env Rscript

# Apply chromosome, minimum-coverage, and extreme-coverage QC to one manifest
# shard. The resulting long tables are inputs to MergeMethylationCohort.R.

suppressPackageStartupMessages(library(optparse))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_file)), "MethylationUtils.R"))

option_list <- list(
    make_option("--InputManifest", type = "character", help = "TSV with sample_id and file_path columns [required]"),
    make_option("--OutputPrefix", type = "character", help = "Prefix for output files [required]"),
    make_option("--MinCoverage", type = "double", default = 10,
                help = "Minimum per-call coverage to retain [default: %default]"),
    make_option("--FilterChroms", type = "character", default = "X|Y|M|_",
                help = "Regex for chromosomes/contigs to remove; use '' to keep all [default: %default]"),
    make_option("--FenceK", type = "double", default = 3,
                help = "Tukey log10-coverage far-out fence multiplier [default: %default]"),
    make_option("--AutosomePrefix", type = "character", default = "chr",
                help = "Prefix used for autosome names, e.g. 'chr' for chr1 or '' for 1 [default: %default]"),
    make_option("--NumThreads", type = "integer", default = 1,
                help = "Threads for bgzip decompression and data.table parsing [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$InputManifest) || is.null(opt$OutputPrefix)) stop("--InputManifest and --OutputPrefix are required")
if (!is.finite(opt$MinCoverage) || opt$MinCoverage < 0) stop("--MinCoverage must be a non-negative number")
if (!is.finite(opt$FenceK) || opt$FenceK < 0) stop("--FenceK must be a non-negative number")
if (is.na(opt$NumThreads) || opt$NumThreads < 1) stop("--NumThreads must be at least 1")
setDTthreads(opt$NumThreads)

manifest <- read_manifest(opt$InputManifest)
n_samples <- nrow(manifest)
message("Processing ", n_samples, " sample(s) in this shard")
autosomes <- paste0(opt$AutosomePrefix, seq_len(22))
all_site_calls_by_chrom <- lapply(autosomes, function(...) vector("list", n_samples))
sample_qc_tables <- vector("list", n_samples)
reference_columns <- NULL

for (i in seq_len(n_samples)) {
    sample_id <- manifest$sample_id[i]
    file_path <- manifest$file_path[i]
    message("[", i, "/", n_samples, "] Loading ", sample_id, ": ", file_path)
    methylation_data <- load_methylation_data(file_path, opt$FilterChroms, opt$FenceK, decomp_threads = opt$NumThreads)
    input_columns <- names(methylation_data)
    if (is.null(reference_columns)) {
        reference_columns <- copy(input_columns)
    } else if (!identical(reference_columns, input_columns)) {
        stop("BED columns in sample ", sample_id, " do not match the first input file. All input BED files must use the same schema.")
    }
    if ("sample_id" %in% input_columns) stop("Input BED files must not already contain a 'sample_id' column")
    duplicate_sites <- methylation_data[, .N, by = site_key][N > 1]
    if (nrow(duplicate_sites) > 0) {
        stop("Found ", nrow(duplicate_sites), " duplicated #chrom/begin/end site(s) in ", sample_id,
             ". Aggregate duplicate calls before merging so a site is counted once per sample.")
    }

    coverage_pass <- !is.na(methylation_data$cov) & methylation_data$cov >= opt$MinCoverage
    extreme_pass <- methylation_data$extreme_cov_flag == "ok"
    n_below_min_coverage <- sum(!coverage_pass)
    n_extreme_coverage <- sum(!extreme_pass)
    n_extreme_after_min_coverage <- sum(coverage_pass & !extreme_pass)
    methylation_data[, `:=`(
        sample_id = sample_id,
        meets_min_coverage = coverage_pass,
        per_sample_qc_pass = coverage_pass & extreme_pass
    )]
    n_retained <- sum(methylation_data$per_sample_qc_pass)
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
        n_passing_per_sample_qc = n_retained
    )
    for (chrom_index in seq_along(autosomes)) {
        all_site_calls_by_chrom[[chrom_index]][[i]] <- methylation_data[`#chrom` == autosomes[[chrom_index]]]
    }
    message("  Input sites: ", attr(methylation_data, "n_input_rows"),
            "; removed by chromosome filter: ", attr(methylation_data, "n_input_rows") - nrow(methylation_data),
            "; evaluated for coverage: ", nrow(methylation_data))
    message("  Per-sample thresholds: ", n_below_min_coverage,
            " fail MinCoverage (<", opt$MinCoverage, "); ", n_extreme_after_min_coverage,
            " fail extreme coverage after MinCoverage; ", n_retained, " pass both thresholds")
}

output_dir <- dirname(opt$OutputPrefix)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
sample_qc_output <- paste0(opt$OutputPrefix, ".methylation.sample_qc.tsv")
fwrite(rbindlist(sample_qc_tables, use.names = TRUE), sample_qc_output, sep = "\t", na = "NA")

for (chrom_index in seq_along(autosomes)) {
    chrom_label <- sprintf("autosome%02d", chrom_index)
    chrom_all_calls_output <- paste0(opt$OutputPrefix, ".methylation.", chrom_label, ".per_sample_qc.long.tsv.gz")
    chrom_call_table <- rbindlist(all_site_calls_by_chrom[[chrom_index]], use.names = TRUE)
    fwrite(chrom_call_table, chrom_all_calls_output, sep = "\t", na = "NA")
}
message("Wrote sample QC: ", sample_qc_output)
message("Wrote autosome-split shard call files for ", length(autosomes),
        " chromosome(s); passing calls are derived from per_sample_qc_pass during cohort merge")

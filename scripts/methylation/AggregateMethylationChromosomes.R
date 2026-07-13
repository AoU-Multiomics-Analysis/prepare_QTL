#!/usr/bin/env Rscript

# Finalize the small sample-QC output and global plots after shell streaming has
# concatenated the chromosome-level tables.

suppressPackageStartupMessages(library(optparse))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_file)), "MethylationUtils.R"))

option_list <- list(
    make_option("--SiteMetadata", type = "character", help = "Aggregated site-metadata table [required]"),
    make_option("--SampleQcList", type = "character", help = "One per-shard sample-QC file path per line [required]"),
    make_option("--TotalSamples", type = "integer", help = "Total input sample count [required]"),
    make_option("--OutputPrefix", type = "character", help = "Prefix for output files [required]"),
    make_option("--SampleQcOutput", type = "character", default = NULL,
                help = "Optional sample-QC output path [default: <OutputPrefix>.methylation.sample_qc.tsv]")
)
opt <- parse_args(OptionParser(option_list = option_list))
required_options <- c("SiteMetadata", "SampleQcList", "OutputPrefix")
if (any(vapply(required_options, function(name) is.null(opt[[name]]), logical(1)))) {
    stop("--SiteMetadata, --SampleQcList, and --OutputPrefix are required")
}
if (is.null(opt$TotalSamples) || opt$TotalSamples <= 0) stop("--TotalSamples must be positive")

output_dir <- dirname(opt$OutputPrefix)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

sample_qc <- read_sample_qc(read_file_list(opt$SampleQcList, "Sample-QC"), opt$TotalSamples)
sample_qc_output <- if (is.null(opt$SampleQcOutput)) {
    paste0(opt$OutputPrefix, ".methylation.sample_qc.tsv")
} else {
    opt$SampleQcOutput
}
fwrite(sample_qc, sample_qc_output, sep = "\t", na = "NA")

plot_columns <- c(
    "failure_reason", "has_missing_or_low_coverage", "has_extreme_coverage_loss",
    "pass_sample_presence_filter", "pass_methylation_mad_filter"
)
site_metadata <- fread(opt$SiteMetadata, select = plot_columns)
write_filter_plots(site_metadata, opt$OutputPrefix)

message("Wrote sample QC: ", sample_qc_output)
message("Created global filter plots from ", nrow(site_metadata), " aggregated methylation sites")

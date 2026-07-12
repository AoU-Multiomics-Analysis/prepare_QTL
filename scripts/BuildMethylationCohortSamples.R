#!/usr/bin/env Rscript

# Build a deterministic cohort sample list from one or more per-sample QC files.

suppressPackageStartupMessages(library(optparse))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_file)), "MethylationUtils.R"))

option_list <- list(
    make_option("--SampleQcList", type = "character", help = "One per-sample or per-shard sample-QC file path per line [required]"),
    make_option("--OutputPrefix", type = "character", help = "Prefix for output files [required]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$SampleQcList) || is.null(opt$OutputPrefix)) {
    stop("--SampleQcList and --OutputPrefix are required")
}

sample_qc_paths <- read_file_list(opt$SampleQcList, "Sample-QC")
sample_qc <- rbindlist(lapply(sample_qc_paths, fread), use.names = TRUE)
if (!("sample_id" %in% names(sample_qc))) stop("Sample-QC files must contain a sample_id column")
sample_qc[, sample_id := as.character(sample_id)]
if (nrow(sample_qc) == 0) stop("Sample-QC files contain no samples")
if (anyNA(sample_qc$sample_id) || any(!nzchar(sample_qc$sample_id))) stop("Sample-QC files contain an empty sample_id")
if (anyDuplicated(sample_qc$sample_id)) stop("A sample appears in more than one Sample-QC file")

output_dir <- dirname(opt$OutputPrefix)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
cohort_samples_output <- paste0(opt$OutputPrefix, ".methylation.cohort_samples.tsv")
cohort_sample_qc_output <- paste0(opt$OutputPrefix, ".methylation.cohort_sample_qc.tsv")
total_samples_output <- paste0(opt$OutputPrefix, ".methylation.total_samples.txt")
fwrite(sample_qc[, .(sample_id)], cohort_samples_output, sep = "\t")
fwrite(sample_qc, cohort_sample_qc_output, sep = "\t", na = "NA")
writeLines(as.character(nrow(sample_qc)), total_samples_output)

message("Built cohort sample list with ", nrow(sample_qc), " sample(s)")

#!/usr/bin/env Rscript

# Concatenate chromosome-level methylation cohort outputs into the final files.

suppressPackageStartupMessages(library(optparse))
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
source(file.path(dirname(normalizePath(script_file)), "MethylationUtils.R"))

option_list <- list(
    make_option("--FilteredCallList", type = "character", help = "One chromosome filtered-call file path per line [required]"),
    make_option("--SiteQcList", type = "character", help = "One chromosome site-QC file path per line [required]"),
    make_option("--SiteMetadataList", type = "character", help = "One chromosome site-metadata file path per line [required]"),
    make_option("--RawBedList", type = "character", help = "One chromosome raw methylation BED path per line [required]"),
    make_option("--IntBedList", type = "character", help = "One chromosome INT methylation BED path per line [required]"),
    make_option("--SampleQcList", type = "character", help = "One per-shard sample-QC file path per line [required]"),
    make_option("--TotalSamples", type = "integer", help = "Total input sample count [required]"),
    make_option("--OutputPrefix", type = "character", help = "Prefix for output files [required]")
)
opt <- parse_args(OptionParser(option_list = option_list))
required_options <- c("FilteredCallList", "SiteQcList", "SiteMetadataList", "RawBedList", "IntBedList", "SampleQcList", "OutputPrefix")
if (any(vapply(required_options, function(name) is.null(opt[[name]]), logical(1)))) {
    stop("--FilteredCallList, --SiteQcList, --SiteMetadataList, --RawBedList, --IntBedList, --SampleQcList, and --OutputPrefix are required")
}
if (is.null(opt$TotalSamples) || opt$TotalSamples <= 0) stop("--TotalSamples must be positive")

read_and_bind <- function(list_path, label) {
    paths <- read_file_list(list_path, label)
    rbindlist(lapply(paths, fread), use.names = TRUE)
}

output_dir <- dirname(opt$OutputPrefix)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

filtered_calls <- read_and_bind(opt$FilteredCallList, "Chromosome filtered-call")
site_qc <- read_and_bind(opt$SiteQcList, "Chromosome site-QC")
site_metadata <- read_and_bind(opt$SiteMetadataList, "Chromosome site-metadata")
raw_bed <- read_and_bind(opt$RawBedList, "Chromosome raw BED")
int_bed <- read_and_bind(opt$IntBedList, "Chromosome INT BED")
sample_qc <- read_sample_qc(read_file_list(opt$SampleQcList, "Sample-QC"), opt$TotalSamples)

if (nrow(filtered_calls) > 0) setorder(filtered_calls, `#chrom`, begin, end, sample_id)
if (nrow(site_qc) > 0) setorder(site_qc, `#chrom`, begin, end)
if (nrow(site_metadata) > 0) setorder(site_metadata, `#chrom`, begin, end)
if (nrow(raw_bed) > 0) setorder(raw_bed, `#chr`, start, end)
if (nrow(int_bed) > 0) setorder(int_bed, `#chr`, start, end)

long_output <- paste0(opt$OutputPrefix, ".methylation.filtered.long.tsv.gz")
site_qc_output <- paste0(opt$OutputPrefix, ".methylation.site_qc.tsv.gz")
site_metadata_output <- paste0(opt$OutputPrefix, ".methylation.site_metadata.tsv.gz")
sample_qc_output <- paste0(opt$OutputPrefix, ".methylation.sample_qc.tsv")
raw_bed_output <- paste0(opt$OutputPrefix, ".methylation.raw.bed.gz")
int_bed_output <- paste0(opt$OutputPrefix, ".methylation.INT.bed.gz")

fwrite(filtered_calls, long_output, sep = "\t", na = "NA")
fwrite(site_qc, site_qc_output, sep = "\t", na = "NA")
fwrite(site_metadata, site_metadata_output, sep = "\t", na = "NA")
fwrite(sample_qc, sample_qc_output, sep = "\t", na = "NA")
fwrite(raw_bed, raw_bed_output, sep = "\t", na = "NA")
fwrite(int_bed, int_bed_output, sep = "\t", na = "NA")
write_filter_plots(site_metadata, opt$OutputPrefix)

message("Aggregated ", nrow(site_qc), " methylation sites across chromosome merge outputs")
message("Wrote filtered long calls: ", long_output)
message("Wrote site QC: ", site_qc_output)
message("Wrote all-site metadata: ", site_metadata_output)
message("Wrote sample QC: ", sample_qc_output)
message("Wrote TensorQTL-compatible raw beta-value BED: ", raw_bed_output)
message("Wrote TensorQTL-compatible inverse-normal BED: ", int_bed_output)

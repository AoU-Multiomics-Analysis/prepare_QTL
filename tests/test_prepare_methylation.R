#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
repo_root <- normalizePath(file.path(dirname(file_arg), ".."), mustWork = TRUE)
script_path <- Sys.getenv(
    "PREPARE_METHYLATION_SCRIPT",
    unset = file.path(repo_root, "scripts", "methylation", "PrepareMethylation.R")
)

assert_true <- function(condition, message) {
    if (!isTRUE(condition)) stop(message, call. = FALSE)
}

assert_equal <- function(actual, expected, message, tolerance = 1e-10) {
    same <- if (is.numeric(actual) && is.numeric(expected)) {
        length(actual) == length(expected) && all(abs(actual - expected) <= tolerance)
    } else {
        identical(actual, expected)
    }
    if (!same) stop(message, call. = FALSE)
}

work_dir <- tempfile("prepare-methylation-test-")
dir.create(work_dir, recursive = TRUE)
on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

set.seed(1729)
sample_ids <- sprintf("sample_%02d", seq_len(40))
sample_columns <- paste0(sample_ids, ".meth_region_stats")
base_signal <- seq(0.05, 0.95, length.out = 47)
values <- vapply(
    seq_len(39),
    function(index) pmin(1, pmax(0, base_signal + rnorm(length(base_signal), sd = 0.002))),
    numeric(length(base_signal))
)
values <- cbind(values, rev(base_signal))
colnames(values) <- sample_columns

input <- data.frame(
    chrom = c(rep("chr1", 40), "chr2", "chr3", "chrX", "chr1_random", "chrY", "chrUn", "chrM"),
    start = seq(100L, by = 10L, length.out = 47),
    end = seq(101L, by = 10L, length.out = 47),
    region = paste0("region_", seq_len(47)),
    values,
    check.names = FALSE
)
input[41, sample_columns[[1]]] <- NA_real_
input[42, sample_columns[1:2]] <- NA_real_

input_path <- file.path(work_dir, "merged_methylation.bed")
sample_list_path <- file.path(work_dir, "samples.tsv")
output_prefix <- file.path(work_dir, "cohort")
write.table(input, input_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(data.frame(sample_id = sample_ids), sample_list_path, sep = "\t", quote = FALSE, row.names = FALSE)

run_prepare <- function(bed_path, samples_path, prefix, missingness_threshold = NULL) {
    command_args <- c(
        shQuote(script_path),
        "--MethylationBed", shQuote(bed_path),
        "--SampleList", shQuote(samples_path),
        "--OutputPrefix", shQuote(prefix)
    )
    if (!is.null(missingness_threshold)) {
        command_args <- c(command_args, "--MissingnessThreshold", as.character(missingness_threshold))
    }
    output <- suppressWarnings(system2("Rscript", command_args, stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    list(status = status, output = output)
}

command_result <- run_prepare(input_path, sample_list_path, output_prefix)
assert_equal(command_result$status, 0L,
             paste(c("PrepareMethylation.R failed:", command_result$output), collapse = "\n"))

read_tsv <- function(path) {
    read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

raw_path <- paste0(output_prefix, ".methylation.raw.bed.gz")
int_path <- paste0(output_prefix, ".methylation.INT.bed.gz")
scaled_path <- paste0(output_prefix, ".methylation.scaled.bed.gz")
int_outlier_path <- paste0(output_prefix, ".methylation.INT.connectivity_outliers.tsv")
scaled_outlier_path <- paste0(output_prefix, ".methylation.scaled.connectivity_outliers.tsv")

expected_paths <- c(raw_path, int_path, scaled_path, int_outlier_path, scaled_outlier_path)
assert_true(all(file.exists(expected_paths)), "One or more expected output files were not created")

raw <- read_tsv(raw_path)
int <- read_tsv(int_path)
scaled <- read_tsv(scaled_path)
int_outliers <- read_tsv(int_outlier_path)
scaled_outliers <- read_tsv(scaled_outlier_path)

assert_equal(names(raw)[1:4], c("#chr", "start", "end", "phenotype_id"), "The BED metadata columns are not TensorQTL-compatible")
assert_equal(nrow(raw), 41L, "The feature missingness or chromosome filter retained the wrong number of rows")
assert_true("region_41" %in% raw$phenotype_id, "A feature below 5% missingness was removed")
assert_true(!("region_42" %in% raw$phenotype_id), "A feature above 5% missingness was retained")
assert_true(!any(raw$`#chr` %in% c("chrX", "chr1_random", "chrY", "chrUn", "chrM")), "A disallowed chromosome was retained")
assert_equal(names(raw)[5:44], sample_ids, "Sample suffixes were not removed or sample order changed")
assert_true("sample_40" %in% names(raw), "The raw BED must keep connectivity outliers")

region_41 <- raw[raw$phenotype_id == "region_41", , drop = FALSE]
expected_imputation <- mean(as.numeric(input[41, sample_columns[-1]]), na.rm = TRUE)
assert_equal(region_41$sample_01, expected_imputation, "Missing methylation was not imputed with the feature mean")

assert_equal(int_outliers$SampleID, "sample_40", "The INT branch did not report the disconnected sample")
assert_equal(scaled_outliers$SampleID, "sample_40", "The scaled branch did not report the disconnected sample")
assert_true(!("sample_40" %in% names(int)), "The INT BED retained a connectivity outlier")
assert_true(!("sample_40" %in% names(scaled)), "The scaled BED retained a connectivity outlier")
assert_equal(names(int)[5:ncol(int)], sample_ids[-40], "The INT BED removed an unexpected sample")
assert_equal(names(scaled)[5:ncol(scaled)], sample_ids[-40], "The scaled BED removed an unexpected sample")

scaled_values <- as.numeric(scaled[1, 5:ncol(scaled)])
assert_true(abs(mean(scaled_values)) < 0.2, "The scaled BED does not contain centered feature values")
assert_true(all(is.finite(as.matrix(int[, 5:ncol(int)]))), "The INT BED contains non-finite values")

headerless_sample_list <- file.path(work_dir, "samples_headerless.tsv")
headerless_prefix <- file.path(work_dir, "headerless")
write.table(sample_ids, headerless_sample_list, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
headerless_result <- run_prepare(input_path, headerless_sample_list, headerless_prefix, missingness_threshold = 6)
assert_equal(headerless_result$status, 0L,
             paste(c("Headerless sample-list run failed:", headerless_result$output), collapse = "\n"))
headerless_raw <- read_tsv(paste0(headerless_prefix, ".methylation.raw.bed.gz"))
assert_equal(names(headerless_raw)[5:44], sample_ids, "A headerless sample list lost or reordered a sample")
assert_true("region_42" %in% headerless_raw$phenotype_id, "The configurable missingness threshold was not applied")

missing_sample_list <- file.path(work_dir, "samples_with_missing_id.tsv")
write.table(data.frame(sample_id = c(sample_ids, "sample_not_in_bed")), missing_sample_list,
            sep = "\t", quote = FALSE, row.names = FALSE)
missing_sample_result <- run_prepare(input_path, missing_sample_list, file.path(work_dir, "missing_sample"))
assert_true(missing_sample_result$status != 0L, "A requested sample absent from the BED did not stop preparation")
assert_true(any(grepl("sample_not_in_bed", missing_sample_result$output, fixed = TRUE)),
            "The absent-sample error does not identify the missing sample")

all_filtered_path <- file.path(work_dir, "all_filtered.bed")
write.table(input[43:47, ], all_filtered_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
all_filtered_result <- run_prepare(all_filtered_path, sample_list_path, file.path(work_dir, "all_filtered"))
assert_true(all_filtered_result$status != 0L, "An input with no retained features did not stop preparation")
assert_true(any(grepl("At least two methylation features", all_filtered_result$output, fixed = TRUE)),
            "The no-feature error does not explain the minimum feature requirement")

missing_chromosome_input <- input
missing_chromosome_input$chrom[[1]] <- NA_character_
missing_chromosome_path <- file.path(work_dir, "missing_chromosome.bed")
write.table(missing_chromosome_input, missing_chromosome_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
missing_chromosome_result <- run_prepare(missing_chromosome_path, sample_list_path,
                                         file.path(work_dir, "missing_chromosome"))
assert_true(missing_chromosome_result$status != 0L, "A missing chromosome value did not stop preparation")
assert_true(any(grepl("missing interval or phenotype metadata", missing_chromosome_result$output, fixed = TRUE)),
            "The missing-metadata error is not specific")

message("PrepareMethylation integration test passed")

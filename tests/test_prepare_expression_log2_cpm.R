#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
repo_root <- normalizePath(file.path(dirname(file_arg), ".."), mustWork = TRUE)
script_path <- Sys.getenv(
    "PREPARE_EXPRESSION_SCRIPT",
    unset = file.path(repo_root, "scripts", "expression", "PrepareExpression.R")
)

assert_true <- function(condition, message) {
    if (!isTRUE(condition)) stop(message, call. = FALSE)
}

assert_equal <- function(actual, expected, message, tolerance = 1e-7) {
    same <- if (is.numeric(actual) && is.numeric(expected)) {
        length(actual) == length(expected) && all(abs(actual - expected) <= tolerance)
    } else {
        identical(actual, expected)
    }
    if (!same) stop(message, call. = FALSE)
}

work_dir <- tempfile("prepare-expression-log2-cpm-test-")
dir.create(work_dir, recursive = TRUE)
on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

input <- data.frame(
    "#chr" = "chr2",
    start = 20L,
    end = 21L,
    gene_id = "gene_b",
    S3 = 3,
    S1 = 0,
    S2 = 1,
    check.names = FALSE
)
input_path <- file.path(work_dir, "cell_type_log2_cpm.bed.gz")
sample_list_path <- file.path(work_dir, "samples.tsv")
output_prefix <- file.path(work_dir, "cohort")

data.table::fwrite(input, input_path, sep = "\t")
write.table(
    c("S2", "S3", "S1"),
    sample_list_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)

command_args <- c(
    shQuote(script_path),
    "--Log2CpmBed", shQuote(input_path),
    "--SampleList", shQuote(sample_list_path),
    "--OutputPrefix", shQuote(output_prefix)
)
command_output <- suppressWarnings(system2(
    "Rscript",
    command_args,
    stdout = TRUE,
    stderr = TRUE
))
command_status <- attr(command_output, "status")
if (is.null(command_status)) command_status <- 0L
assert_equal(
    command_status,
    0L,
    paste(c("PrepareExpression.R failed:", command_output), collapse = "\n")
)

read_bed <- function(suffix) {
    read.delim(
        paste0(output_prefix, suffix),
        check.names = FALSE,
        stringsAsFactors = FALSE
    )
}

raw <- read_bed(".expression.raw.bed.gz")
int <- read_bed(".expression.INT.bed.gz")
scaled <- read_bed(".expression.scaled.bed.gz")

expected_columns <- c("#chr", "start", "end", "gene_id", "S2", "S3", "S1")
assert_equal(names(raw), expected_columns, "The raw BED did not preserve metadata or sample-list order")
assert_equal(raw[1, 1:4], input[1, 1:4], "The raw BED changed input coordinates or the gene identifier")
assert_equal(as.numeric(raw[1, 5:7]), c(1, 3, 0), "The raw BED changed supplied log2 CPM values")
assert_equal(
    as.numeric(int[1, 5:7]),
    c(0, 0.8694238, -0.8694238),
    "The INT BED did not rank-normalize the supplied log2 CPM values"
)
assert_equal(
    as.numeric(scaled[1, 5:7]),
    c(-0.2182179, 1.0910895, -0.8728716),
    "The scaled BED did not scale log2 CPM directly"
)
assert_true(
    file.exists(paste0(output_prefix, ".expression.INT.connectivity_outliers.tsv")),
    "The INT connectivity-outlier report was not written"
)
assert_true(
    file.exists(paste0(output_prefix, ".expression.scaled.connectivity_outliers.tsv")),
    "The scaled connectivity-outlier report was not written"
)

run_invalid <- function(input_args, prefix_name) {
    invalid_args <- c(
        shQuote(script_path),
        input_args,
        "--SampleList", shQuote(sample_list_path),
        "--OutputPrefix", shQuote(file.path(work_dir, prefix_name))
    )
    output <- suppressWarnings(system2(
        "Rscript",
        invalid_args,
        stdout = TRUE,
        stderr = TRUE
    ))
    status <- attr(output, "status")
    if (is.null(status)) status <- 0L
    list(status = status, output = output)
}

neither_mode <- run_invalid(character(), "neither")
assert_true(neither_mode$status != 0L, "Missing expression inputs did not stop preparation")
assert_true(
    any(grepl("Provide exactly one of --CountGCT or --Log2CpmBed", neither_mode$output, fixed = TRUE)),
    "The missing-input error did not explain the mutually exclusive modes"
)

both_modes <- run_invalid(c(
    "--CountGCT", shQuote(input_path),
    "--AnnotationGTF", shQuote(input_path),
    "--Log2CpmBed", shQuote(input_path)
), "both")
assert_true(both_modes$status != 0L, "Supplying both expression inputs did not stop preparation")
assert_true(
    any(grepl("Provide exactly one of --CountGCT or --Log2CpmBed", both_modes$output, fixed = TRUE)),
    "The multiple-input error did not explain the mutually exclusive modes"
)

missing_annotation <- run_invalid(c(
    "--CountGCT", shQuote(input_path)
), "missing-annotation")
assert_true(missing_annotation$status != 0L, "Raw counts without a GTF did not stop preparation")
assert_true(
    any(grepl("--AnnotationGTF is required with --CountGCT", missing_annotation$output, fixed = TRUE)),
    "The raw-count error did not identify the required GTF"
)

message("PrepareExpression log2 CPM BED integration test passed")

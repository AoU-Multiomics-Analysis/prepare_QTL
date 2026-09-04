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
    "#chr" = c("chr2", "chr6"),
    start = c(20L, 30L),
    end = c(21L, 31L),
    gene_id = c("00123", "HLA-DRA"),
    "S-1" = c(3, 4),
    "X100" = c(0, 1),
    "1002" = c(1, 2),
    check.names = FALSE
)
input_path <- file.path(work_dir, "cell_type_log2_cpm.bed.gz")
sample_list_path <- file.path(work_dir, "samples.tsv")
output_prefix <- file.path(work_dir, "cohort")

data.table::fwrite(input, input_path, sep = "\t")
write.table(
    c("1002", "S-1", "X100"),
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

expected_columns <- c("#chr", "start", "end", "gene_id", "1002", "S-1", "X100")
for (bed_output in list(raw = raw, int = int, scaled = scaled)) {
    assert_equal(
        names(bed_output),
        expected_columns,
        "An expression BED did not preserve metadata or sample-list order"
    )
    assert_equal(
        bed_output[, 1:4],
        input[, 1:4],
        "An expression BED changed coordinates, gene identifiers, or gene order"
    )
}
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

# A missing or double log transform must change these expected scaled values.
# Linear values 1, 7, 0 become 1, 3, 0 after log2(CPM + 1).
cpm_input <- input |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c("S-1", "X100", "1002")), ~ 2^.x - 1))
cpm_input_path <- file.path(work_dir, "linear CPM.bed.gz")
cpm_prefix <- file.path(work_dir, "linear-cpm-cohort")
data.table::fwrite(cpm_input, cpm_input_path, sep = "\t")
cpm_output <- suppressWarnings(system2(
    "Rscript",
    c(shQuote(script_path), "--CpmBed", shQuote(cpm_input_path),
      "--SampleList", shQuote(sample_list_path), "--OutputPrefix", shQuote(cpm_prefix)),
    stdout = TRUE, stderr = TRUE
))
cpm_status <- attr(cpm_output, "status")
if (is.null(cpm_status)) cpm_status <- 0L
assert_equal(cpm_status, 0L, paste(c("CPM BED preparation failed:", cpm_output), collapse = "\n"))
for (suffix in c("raw", "INT", "scaled")) {
    cpm_bed <- read.delim(paste0(cpm_prefix, ".expression.", suffix, ".bed.gz"),
                         check.names = FALSE, colClasses = c(gene_id = "character"))
    assert_equal(names(cpm_bed), expected_columns, "CPM BED changed the sample order")
    assert_equal(cpm_bed[, 1:4], input[, 1:4], "CPM BED changed coordinates or gene IDs")
    expected_values <- switch(suffix,
        raw = c(1, 7, 0),
        INT = c(0, 0.8694238, -0.8694238),
        scaled = c(-0.2182179, 1.0910895, -0.8728716))
    assert_equal(as.numeric(cpm_bed[1, 5:7]), expected_values,
                 paste("CPM BED has the wrong transformation for", suffix))
}
for (suffix in c("INT", "scaled")) {
    assert_true(file.exists(paste0(cpm_prefix, ".expression.", suffix, ".connectivity_outliers.tsv")),
                paste("CPM BED is missing its", suffix, "outlier report"))
}

numeric_gene_input <- data.frame(
    "#chr" = c("chr2", "chr6"),
    start = c(20L, 30L),
    end = c(21L, 31L),
    gene_id = c("00123", "00456"),
    "1001" = c(0, 1),
    "1002" = c(1, 2),
    "1003" = c(3, 4),
    check.names = FALSE
)
numeric_gene_input_path <- file.path(work_dir, "numeric_gene_ids.bed.gz")
numeric_sample_list_path <- file.path(work_dir, "numeric_samples.tsv")
numeric_gene_prefix <- file.path(work_dir, "numeric-gene-cohort")
data.table::fwrite(numeric_gene_input, numeric_gene_input_path, sep = "\t", quote = FALSE)
write.table(
    c("1002", "1003", "1001"),
    numeric_sample_list_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
)
numeric_gene_output <- suppressWarnings(system2(
    "Rscript",
    c(
        shQuote(script_path),
        "--Log2CpmBed", shQuote(numeric_gene_input_path),
        "--SampleList", shQuote(numeric_sample_list_path),
        "--OutputPrefix", shQuote(numeric_gene_prefix)
    ),
    stdout = TRUE,
    stderr = TRUE
))
numeric_gene_status <- attr(numeric_gene_output, "status")
if (is.null(numeric_gene_status)) numeric_gene_status <- 0L
assert_equal(
    numeric_gene_status,
    0L,
    paste(c("PrepareExpression.R failed for numeric gene IDs:", numeric_gene_output), collapse = "\n")
)
numeric_gene_raw <- read.delim(
    paste0(numeric_gene_prefix, ".expression.raw.bed.gz"),
    check.names = FALSE,
    colClasses = c(gene_id = "character"),
    stringsAsFactors = FALSE
)
assert_equal(
    numeric_gene_raw$gene_id,
    numeric_gene_input$gene_id,
    "The raw BED did not preserve leading zeros in numeric gene identifiers"
)
assert_equal(
    names(numeric_gene_raw)[5:7],
    c("1002", "1003", "1001"),
    "The raw BED did not preserve numeric sample IDs or sample-list order"
)

count_gct_path <- file.path(work_dir, "counts.gct")
annotation_gtf_path <- file.path(work_dir, "genes.gtf")
raw_count_prefix <- file.path(work_dir, "raw-count-cohort")
writeLines(c(
    "#1.2",
    "2\t3",
    "Name\tDescription\tX100\t1002\tS-1",
    "00123\tgene one\t10\t20\t30",
    "HLA-DRA\tgene two\t40\t20\t10"
), count_gct_path)
writeLines(c(
    "chr2\ttest\tgene\t21\t30\t.\t+\t.\tgene_id \"00123\"; gene_name \"G1\";",
    "chr6\ttest\tgene\t31\t40\t.\t-\t.\tgene_id \"HLA-DRA\"; gene_name \"G2\";"
), annotation_gtf_path)
raw_count_output <- suppressWarnings(system2(
    "Rscript",
    c(
        shQuote(script_path),
        "--CountGCT", shQuote(count_gct_path),
        "--AnnotationGTF", shQuote(annotation_gtf_path),
        "--SampleList", shQuote(sample_list_path),
        "--OutputPrefix", shQuote(raw_count_prefix)
    ),
    stdout = TRUE,
    stderr = TRUE
))
raw_count_status <- attr(raw_count_output, "status")
if (is.null(raw_count_status)) raw_count_status <- 0L
assert_equal(
    raw_count_status,
    0L,
    paste(c("PrepareExpression.R raw-count mode failed:", raw_count_output), collapse = "\n")
)
raw_count_bed <- read.delim(
    paste0(raw_count_prefix, ".expression.raw.bed.gz"),
    check.names = FALSE,
    stringsAsFactors = FALSE
)
assert_equal(
    names(raw_count_bed),
    expected_columns,
    "Raw-count mode changed its BED metadata or sample-list order"
)
assert_true(
    all(is.finite(as.matrix(raw_count_bed[, 5:7]))),
    "Raw-count mode wrote non-finite CPM values"
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
    any(grepl("Provide exactly one of --CountGCT, --CpmBed, or --Log2CpmBed", neither_mode$output, fixed = TRUE)),
    "The missing-input error did not explain the mutually exclusive modes"
)

both_modes <- run_invalid(c(
    "--CountGCT", shQuote(input_path),
    "--AnnotationGTF", shQuote(input_path),
    "--Log2CpmBed", shQuote(input_path)
), "both")
assert_true(both_modes$status != 0L, "Supplying both expression inputs did not stop preparation")
assert_true(
    any(grepl("Provide exactly one of --CountGCT, --CpmBed, or --Log2CpmBed", both_modes$output, fixed = TRUE)),
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

log2_with_annotation <- run_invalid(c(
    "--Log2CpmBed", shQuote(input_path),
    "--AnnotationGTF", shQuote(input_path)
), "log2-with-annotation")
assert_true(
    log2_with_annotation$status != 0L,
    "Supplying a GTF with a log2 CPM BED did not stop preparation"
)
assert_true(
    any(grepl("--AnnotationGTF cannot be used with --Log2CpmBed", log2_with_annotation$output, fixed = TRUE)),
    "The log2 CPM mode error did not reject the GTF input"
)

constant_gene_input <- input
constant_gene_input[1, c("S-1", "X100", "1002")] <- 2
constant_gene_input_path <- file.path(work_dir, "constant_gene.bed.gz")
data.table::fwrite(constant_gene_input, constant_gene_input_path, sep = "\t", quote = FALSE)
constant_gene <- run_invalid(c(
    "--Log2CpmBed", shQuote(constant_gene_input_path)
), "constant-gene")
assert_true(
    constant_gene$status != 0L,
    "A zero-variance gene in the log2 CPM BED did not stop preparation"
)
assert_true(
    any(grepl("Log2 CPM BED contains genes with zero variance", constant_gene$output, fixed = TRUE)),
    "The zero-variance error did not identify the invalid log2 CPM input"
)

for (other in c("--CountGCT", "--Log2CpmBed")) {
    conflicting <- run_invalid(c("--CpmBed", shQuote(cpm_input_path), other, shQuote(input_path)),
                               paste0("conflict-", sub("--", "", other)))
    assert_true(conflicting$status != 0L &&
                    any(grepl("Provide exactly one", conflicting$output, fixed = TRUE)),
                paste("CPM BED was not rejected together with", other))
}
cpm_with_gtf <- run_invalid(c("--CpmBed", shQuote(cpm_input_path),
                              "--AnnotationGTF", shQuote(annotation_gtf_path)), "cpm-with-gtf")
assert_true(cpm_with_gtf$status != 0L &&
                any(grepl("--AnnotationGTF cannot be used with --CpmBed", cpm_with_gtf$output, fixed = TRUE)),
            "CPM BED incorrectly accepted a GTF")

# Validate CPM before any outputs, but keep valid negative log2 values unchanged.
for (value in c(-0.1, -2, Inf, NA_real_)) {
    invalid_cpm <- cpm_input
    invalid_cpm[["S-1"]][[1L]] <- value
    path <- file.path(work_dir, paste0("invalid-", value, ".bed.gz"))
    data.table::fwrite(invalid_cpm, path, sep = "\t")
    prefix <- paste0("invalid-cpm-", value)
    invalid <- run_invalid(c("--CpmBed", shQuote(path)), prefix)
    assert_true(invalid$status != 0L, "Invalid CPM values did not stop preparation")
    expected_error <- if (is.finite(value)) "negative" else "finite numeric"
    assert_true(any(grepl(expected_error, invalid$output, fixed = TRUE)),
                paste("CPM validation did not explain the invalid value:", value))
    assert_true(!file.exists(file.path(work_dir, paste0(prefix, ".expression.INT.bed.gz"))),
                "Invalid CPM values produced partial expression outputs")
}
negative_log <- input |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c("S-1", "X100", "1002")), ~ .x - 5))
negative_log_path <- file.path(work_dir, "negative-log.bed.gz")
data.table::fwrite(negative_log, negative_log_path, sep = "\t")
negative_result <- run_invalid(c("--Log2CpmBed", shQuote(negative_log_path)), "negative-log")
assert_equal(negative_result$status, 0L, "Valid negative log2 CPM values were rejected")
negative_scaled <- read.delim(file.path(work_dir, "negative-log.expression.scaled.bed.gz"),
                             check.names = FALSE)
assert_equal(as.numeric(negative_scaled[1, 5:7]), c(-0.2182179, 1.0910895, -0.8728716),
             "Negative log2 CPM values were clipped or transformed again")
cpm_constant <- run_invalid(c("--CpmBed", shQuote(constant_gene_input_path)), "constant-cpm")
assert_true(cpm_constant$status != 0L && any(grepl("zero variance", cpm_constant$output, fixed = TRUE)),
            "Constant CPM genes did not stop preparation")

message("PrepareExpression CPM and log2 CPM BED integration tests passed")

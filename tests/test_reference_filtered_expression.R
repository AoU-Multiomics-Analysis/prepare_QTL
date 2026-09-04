#!/usr/bin/env Rscript

# Real local integration: per-cell filtering -> unchanged CPM -> eQTL transforms.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
path <- gsub("~+~", " ", sub("^--file=", "", file_arg), fixed = TRUE)
root <- normalizePath(file.path(dirname(path), ".."))
script_root <- file.path(root, "scripts/cell_type_specific_expression")
Sys.setenv(CELL_TYPE_SPECIFIC_EXPRESSION_ROOT = script_root)
source(file.path(script_root, "bootstrap.R"))
directory <- tempfile("filtered-expression-integration-")
dir.create(directory)

check <- function(condition, label) {
  if (!isTRUE(condition)) stop(label, call. = FALSE)
}

coordinates <- tibble::tibble(`#chr` = "chr1", start = 1:4, end = 2:5,
  gene_id = paste0("g", 1:4))
values <- matrix(c(0.12345678901234567, 1, 7, -1, 1, 3, 0, 3, 15, 2, 3, 6),
  nrow = 4, byrow = TRUE, dimnames = list(coordinates$gene_id, c("S1", "S2", "S3")))
source_beds <- write_cell_type_beds(list("CD4 T cells" = values), coordinates,
  file.path(directory, "source"))
reference <- tibble::tibble(gene_id = paste0("g", 1:4), cell_type = "CD4 T cells",
  n_samples = 2L, mean_log2_cpm1 = c(4, 4, 0.001, 4), median_log2_cpm1 = c(4, 4, 0.001, 4))
filter_directory <- file.path(directory, "filtered")
filter_cell_type_beds(source_beds$inventory, unname(source_beds$paths), filter_directory,
  reference_summary = reference, chunk_size = 1L)
filtered_path <- readLines(file.path(filter_directory, "filtered_beds.txt"))
filtered <- readr::read_tsv(filtered_path, show_col_types = FALSE)
check(identical(filtered$gene_id, c("g1", "g4")), "Wrong negative/reference-filtered genes")
check(identical(as.matrix(filtered[-(1:4)]), unname(values[c(1, 4), , drop = FALSE])) ||
  isTRUE(all.equal(unname(as.matrix(filtered[-(1:4)])),
    unname(values[c(1, 4), , drop = FALSE]), tolerance = 0)), "CPM values changed during filtering")

sample_path <- file.path(directory, "samples.txt")
writeLines(c("S3", "S1", "S2"), sample_path)
output_prefix <- file.path(directory, "qtl")
prepare_script <- Sys.getenv("PREPARE_EXPRESSION_SCRIPT",
  unset = file.path(root, "scripts/expression/PrepareExpression.R"))
logs <- suppressWarnings(system2("Rscript", shQuote(c(prepare_script,
  "--CpmBed", filtered_path, "--SampleList", sample_path, "--OutputPrefix", output_prefix)),
  stdout = TRUE, stderr = TRUE))
check(is.null(attr(logs, "status")), paste(logs, collapse = "\n"))
raw <- readr::read_tsv(paste0(output_prefix, ".expression.raw.bed.gz"), show_col_types = FALSE)
scaled <- readr::read_tsv(paste0(output_prefix, ".expression.scaled.bed.gz"), show_col_types = FALSE)
int <- readr::read_tsv(paste0(output_prefix, ".expression.INT.bed.gz"), show_col_types = FALSE)
check(identical(raw$gene_id, c("g1", "g4")) && identical(scaled$gene_id, raw$gene_id) &&
  identical(int$gene_id, raw$gene_id), "QTL preparation changed filtered gene membership")
check(identical(names(raw)[-(1:4)], c("S3", "S1", "S2")), "Sample-list ordering changed")
expected_cpm <- values[c(1, 4), c("S3", "S1", "S2"), drop = FALSE]
check(isTRUE(all.equal(unname(as.matrix(raw[-(1:4)])), unname(expected_cpm), tolerance = 1e-12)),
  "QTL raw output recalculated CPM")
expected_scaled <- t(apply(log2(expected_cpm + 1), 1, function(x) as.numeric(scale(x))))
check(isTRUE(all.equal(unname(as.matrix(scaled[-(1:4)])), unname(expected_scaled), tolerance = 1e-7)),
  "QTL scaled output must apply log2(CPM + 1) exactly once")
check(isTRUE(all.equal(as.numeric(int[1, -(1:4)]), c(0.8694238, -0.8694238, 0), tolerance = 1e-7)),
  "INT ranks did not match retained CPM values")
message("Reference-filtered CPM -> PrepareExpression integration passed")

#!/usr/bin/env Rscript

# Synthetic counts for a plumbing test, not a biological validation reference.
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: generate_reference_fixture.R EXPRESSION_BED OUTPUT_COUNTS_TSV_GZ", call. = FALSE)
}
bed <- readr::read_tsv(arguments[[1]], show_col_types = FALSE, progress = FALSE)
populations <- c("NveB", "MemB", "CD4T", "CD8T", "NK", "Mono", "MonoNonClassical",
  "Neut", "Eo", "myDC", "myDC123", "pDC")
samples <- unlist(purrr::map(populations, function(population) {
  paste(population, seq_len(if (population == "NveB") 3L else 2L), sep = ".")
}), use.names = FALSE)
counts <- tibble::tibble(gene_id = bed$gene_id)
gene_index <- seq_len(nrow(bed))
for (sample_index in seq_along(samples)) {
  counts[[samples[[sample_index]]]] <- as.integer(
    100 + gene_index * 7 + (gene_index * sample_index * 13) %% 179
  )
}
dir.create(dirname(arguments[[2]]), recursive = TRUE, showWarnings = FALSE)
readr::write_tsv(counts, arguments[[2]])
message(sprintf("Reference smoke fixture: genes=%d samples=%d", nrow(counts), length(samples)))

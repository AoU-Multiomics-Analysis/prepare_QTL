#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L || !nzchar(arguments[[1L]])) {
  stop(
    "Usage: generate_synthetic_fixture.R OUTPUT_DIRECTORY",
    call. = FALSE
  )
}

output_directory <- arguments[[1L]]
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

seed <- 20260901L
set.seed(seed)

# Eight digits after the decimal are more than 1,000 times coarser than the
# observed cross-BLAS accumulation noise of approximately 7.3e-12. This
# precision also keeps rounded 22-part proportion row sums within 1e-6.
fixture_decimal_digits <- 8L
format_fixture_numeric <- function(values) {
  formatC(
    values,
    format = "f",
    digits = fixture_decimal_digits,
    decimal.mark = "."
  )
}
stabilize_fixture_table <- function(table) {
  table |>
    dplyr::mutate(
      dplyr::across(tidyselect::where(is.double), format_fixture_numeric)
    )
}
write_fixture_tsv <- function(table, path) {
  table |>
    stabilize_fixture_table() |>
    readr::write_tsv(path, na = "")
}

cell_types <- c(
  "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
  "T cells CD4 naive", "T cells CD4 memory resting",
  "T cells CD4 memory activated", "T cells follicular helper",
  "T cells regulatory (Tregs)", "T cells gamma delta",
  "NK cells resting", "NK cells activated", "Monocytes",
  "Macrophages M0", "Macrophages M1", "Macrophages M2",
  "Dendritic cells resting", "Dendritic cells activated",
  "Mast cells resting", "Mast cells activated", "Eosinophils", "Neutrophils"
)
cell_groups <- c(
  "B cells", "CD4 T cells", "CD8 T cells", "Gamma-delta T cells",
  "NK cells", "Monocyte/myeloid", "Neutrophils", "Eosinophils",
  "Dendritic cells", "Mast cells"
)
sample_ids <- sprintf("SAMPLE_%02d", seq_len(12L))
signature_gene_names <- sprintf("SYN_GENE_%03d", seq_len(66L))
signature_gene_ids <- sprintf("ENSGSYN%06d", seq_len(66L))
extra_gene_names <- sprintf("EXTRA_GENE_%02d", seq_len(6L))
extra_gene_ids <- sprintf("ENSGEXTRA%06d", seq_len(6L))

signature <- matrix(
  stats::runif(66L * 22L, min = 0.25, max = 1.25),
  nrow = 66L,
  ncol = 22L,
  dimnames = list(signature_gene_names, cell_types)
)
for (cell_type_index in seq_along(cell_types)) {
  marker_indices <- (3L * cell_type_index - 2L):(3L * cell_type_index)
  signature[marker_indices, cell_type_index] <-
    signature[marker_indices, cell_type_index] + 60
}

weights <- matrix(
  stats::rgamma(length(sample_ids) * length(cell_types), shape = 2, rate = 1),
  nrow = length(sample_ids),
  ncol = length(cell_types),
  dimnames = list(sample_ids, cell_types)
)
weights <- sweep(weights, 1L, rowSums(weights), "/")

extra_gene_expression <- matrix(
  stats::runif(length(extra_gene_names) * length(sample_ids), min = 0.5, max = 4),
  nrow = length(extra_gene_names),
  ncol = length(sample_ids),
  dimnames = list(extra_gene_names, sample_ids)
)
duplicate_gene_id <- "ENSGDUP000001"
duplicate_gene_name <- signature_gene_names[[1L]]
duplicate_linear <- matrix(
  stats::runif(length(sample_ids), min = 0.75, max = 3.25),
  nrow = 1L,
  dimnames = list(duplicate_gene_name, sample_ids)
)
bulk_linear <- rbind(
  signature %*% t(weights),
  extra_gene_expression,
  duplicate_linear
)
all_gene_ids <- c(signature_gene_ids, extra_gene_ids, duplicate_gene_id)
all_gene_names <- c(signature_gene_names, extra_gene_names, duplicate_gene_name)
bulk_cpm <- sweep(bulk_linear, 2L, colSums(bulk_linear), "/") * 1e6
stopifnot(
  all(bulk_cpm > 0),
  max(abs(colSums(bulk_cpm) - 1e6)) < 1e-6
)

write_matrix <- function(matrix_value, path, id_column) {
  matrix_value |>
    as.data.frame(check.names = FALSE) |>
    tibble::rownames_to_column(var = id_column) |>
    write_fixture_tsv(path)
}

expression_cpm <- bulk_cpm
rownames(expression_cpm) <- all_gene_ids
coordinates <- tibble::tibble(
  `#chr` = rep(c("chr1", "chr2"), length.out = nrow(expression_cpm)),
  start = as.integer((seq_len(nrow(expression_cpm)) - 1L) * 1000L),
  end = as.integer((seq_len(nrow(expression_cpm)) - 1L) * 1000L + 500L),
  gene_id = rownames(expression_cpm)
)
expression_bed <- dplyr::bind_cols(
  coordinates,
  tibble::as_tibble(expression_cpm, .name_repair = "minimal")
)
write_fixture_tsv(
  expression_bed,
  file.path(output_directory, "synthetic_expression.bed")
)
expression_with_zero <- expression_cpm
expression_with_zero[signature_gene_ids[[1L]], sample_ids[[1L]]] <- 0
expression_with_zero_bed <- dplyr::bind_cols(
  coordinates,
  tibble::as_tibble(expression_with_zero, .name_repair = "minimal")
)
write_fixture_tsv(
  expression_with_zero_bed,
  file.path(output_directory, "synthetic_expression_with_zero.bed")
)
write_matrix(
  signature,
  file.path(output_directory, "synthetic_signature.tsv"),
  "gene_symbol"
)
write_matrix(
  weights,
  file.path(output_directory, "precomputed_proportions.tsv"),
  "sample_id"
)

additional_covariates <- tibble::tibble(
  sample_id = sample_ids,
  batch_indicator = rep(c(0, 1), each = 6L),
  genotype_pc1 = seq(-1.1, 1.1, length.out = length(sample_ids))
)
write_fixture_tsv(
  additional_covariates,
  file.path(output_directory, "additional_covariates.tsv")
)
writeLines(sample_ids, file.path(output_directory, "samples.tsv"))

gene_types <- ifelse(
  seq_along(all_gene_ids) %% 5L == 0L,
  "lncRNA",
  "protein_coding"
)
gene_types[all_gene_ids == "ENSGEXTRA000006"] <- "processed_pseudogene"
gtf_lines <- purrr::pmap_chr(
  list(
    gene_id = all_gene_ids,
    gene_name = all_gene_names,
    gene_type = gene_types,
    chromosome = coordinates[["#chr"]],
    start = coordinates$start + 1L,
    end = coordinates$end
  ),
  function(gene_id, gene_name, gene_type, chromosome, start, end) {
    attributes <- sprintf(
      'gene_id "%s"; gene_name "%s"; gene_type "%s";',
      gene_id,
      gene_name,
      gene_type
    )
    paste(
      chromosome, "synthetic", "gene", start, end, ".", "+", ".",
      attributes,
      sep = "\t"
    )
  }
)
writeLines(gtf_lines, file.path(output_directory, "synthetic.gtf"))

writeLines(cell_groups, file.path(output_directory, "expected_groups.txt"))

fixture_root <- "tests/cell_type_specific_expression/fixtures"
common_inputs <- list(
  "PrepareCellTypeEqtlWorkflow.gtf" = file.path(fixture_root, "synthetic.gtf"),
  "PrepareCellTypeEqtlWorkflow.SampleList" =
    file.path(fixture_root, "samples.tsv"),
  "PrepareCellTypeEqtlWorkflow.AdditionalCovariates" =
    file.path(fixture_root, "additional_covariates.tsv"),
  "PrepareCellTypeEqtlWorkflow.deconvolution_docker_image" =
    "cell-type-specific-expression:test",
  "PrepareCellTypeEqtlWorkflow.qtl_docker_image" = "prepare-qtl:test",
  "PrepareCellTypeEqtlWorkflow.preemptible_attempts" = 0L,
  "PrepareCellTypeEqtlWorkflow.max_retries" = 0L,
  "PrepareCellTypeEqtlWorkflow.min_lm22_overlap" = 0.80,
  "PrepareCellTypeEqtlWorkflow.dtangle_marker_fraction" = 0.10,
  "PrepareCellTypeEqtlWorkflow.dtangle_quantile_normalize" = FALSE,
  "PrepareCellTypeEqtlWorkflow.group_mean_threshold" = 0.0001,
  "PrepareCellTypeEqtlWorkflow.zero_floor" = 0.000001,
  "PrepareCellTypeEqtlWorkflow.tca_max_iters" = 10L,
  "PrepareCellTypeEqtlWorkflow.tca_parallel" = FALSE,
  "PrepareCellTypeEqtlWorkflow.gene_type" = c("protein_coding", "lncRNA"),
  "PrepareCellTypeEqtlWorkflow.random_seed" = seed,
  "PrepareCellTypeEqtlWorkflow.dtangle_cpu" = 2L,
  "PrepareCellTypeEqtlWorkflow.dtangle_memory" = "8 GB",
  "PrepareCellTypeEqtlWorkflow.dtangle_disk_gb" = 20L,
  "PrepareCellTypeEqtlWorkflow.proportions_cpu" = 1L,
  "PrepareCellTypeEqtlWorkflow.proportions_memory" = "4 GB",
  "PrepareCellTypeEqtlWorkflow.proportions_disk_gb" = 10L,
  "PrepareCellTypeEqtlWorkflow.fit_cpu" = 2L,
  "PrepareCellTypeEqtlWorkflow.fit_memory" = "8 GB",
  "PrepareCellTypeEqtlWorkflow.fit_disk_gb" = 20L,
  "PrepareCellTypeEqtlWorkflow.export_cpu" = 2L,
  "PrepareCellTypeEqtlWorkflow.export_memory" = "8 GB",
  "PrepareCellTypeEqtlWorkflow.export_disk_gb" = 20L,
  "PrepareCellTypeEqtlWorkflow.manifest_cpu" = 1L,
  "PrepareCellTypeEqtlWorkflow.manifest_memory" = "4 GB",
  "PrepareCellTypeEqtlWorkflow.manifest_disk_gb" = 10L,
  "PrepareCellTypeEqtlWorkflow.scatter_cpu" = 1L,
  "PrepareCellTypeEqtlWorkflow.scatter_memory" = "4 GB",
  "PrepareCellTypeEqtlWorkflow.scatter_disk_gb" = 10L,
  "PrepareCellTypeEqtlWorkflow.eqtl_cpu" = 2L,
  "PrepareCellTypeEqtlWorkflow.eqtl_memory" = 8L,
  "PrepareCellTypeEqtlWorkflow.eqtl_disk_gb" = 20L
)
dtangle_inputs <- append(
  common_inputs,
  list(
    "PrepareCellTypeEqtlWorkflow.expression" =
      file.path(fixture_root, "synthetic_expression.bed"),
    "PrepareCellTypeEqtlWorkflow.lm22" =
      file.path(fixture_root, "synthetic_signature.tsv"),
    "PrepareCellTypeEqtlWorkflow.OutputPrefix" = "synthetic.dtangle",
    "PrepareCellTypeEqtlWorkflow.log2_pseudocount" = 0.0
  ),
  after = 0L
)
precomputed_inputs <- append(
  common_inputs,
  list(
    "PrepareCellTypeEqtlWorkflow.expression" =
      file.path(fixture_root, "synthetic_expression_with_zero.bed"),
    "PrepareCellTypeEqtlWorkflow.lm22" =
      file.path(fixture_root, "synthetic_signature.tsv"),
    "PrepareCellTypeEqtlWorkflow.precomputed_proportions" =
      file.path(fixture_root, "precomputed_proportions.tsv"),
    "PrepareCellTypeEqtlWorkflow.OutputPrefix" = "synthetic.precomputed",
    "PrepareCellTypeEqtlWorkflow.log2_pseudocount" = 1.0
  ),
  after = 0L
)
jsonlite::write_json(
  dtangle_inputs,
  file.path(output_directory, "dtangle-e2e.inputs.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)
jsonlite::write_json(
  precomputed_inputs,
  file.path(output_directory, "precomputed-e2e.inputs.json"),
  auto_unbox = TRUE,
  pretty = TRUE
)

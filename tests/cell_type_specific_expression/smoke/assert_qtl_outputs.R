#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2L || length(arguments) > 3L) {
  stop(
    "Usage: assert_qtl_outputs.R OUTPUTS_JSON INPUTS_JSON [FIXTURE_DIRECTORY]",
    call. = FALSE
  )
}

outputs_path <- arguments[[1L]]
inputs_path <- arguments[[2L]]
fixture_directory <- if (length(arguments) == 3L) {
  arguments[[3L]]
} else {
  "tests/cell_type_specific_expression/fixtures"
}
workflow_name <- "PrepareCellTypeEqtlWorkflow"
outputs <- jsonlite::read_json(outputs_path, simplifyVector = TRUE)
inputs <- jsonlite::read_json(inputs_path, simplifyVector = TRUE)

require_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

input_value <- function(name) {
  key <- paste(workflow_name, name, sep = ".")
  value <- inputs[[key]]
  if (is.null(value) || length(value) == 0L) {
    stop(sprintf("The required workflow input is absent: %s", key), call. = FALSE)
  }
  value
}

output_value <- function(name) {
  key <- paste(workflow_name, name, sep = ".")
  value <- outputs[[key]]
  if (is.null(value) || length(value) == 0L) {
    stop(sprintf("The required workflow output is absent: %s", key), call. = FALSE)
  }
  value
}

expected_groups <- readLines(
  file.path(fixture_directory, "expected_groups.txt"),
  warn = FALSE
)
expected_samples <- readLines(
  file.path(fixture_directory, "samples.tsv"),
  warn = FALSE
)
manifest_columns <- c(
  "entity:cell_type_id", "cell_type", "cell_type_slug", "int_bed", "scaled_bed",
  "int_phenotype_pcs", "int_phenotype_pcs_all",
  "scaled_phenotype_pcs", "scaled_phenotype_pcs_all",
  "int_merged_covariates", "scaled_merged_covariates",
  "int_connectivity_outliers", "scaled_connectivity_outliers"
)
manifest <- readr::read_tsv(
  output_value("cell_type_qtl_manifest"),
  col_types = readr::cols(.default = readr::col_character()),
  name_repair = "minimal",
  show_col_types = FALSE,
  progress = FALSE
)
require_true(
  identical(names(manifest), manifest_columns),
  "The QTL manifest must have the exact 13-column schema"
)
require_true(
  identical(manifest$cell_type, expected_groups),
  "The QTL manifest must have one ordered row for each expected group"
)
require_true(
  nrow(manifest) == length(expected_groups),
  "The QTL manifest row count is incorrect"
)
require_true(
  identical(output_value("cell_types"), manifest$cell_type),
  "The authoritative cell-type array does not match the QTL manifest"
)
require_true(
  identical(output_value("cell_type_slugs"), manifest$cell_type_slug),
  "The authoritative cell-type slug array does not match the QTL manifest"
)
require_true(
  all(grepl("^[a-z0-9]+(_[a-z0-9]+)*$", manifest$cell_type_slug)) &&
    anyDuplicated(manifest$cell_type_slug) == 0L,
  "The QTL manifest slugs must be safe and unique"
)
require_true(
  identical(manifest[["entity:cell_type_id"]], manifest$cell_type_slug),
  "The Terra entity IDs must match the ordered cell-type slugs"
)
file_columns <- manifest_columns[-c(1L, 2L, 3L)]
file_values <- unlist(manifest[file_columns], use.names = FALSE)
require_true(
  all(nzchar(file_values)) &&
    all(grepl("^(/|gs://[^/]+/)", file_values)),
  "The QTL manifest file fields must contain full paths"
)
require_true(
  !any(grepl("residual", names(manifest), ignore.case = TRUE)),
  "The QTL manifest must not contain residualized file columns"
)

array_outputs <- c(
  int_bed = "int_beds",
  scaled_bed = "scaled_beds",
  int_phenotype_pcs = "int_phenotype_pcs",
  int_phenotype_pcs_all = "int_phenotype_pcs_all",
  scaled_phenotype_pcs = "scaled_phenotype_pcs",
  scaled_phenotype_pcs_all = "scaled_phenotype_pcs_all",
  int_merged_covariates = "int_merged_covariates",
  scaled_merged_covariates = "scaled_merged_covariates",
  int_connectivity_outliers = "int_connectivity_outliers",
  scaled_connectivity_outliers = "scaled_connectivity_outliers"
)
authoritative_paths <- purrr::imap(array_outputs, function(output_name, column_name) {
  paths <- output_value(output_name)
  require_true(
    length(paths) == length(expected_groups),
    sprintf("The %s authoritative array has an incorrect length", output_name)
  )
  require_true(
    identical(
      normalizePath(paths, mustWork = TRUE),
      normalizePath(manifest[[column_name]], mustWork = TRUE)
    ),
    sprintf("The %s authoritative array does not match the manifest", output_name)
  )
  require_true(
    all(file.exists(paths)) && anyDuplicated(paths) == 0L,
    sprintf("The %s authoritative array has an absent or duplicate file", output_name)
  )
  paths
})

scatter_prefixes <- paste(
  as.character(input_value("OutputPrefix")),
  manifest$cell_type_slug,
  sep = "."
)
expected_manifest_files <- list(
  int_bed = paste0(scatter_prefixes, ".expression.INT.bed.gz"),
  scaled_bed = paste0(scatter_prefixes, ".expression.scaled.bed.gz"),
  int_phenotype_pcs = paste0(scatter_prefixes, ".expression_phenotype_PCs.INT.tsv"),
  int_phenotype_pcs_all = paste0(
    scatter_prefixes,
    ".expression_phenotype_PCs.INT.all.tsv"
  ),
  scaled_phenotype_pcs = paste0(
    scatter_prefixes,
    ".expression_phenotype_PCs.scaled.tsv"
  ),
  scaled_phenotype_pcs_all = paste0(
    scatter_prefixes,
    ".expression_phenotype_PCs.scaled.all.tsv"
  ),
  int_merged_covariates = paste0(
    scatter_prefixes,
    ".expression_QTL_covariates.INT.tsv"
  ),
  scaled_merged_covariates = paste0(
    scatter_prefixes,
    ".expression_QTL_covariates.scaled.tsv"
  ),
  int_connectivity_outliers = paste0(
    scatter_prefixes,
    ".expression.INT.connectivity_outliers.tsv"
  ),
  scaled_connectivity_outliers = paste0(
    scatter_prefixes,
    ".expression.scaled.connectivity_outliers.tsv"
  )
)
purrr::iwalk(expected_manifest_files, function(expected_files, column_name) {
  require_true(
    identical(basename(manifest[[column_name]]), expected_files),
    sprintf("The %s manifest basenames do not match the output contract", column_name)
  )
})

read_qtl_bed <- function(path, label) {
  bed <- readr::read_tsv(
    path,
    col_types = readr::cols(
      `#chr` = readr::col_character(),
      start = readr::col_integer(),
      end = readr::col_integer(),
      gene_id = readr::col_character(),
      .default = readr::col_double()
    ),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  metadata_columns <- c("#chr", "start", "end", "gene_id")
  require_true(
    identical(names(bed)[seq_along(metadata_columns)], metadata_columns),
    sprintf("The %s BED has an invalid schema: %s", label, basename(path))
  )
  sample_columns <- names(bed)[-seq_along(metadata_columns)]
  require_true(
    nrow(bed) > 0L && length(sample_columns) > 0L &&
      all(sample_columns %in% expected_samples),
    sprintf("The %s BED has invalid dimensions or samples: %s", label, basename(path))
  )
  require_true(
    !anyNA(bed$gene_id) && all(nzchar(bed$gene_id)) &&
      anyDuplicated(bed$gene_id) == 0L,
    sprintf("The %s BED has invalid gene IDs: %s", label, basename(path))
  )
  values <- bed |>
    dplyr::select(dplyr::all_of(sample_columns)) |>
    as.matrix()
  require_true(
    all(is.finite(values)),
    sprintf("The %s BED has non-finite values: %s", label, basename(path))
  )
  bed
}

int_bed_tables <- purrr::map(authoritative_paths$int_bed, read_qtl_bed, "INT")
scaled_bed_tables <- purrr::map(
  authoritative_paths$scaled_bed,
  read_qtl_bed,
  "scaled"
)
purrr::walk2(int_bed_tables, scaled_bed_tables, function(int_bed, scaled_bed) {
  require_true(
    identical(int_bed$gene_id, scaled_bed$gene_id),
    "The paired INT and scaled BED gene orders differ"
  )
})

# Compare the final QTL values with the exported linear-CPM matrices. Compute
# scaling on all selected samples first, then select the non-outlier samples.
# This catches a pipeline that passes CPM through Log2CpmBed and skips the log.
tca_paths <- output_value("cell_type_beds")
purrr::walk(seq_len(nrow(manifest)), function(index) {
  tca_path <- tca_paths[basename(tca_paths) == paste0(manifest$cell_type_slug[[index]], ".bed.gz")]
  require_true(length(tca_path) == 1L, "Each QTL cell type must match one exported TCA BED")
  tca_bed <- read_qtl_bed(tca_path, "TCA CPM")
  scaled_bed <- scaled_bed_tables[[index]]
  require_true(identical(tca_bed$gene_id, scaled_bed$gene_id),
               "The QTL scaled BED changed the TCA gene order")
  cpm <- tca_bed |>
    dplyr::select(dplyr::all_of(expected_samples)) |>
    as.matrix()
  require_true(all(cpm >= 0), "CPM inputs must not contain negative TCA estimates")
  logged <- log2(cpm + 1)
  centered <- sweep(logged, 1L, rowMeans(logged), "-")
  deviations <- sqrt(rowSums(centered^2) / (ncol(centered) - 1L))
  expected <- sweep(centered, 1L, deviations, "/")
  kept_samples <- names(scaled_bed)[-(1:4)]
  observed <- scaled_bed |>
    dplyr::select(dplyr::all_of(kept_samples)) |>
    as.matrix()
  require_true(
    isTRUE(all.equal(unname(observed), unname(expected[, kept_samples, drop = FALSE]),
                     tolerance = 1e-7)),
    "The QTL scaled BED must use log2(CPM + 1) before centering and scaling"
  )
})

read_pc_table <- function(path, expected_bed, label) {
  table <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  require_true(
    ncol(table) >= 1L && identical(names(table)[[1L]], "ID"),
    sprintf("The %s PC file must start with ID: %s", label, basename(path))
  )
  bed_samples <- names(expected_bed)[-(1:4)]
  require_true(
    identical(table$ID, bed_samples),
    sprintf("The %s PC sample order differs from its BED: %s", label, basename(path))
  )
  if (ncol(table) > 1L) {
    require_true(
      all(grepl("^PC[0-9]+$", names(table)[-1L])),
      sprintf("The %s PC names are invalid: %s", label, basename(path))
    )
    values <- table[-1L] |>
      dplyr::mutate(dplyr::across(dplyr::everything(), readr::parse_double)) |>
      as.matrix()
    require_true(
      all(is.finite(values)),
      sprintf("The %s PC values are not finite: %s", label, basename(path))
    )
  }
  table
}

check_pc_pair <- function(selected_path, all_path, bed, label) {
  selected <- read_pc_table(selected_path, bed, paste(label, "selected"))
  complete <- read_pc_table(all_path, bed, paste(label, "complete"))
  require_true(
    all(names(selected) %in% names(complete)),
    sprintf("The selected %s PCs are not a subset of the complete PCs", label)
  )
  require_true(
    isTRUE(all.equal(
      as.data.frame(selected),
      as.data.frame(complete[names(selected)]),
      check.attributes = FALSE
    )),
    sprintf("The selected %s PC values differ from the complete PCs", label)
  )
}

purrr::pwalk(
  list(
    authoritative_paths$int_phenotype_pcs,
    authoritative_paths$int_phenotype_pcs_all,
    int_bed_tables
  ),
  ~ check_pc_pair(..1, ..2, ..3, "INT")
)
purrr::pwalk(
  list(
    authoritative_paths$scaled_phenotype_pcs,
    authoritative_paths$scaled_phenotype_pcs_all,
    scaled_bed_tables
  ),
  ~ check_pc_pair(..1, ..2, ..3, "scaled")
)

check_merged_covariates <- function(path, bed, label) {
  table <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  require_true(
    nrow(table) >= 2L && identical(names(table)[[1L]], "ID"),
    sprintf("The %s merged covariates have an invalid schema", label)
  )
  bed_samples <- names(bed)[-(1:4)]
  require_true(
    identical(names(table)[-1L], bed_samples),
    sprintf("The %s merged covariate sample order differs from its BED", label)
  )
  require_true(
    all(c("batch_indicator", "genotype_pc1") %in% table$ID),
    sprintf("The %s merged covariates omit required fixture covariates", label)
  )
  values <- table[-1L] |>
    dplyr::mutate(dplyr::across(dplyr::everything(), readr::parse_double)) |>
    as.matrix()
  require_true(
    all(is.finite(values)),
    sprintf("The %s merged covariate values are not finite", label)
  )
}

purrr::pwalk(
  list(authoritative_paths$int_merged_covariates, int_bed_tables),
  ~ check_merged_covariates(..1, ..2, "INT")
)
purrr::pwalk(
  list(authoritative_paths$scaled_merged_covariates, scaled_bed_tables),
  ~ check_merged_covariates(..1, ..2, "scaled")
)

check_outlier_report <- function(path, bed, label) {
  table <- readr::read_tsv(
    path,
    col_types = readr::cols(
      SampleID = readr::col_character(),
      Z_score = readr::col_double()
    ),
    show_col_types = FALSE,
    progress = FALSE
  )
  require_true(
    identical(names(table), c("SampleID", "Z_score")),
    sprintf("The %s outlier report has an invalid schema", label)
  )
  require_true(
    all(table$SampleID %in% expected_samples) && all(is.finite(table$Z_score)),
    sprintf("The %s outlier report has invalid values", label)
  )
  require_true(
    anyDuplicated(table$SampleID) == 0L,
    sprintf("The %s outlier report has duplicate sample IDs", label)
  )
  retained_samples <- names(bed)[-(1:4)]
  require_true(
    length(intersect(table$SampleID, retained_samples)) == 0L,
    sprintf("The %s outlier report contains a retained sample", label)
  )
  require_true(
    identical(sort(c(table$SampleID, retained_samples)), sort(expected_samples)),
    sprintf("The %s retained and outlier samples do not partition the cohort", label)
  )
}

purrr::pwalk(
  list(authoritative_paths$int_connectivity_outliers, int_bed_tables),
  ~ check_outlier_report(..1, ..2, "INT")
)
purrr::pwalk(
  list(authoritative_paths$scaled_connectivity_outliers, scaled_bed_tables),
  ~ check_outlier_report(..1, ..2, "scaled")
)

require_true(
  all(manifest$int_bed != manifest$scaled_bed),
  "The INT and scaled BED basenames must differ"
)
require_true(
  all(manifest$int_phenotype_pcs != manifest$scaled_phenotype_pcs) &&
    all(manifest$int_phenotype_pcs_all != manifest$scaled_phenotype_pcs_all) &&
    all(manifest$int_merged_covariates != manifest$scaled_merged_covariates) &&
    all(manifest$int_connectivity_outliers != manifest$scaled_connectivity_outliers),
  "The INT and scaled QTL output basenames must differ"
)

expression <- readr::read_tsv(
  input_value("expression"),
  show_col_types = FALSE,
  progress = FALSE
)
expression_values <- expression |>
  dplyr::select(-dplyr::all_of(c("#chr", "start", "end", "gene_id"))) |>
  as.matrix()
pseudocount <- as.numeric(input_value("log2_pseudocount"))
effective_parameters <- jsonlite::read_json(
  output_value("effective_parameters_file"),
  simplifyVector = TRUE
)
require_true(
  as.numeric(effective_parameters$log2_pseudocount) == pseudocount,
  "The effective pseudocount does not match the workflow input"
)
if (pseudocount > 0) {
  require_true(
    any(expression_values == 0) && all(expression_values >= 0),
    "The positive-pseudocount smoke run must accept a zero CPM input"
  )
} else {
  require_true(
    all(expression_values > 0),
    "The zero-pseudocount smoke run must use positive CPM input"
  )
}

required_upstream_outputs <- c(
  "proportion_mode_validation_log", "proportions_lm22",
  "proportions_combined", "tca_weights", "cell_group_filter_report",
  "proportions_log", "tca_model", "tca_model_unfiltered", "tca_numerical_excluded_genes",
  "tca_cleanup_log", "tca_model_log",
  "tca_excluded_genes", "fit_tca_log", "cell_type_bed_inventory",
  "reconstruction_by_sample", "qc_summary", "qc_plots", "export_log",
  "export_detail_log", "output_manifest", "output_inventory", "manifest_log",
  "effective_parameters_file", "cell_type_qtl_manifest_log"
)
purrr::walk(required_upstream_outputs, function(name) {
  path <- output_value(name)
  require_true(file.exists(path), sprintf("The upstream output is absent: %s", name))
  require_true(file.info(path)$size > 0, sprintf("The upstream output is empty: %s", name))
})
upstream_beds <- output_value("cell_type_beds")
require_true(
  length(upstream_beds) == length(expected_groups) && all(file.exists(upstream_beds)),
  "The upstream cell-type BED array is incomplete"
)

message(sprintf(
  "QTL smoke assertions passed: groups=%d manifest_columns=%d pseudocount=%g",
  length(expected_groups),
  length(manifest_columns),
  pseudocount
))

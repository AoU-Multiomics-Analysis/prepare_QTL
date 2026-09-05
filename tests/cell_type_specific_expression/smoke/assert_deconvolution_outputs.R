#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) < 2L || length(arguments) > 3L) {
  stop(
    paste(
      "Usage: assert_deconvolution_outputs.R",
      "OUTPUTS_JSON INPUTS_JSON [FIXTURE_DIRECTORY]"
    ),
    call. = FALSE
  )
}

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_argument) != 1L) {
  stop("The assertion script path is unavailable", call. = FALSE)
}
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
source(file.path(dirname(script_path), "deconvolution_assertion_helpers.R"))

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

expected_tca_parallel <- input_value("tca_parallel")
require_true(
  is.logical(expected_tca_parallel) &&
    length(expected_tca_parallel) == 1L &&
    !is.na(expected_tca_parallel),
  "The tca_parallel workflow input must be true or false"
)

read_matrix_path <- function(path, id_column, label) {
  table <- readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    show_col_types = FALSE,
    progress = FALSE
  )
  require_true(
    identical(names(table)[[1L]], id_column),
    sprintf("The %s table must start with %s", label, id_column)
  )
  values <- table[-1L] |>
    dplyr::mutate(dplyr::across(dplyr::everything(), readr::parse_double)) |>
    as.matrix()
  require_true(all(is.finite(values)), sprintf("The %s table must be finite", label))
  rownames(values) <- table[[id_column]]
  values
}

read_matrix_table <- function(name, id_column) {
  read_matrix_path(output_value(name), id_column, name)
}

expression_path <- input_value("expression")
expression <- readr::read_tsv(
  expression_path,
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
  identical(names(expression)[seq_along(metadata_columns)], metadata_columns),
  "The fixture expression BED has an invalid schema"
)
expected_samples <- readLines(
  file.path(fixture_directory, "samples.tsv"),
  warn = FALSE
)
gene_type_filter_report <- readr::read_tsv(
  output_value("gene_type_filter_report"),
  show_col_types = FALSE,
  progress = FALSE
)
require_true(
  identical(
    names(gene_type_filter_report),
    c("gene_id", "gene_name", "gene_type", "retained", "filter_reason")
  ),
  "The gene-type filter report has an invalid schema"
)
require_true(
  identical(gene_type_filter_report$gene_id, expression$gene_id),
  "The gene-type filter report changed the input gene order"
)
selected_gene_types <- input_value("gene_type")
expected_retained <- gene_type_filter_report$gene_type %in% selected_gene_types
require_true(
  identical(gene_type_filter_report$retained, expected_retained),
  "The gene-type filter report retained an unexpected gene"
)
require_true(
  any(!gene_type_filter_report$retained) &&
    all(
      gene_type_filter_report$filter_reason[!gene_type_filter_report$retained] ==
        "gene_type_not_selected"
    ),
  "The smoke fixture must exclude at least one unselected gene type"
)
expected_gene_ids <- gene_type_filter_report$gene_id[expected_retained]
expected_coordinates <- expression |>
  dplyr::filter(.data$gene_id %in% expected_gene_ids) |>
  dplyr::select(dplyr::all_of(metadata_columns))
filtered_expression <- readr::read_tsv(
  output_value("filtered_expression"),
  show_col_types = FALSE,
  progress = FALSE
)
require_true(
  identical(filtered_expression$gene_id, expected_gene_ids),
  "The filtered expression BED has an incorrect gene order"
)
expected_groups <- readLines(
  file.path(fixture_directory, "expected_groups.txt"),
  warn = FALSE
)
expression_values <- expression |>
  dplyr::select(-dplyr::all_of(metadata_columns)) |>
  as.matrix()
log2_pseudocount <- as.numeric(input_value("log2_pseudocount"))
require_true(
  identical(names(expression)[-seq_along(metadata_columns)], expected_samples),
  "The expression sample order does not match samples.tsv"
)
require_true(all(is.finite(expression_values)), "The input CPM values must be finite")
require_true(all(expression_values >= 0), "The input CPM values must be non-negative")
lm22_key <- paste(workflow_name, "lm22", sep = ".")
precomputed_key <- paste(workflow_name, "precomputed_proportions", sep = ".")
has_lm22 <- !is.null(inputs[[lm22_key]])
has_precomputed <- !is.null(inputs[[precomputed_key]])
require_true(has_lm22, "The fixture must provide the required LM22 reference")
expected_proportion_mode <- if (has_precomputed) "precomputed" else "hspe"
if (identical(expected_proportion_mode, "hspe") && log2_pseudocount == 0) {
  require_true(all(expression_values > 0), "Pseudocount zero requires positive CPM")
}

proportions <- read_matrix_table("proportions_lm22", "sample_id")
combined <- read_matrix_table("proportions_combined", "sample_id")
tca_weights <- read_matrix_table("tca_weights", "sample_id")
proportion_value_tolerance <- 1e-10
derived_value_tolerance <- 1e-10
authoritative_proportion_path <- if (identical(expected_proportion_mode, "hspe")) {
  output_value("estimated_proportions")
} else {
  input_value("precomputed_proportions")
}
authoritative_proportions <- read_matrix_path(
  authoritative_proportion_path,
  "sample_id",
  paste(expected_proportion_mode, "authoritative proportions")
)
require_matrix_equal(
  proportions,
  authoritative_proportions,
  proportion_value_tolerance,
  "LM22 proportions"
)
mean_threshold <- as.numeric(input_value("group_mean_threshold"))
zero_floor <- as.numeric(input_value("zero_floor"))
expected_proportion_outputs <- derive_expected_proportion_outputs(
  authoritative_proportions,
  mean_threshold,
  zero_floor
)
require_matrix_equal(
  combined,
  expected_proportion_outputs$combined,
  derived_value_tolerance,
  "combined proportions"
)
require_matrix_equal(
  tca_weights,
  expected_proportion_outputs$tca_weights,
  derived_value_tolerance,
  "TCA weights"
)
signature_columns <- readr::read_tsv(
  file.path(fixture_directory, "synthetic_signature.tsv"),
  n_max = 0L,
  show_col_types = FALSE,
  progress = FALSE
) |>
  names()
expected_lm22_types <- signature_columns[-1L]
require_true(all(proportions >= 0), "The LM22 proportions must be non-negative")
require_true(all(combined >= 0), "The combined proportions must be non-negative")
require_true(all(tca_weights > 0), "The TCA weights must be positive")
require_row_sums_within_tolerance(
  proportions,
  workflow_input_proportion_row_sum_tolerance,
  "LM22 proportion"
)
require_row_sums_within_tolerance(
  combined,
  workflow_input_proportion_row_sum_tolerance,
  "combined proportion"
)
require_true(
  max(abs(rowSums(tca_weights) - 1)) < 1e-8,
  "The TCA weight rows must sum to one"
)
require_true(
  identical(colnames(proportions), expected_lm22_types),
  "The LM22 columns are not authoritative"
)
purrr::walk(
  list(proportions, combined, tca_weights),
  ~ require_true(
    identical(rownames(.x), expected_samples),
    "The proportion sample order does not match samples.tsv"
  )
)
require_true(
  identical(colnames(combined), names(canonical_lm22_group_map())),
  "The combined group order does not match the canonical LM22 mapping"
)
require_true(
  identical(expected_proportion_outputs$retained_groups, expected_groups),
  "The configured threshold did not retain the expected group order"
)
require_true(
  identical(colnames(tca_weights), expected_groups),
  "The TCA group order does not match expected_groups.txt"
)
require_true(
  identical(colnames(combined), expected_groups),
  "The combined group order does not match expected_groups.txt"
)

filter_report <- readr::read_tsv(
  output_value("cell_group_filter_report"),
  show_col_types = FALSE,
  progress = FALSE
)
require_true(
  identical(filter_report$cell_group, names(canonical_lm22_group_map())),
  "The filter report group order is not canonical"
)
require_true(
  max(abs(filter_report$cohort_mean - expected_proportion_outputs$cohort_means)) <=
    derived_value_tolerance,
  "The filter report cohort means are incorrect"
)
require_true(
  all(filter_report$threshold == mean_threshold) &&
    identical(filter_report$retained, unname(expected_proportion_outputs$retained)) &&
    identical(
      as.integer(filter_report$zero_count_before),
      as.integer(expected_proportion_outputs$zero_counts)
    ) &&
    all(filter_report$zero_floor == zero_floor),
  "The filter report threshold or zero-floor values are incorrect"
)
retained_groups <- filter_report |>
  dplyr::filter(.data$retained) |>
  dplyr::pull(.data$cell_group)
require_true(
  identical(retained_groups, expected_groups),
  "The retained groups do not match expected_groups.txt"
)

reconstruction <- readr::read_tsv(
  output_value("reconstruction_by_sample"),
  show_col_types = FALSE,
  progress = FALSE
)
require_true(
  identical(reconstruction$sample_id, expected_samples),
  "The reconstruction sample order does not match samples.tsv"
)
require_true(
  all(reconstruction$gene_count == length(expected_gene_ids)),
  "The reconstruction gene count is incorrect"
)
require_true(
  all(is.finite(reconstruction$correlation)) && all(is.finite(reconstruction$rmse)),
  "The reconstruction metrics must be finite"
)

cell_type_bed_paths <- normalizePath(output_value("cell_type_beds"))
require_true(
  length(cell_type_bed_paths) == length(expected_groups),
  "The cell-type BED count is incorrect"
)
require_true(
  all(grepl("[.]bed[.]gz$", cell_type_bed_paths)) &&
    anyDuplicated(cell_type_bed_paths) == 0L,
  "The cell-type BED paths must be unique BED gzip files"
)
purrr::walk(cell_type_bed_paths, function(path) {
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
  require_true(
    identical(names(bed)[seq_along(metadata_columns)], metadata_columns),
    sprintf("The cell-type BED has an invalid schema: %s", basename(path))
  )
  require_true(
    identical(bed$gene_id, expected_gene_ids) &&
      identical(names(bed)[-seq_along(metadata_columns)], expected_samples),
    sprintf("The cell-type BED has an invalid identity order: %s", basename(path))
  )
  require_true(
    identical(bed[["#chr"]], expected_coordinates[["#chr"]]) &&
      identical(bed$start, expected_coordinates$start) &&
      identical(bed$end, expected_coordinates$end) &&
      identical(bed$gene_id, expected_coordinates$gene_id),
    sprintf("The cell-type BED has changed coordinates: %s", basename(path))
  )
  values <- bed |>
    dplyr::select(-dplyr::all_of(metadata_columns)) |>
    as.matrix()
  require_true(
    all(is.finite(values)),
    sprintf("The cell-type BED contains a non-finite value: %s", basename(path))
  )
})

inventory <- readr::read_tsv(
  output_value("cell_type_bed_inventory"),
  show_col_types = FALSE,
  progress = FALSE
)
inventory_basenames <- inventory$path
require_true(nrow(inventory) == length(expected_groups), "The BED inventory row count is incorrect")
require_true(
  identical(inventory$cell_group, expected_groups),
  "The BED inventory group order is incorrect"
)
require_true(
  all(inventory$n_genes == length(expected_gene_ids)) &&
    all(inventory$n_samples == length(expected_samples)) &&
    all(inventory$scale == "cpm"),
  "The BED inventory dimensions or scale are incorrect"
)
require_true(
  all(nzchar(inventory$slug)) && anyDuplicated(inventory$slug) == 0L,
  "The BED inventory slugs must be non-empty and unique"
)
require_true(
  all(grepl("^[^/]+[.]bed[.]gz$", inventory_basenames)) &&
    anyDuplicated(inventory_basenames) == 0L,
  "The BED inventory paths must be stable unique basenames"
)
require_true(
  setequal(basename(cell_type_bed_paths), inventory_basenames),
  "The BED output array does not match the authoritative inventory"
)
require_true(
  all(grepl("^[0-9a-f]{64}$", inventory$sha256)),
  "The BED inventory must contain valid SHA-256 checksums"
)
purrr::pwalk(
  list(inventory_basenames, inventory$sha256),
  function(output_basename, expected_sha256) {
    output_path <- cell_type_bed_paths[
      match(output_basename, basename(cell_type_bed_paths))
    ]
    observed_sha256 <- digest::digest(
      file = output_path,
      algo = "sha256",
      serialize = FALSE
    )
    require_true(
      identical(observed_sha256, expected_sha256),
      sprintf("The BED checksum is incorrect: %s", output_basename)
    )
  }
)
gene_summary <- readr::read_tsv(output_value("cell_type_gene_summary"), show_col_types = FALSE)
require_true(
  identical(gene_summary$cell_type, rep(expected_groups, each = length(expected_gene_ids))) &&
    identical(gene_summary$gene_id, rep(expected_gene_ids, length(expected_groups))) &&
    all(gene_summary$scale == "cpm") && all(gene_summary$n_samples == length(expected_samples)),
  "Gene summary identities, sample counts, or CPM scale are incorrect"
)
purrr::walk(seq_len(nrow(inventory)), function(index) {
  path <- cell_type_bed_paths[match(inventory$path[[index]], basename(cell_type_bed_paths))]
  bed <- readr::read_tsv(path, show_col_types = FALSE)
  values <- as.matrix(bed[, expected_samples])
  observed <- gene_summary |> dplyr::filter(.data$cell_type == inventory$cell_group[[index]])
  require_true(
    identical(observed[["#chr"]], bed[["#chr"]]) &&
      identical(observed$start, bed$start) && identical(observed$end, bed$end),
    "Gene summary coordinates differ from the exported BED"
  )
  expected <- list(
    mean_cpm = rowMeans(values), median_cpm = apply(values, 1, stats::median),
    sd_cpm = apply(values, 1, stats::sd),
    se_mean_cpm = apply(values, 1, stats::sd) / sqrt(ncol(values)),
    q1_cpm = apply(values, 1, stats::quantile, probs = 0.25, type = 7),
    q3_cpm = apply(values, 1, stats::quantile, probs = 0.75, type = 7),
    iqr_cpm = apply(values, 1, stats::IQR, type = 7)
  )
  purrr::iwalk(expected, function(value, name) {
    require_true(
      isTRUE(all.equal(observed[[name]], unname(value), tolerance = 1e-10)),
      sprintf("Gene summary %s differs from exported CPM values", name)
    )
  })
})

public_inventory <- readr::read_tsv(
  output_value("output_inventory"),
  show_col_types = FALSE,
  progress = FALSE
)
require_true(
  isTRUE(all.equal(
    as.data.frame(public_inventory),
    as.data.frame(inventory),
    check.attributes = FALSE
  )),
  "The public inventory differs from the BED inventory"
)

manifest <- jsonlite::read_json(output_value("output_manifest"), simplifyVector = FALSE)
parameters <- manifest$parameters
required_parameter_names <- c(
  "proportion_mode", "log2_pseudocount", "min_lm22_overlap",
  "hspe_marker_fraction", "hspe_marker_method",
  "hspe_quantile_normalize", "group_mean_threshold", "zero_floor",
  "tca_max_iters", "tca_parallel", "gene_type", "random_seed", "scale"
)
numeric_parameter <- function(name) as.numeric(parameters[[name]])
require_true(
  setequal(names(parameters), required_parameter_names),
  "The deconvolution manifest has an invalid parameter schema"
)
require_true(
  identical(parameters$proportion_mode, expected_proportion_mode),
  "The deconvolution manifest has an incorrect proportion mode"
)
require_true(
  numeric_parameter("log2_pseudocount") == log2_pseudocount,
  "The deconvolution manifest has an incorrect log2 pseudocount"
)
require_true(
  numeric_parameter("min_lm22_overlap") == 0.80 &&
    numeric_parameter("hspe_marker_fraction") == 0.10 &&
    identical(parameters$hspe_marker_method, "ratio") &&
    identical(parameters$hspe_quantile_normalize, FALSE) &&
    numeric_parameter("group_mean_threshold") == 0.0001 &&
    numeric_parameter("zero_floor") == 0.000001 &&
    numeric_parameter("tca_max_iters") == 10 &&
    identical(parameters$tca_parallel, expected_tca_parallel) &&
    identical(
      unlist(parameters$gene_type, use.names = FALSE),
      input_value("gene_type")
    ) &&
    numeric_parameter("random_seed") == as.numeric(input_value("random_seed")) &&
    identical(parameters$scale, "cpm") &&
    identical(manifest$tca_version, "1.2.1"),
  "The deconvolution manifest parameters are incorrect"
)
require_true(
  identical(manifest$container_image, input_value("downstream_docker_image")),
  "The manifest did not record its downstream image"
)
require_true(
  length(manifest$outputs) == nrow(inventory),
  "The deconvolution manifest output count is incorrect"
)
purrr::walk(manifest$outputs, function(entry) {
  require_true(
    identical(entry$path, entry$file_name) &&
      grepl("^[^/]+[.]bed[.]gz$", entry$path) &&
      identical(entry$scale, "cpm") &&
      length(entry$sha256) == 1L && grepl("^[0-9a-f]{64}$", entry$sha256),
    "A deconvolution manifest output is not portable"
  )
})

effective_parameters <- jsonlite::read_json(
  output_value("effective_parameters_file"),
  simplifyVector = TRUE
)
require_true(
  identical(sort(names(effective_parameters)), sort(required_parameter_names)),
  "The effective-parameter file has an invalid schema"
)
require_true(
  identical(effective_parameters$proportion_mode, expected_proportion_mode) &&
    as.numeric(effective_parameters$log2_pseudocount) == log2_pseudocount &&
    identical(effective_parameters$tca_parallel, expected_tca_parallel),
  "The effective-parameter file does not match the fixture"
)

tca_model <- readRDS(output_value("tca_model"))
unfiltered_model <- readRDS(output_value("tca_model_unfiltered"))
numerical_exclusions <- readr::read_tsv(output_value("tca_numerical_excluded_genes"),
                                       show_col_types = FALSE, col_types = readr::cols(gene_id = readr::col_character()))
retained_genes <- rownames(tca_model$mus_hat)
require_true(
  identical(tca_model$gene_filter$original_gene_ids, rownames(unfiltered_model$mus_hat)) &&
    identical(tca_model$gene_filter$excluded_gene_ids, numerical_exclusions$gene_id) &&
    identical(retained_genes, setdiff(rownames(unfiltered_model$mus_hat), numerical_exclusions$gene_id)),
  "The final TCA model must contain exactly the numerically retained genes"
)
require_true(
  identical(tca_model$W, unfiltered_model$W) &&
    identical(tca_model$tau_hat, unfiltered_model$tau_hat) &&
    identical(tca_model$mus_hat, unfiltered_model$mus_hat[retained_genes, , drop = FALSE]) &&
    identical(tca_model$sigmas_hat, unfiltered_model$sigmas_hat[retained_genes, , drop = FALSE]),
  "TCA cleanup must preserve the fitted parameters without refitting"
)
for (bed_path in output_value("cell_type_beds")) {
  bed_genes <- readr::read_tsv(bed_path, col_types = readr::cols(.default = readr::col_character()),
                              show_col_types = FALSE)[[4L]]
  require_true(identical(bed_genes, retained_genes),
               "Every exported BED must match the final model gene order")
}
require_true(
  identical(tca_model$tca_parallel, expected_tca_parallel),
  "The TCA model metadata has an incorrect parallel setting"
)
parallel_log_token <- sprintf(
  "parallel=%s",
  tolower(as.character(expected_tca_parallel))
)
fit_log_lines <- readLines(output_value("tca_model_log"), warn = FALSE)
export_log_lines <- readLines(output_value("export_detail_log"), warn = FALSE)
require_true(
  any(grepl(parallel_log_token, fit_log_lines, fixed = TRUE)),
  "The TCA fit log has an incorrect parallel setting"
)
require_true(
  any(grepl(parallel_log_token, export_log_lines, fixed = TRUE)),
  "The TCA export log has an incorrect parallel setting"
)

qc_summary <- readr::read_tsv(
  output_value("qc_summary"),
  show_col_types = FALSE,
  progress = FALSE
)
required_qc_metrics <- c(
  "gene_count", "sample_count", "cell_group_count",
  "excluded_constant_gene_count", "correlation_min", "correlation_median",
  "correlation_mean", "correlation_max", "rmse_min", "rmse_median",
  "rmse_mean", "rmse_max", "lm22_gene_count", "lm22_cell_type_count",
  "lm22_value_validation", "lm22_value_min", "lm22_value_max",
  "input_proportion_max_row_sum_error",
  "combined_proportion_max_row_sum_error", "adjusted_weight_max_row_sum_error",
  "normalization_adjustment_max_abs", "zero_values_adjusted",
  "tca_internal_iterations", "tca_max_internal_iterations", "tca_convergence",
  "tca_tau_hat"
)
require_true(
  setequal(required_qc_metrics, qc_summary$metric),
  "The deconvolution QC summary has an invalid metric set"
)
qc_status <- stats::setNames(qc_summary$status, qc_summary$metric)
expected_lm22_status <- if (expected_proportion_mode == "hspe") {
  "passed"
} else {
  "not_applicable_precomputed_mode"
}
require_true(
  identical(qc_status[["lm22_value_validation"]], expected_lm22_status),
  "The LM22 QC status does not match the proportion mode"
)
require_true(
  qc_status[["tca_convergence"]] %in% c("converged", "max_iterations_reached"),
  "The TCA convergence status is invalid"
)

required_file_outputs <- c(
  "filtered_expression", "gene_type_filter_report", "gene_type_filter_log",
  "proportion_mode_validation_log", "proportions_lm22",
  "proportions_combined", "tca_weights", "cell_group_filter_report",
  "proportions_log", "tca_model", "tca_model_unfiltered", "tca_numerical_excluded_genes",
  "tca_cleanup_log", "tca_model_log",
  "tca_excluded_genes", "fit_tca_log", "cell_type_bed_inventory",
  "reconstruction_by_sample", "qc_summary", "qc_plots", "export_log",
  "export_detail_log", "output_manifest", "output_inventory", "manifest_log",
  "effective_parameters_file", "cell_type_gene_summary", "gene_summary_log"
)
if (identical(expected_proportion_mode, "hspe")) {
  hspe_metadata <- jsonlite::read_json(output_value("hspe_metadata"))
  require_true(
    identical(hspe_metadata$hspe_version, "0.1") &&
      identical(hspe_metadata$optimizer, "DEoptimR") &&
      as.numeric(hspe_metadata$random_seed) == as.numeric(input_value("random_seed")),
    "HSPE metadata must record the source package, optimizer, and workflow seed"
  )
  require_true(
    as.numeric(hspe_metadata$batch_size) == as.numeric(input_value("hspe_batch_size")) &&
      as.numeric(hspe_metadata$batch_count) == ceiling(length(expected_samples) /
        as.numeric(input_value("hspe_batch_size"))) &&
      as.numeric(hspe_metadata$batch_count) > 1,
    "The HSPE smoke run must fit and merge multiple batches, including a partial batch"
  )
  diagnostics <- readr::read_tsv(output_value("hspe_sample_diagnostics"),
                                show_col_types = FALSE)
  require_true(identical(diagnostics$sample_id, expected_samples),
               "HSPE diagnostics must preserve every sample in the original order")
  require_true(all(diagnostics$random_seed > 0 & diagnostics$iterations > 0 &
                     diagnostics$convergence %in% c(0, 1)),
               "HSPE diagnostics must contain valid seeds and optimizer status")
  required_file_outputs <- c(
    required_file_outputs,
    "estimated_proportions", "hspe_markers", "hspe_metadata",
    "hspe_overlap_report", "transformed_lm22", "hspe_log", "hspe_sample_diagnostics"
  )
} else {
  require_true(is.null(outputs[[paste(workflow_name, "hspe_sample_diagnostics", sep = ".")]]),
               "Precomputed mode must skip HSPE batches")
}
purrr::walk(required_file_outputs, function(name) {
  path <- output_value(name)
  require_true(file.exists(path), sprintf("The required output file is absent: %s", name))
  require_true(file.info(path)$size > 0, sprintf("The required output file is empty: %s", name))
})

message(sprintf(
  paste(
    "Deconvolution smoke assertions passed:",
    "mode=%s genes=%d samples=%d groups=%d pseudocount=%g"
  ),
  expected_proportion_mode,
  length(expected_gene_ids),
  length(expected_samples),
  length(expected_groups),
  log2_pseudocount
))

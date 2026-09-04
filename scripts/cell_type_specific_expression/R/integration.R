scatter_inventory_columns <- c(
  "logical_name", "path", "sha256", "n_genes", "n_samples", "scale",
  "cell_group", "slug"
)

validate_scatter_nonempty_unique <- function(values, label) {
  if (!is.character(values) || anyNA(values) || any(!nzchar(values))) {
    stop(sprintf("Inventory %s must be non-empty", label), call. = FALSE)
  }
  if (anyDuplicated(values) > 0L) {
    stop(sprintf("Inventory %s must be unique", label), call. = FALSE)
  }
  invisible(TRUE)
}

validate_scatter_positive_dimension <- function(values, label) {
  valid <- is.numeric(values) && !anyNA(values) && all(is.finite(values)) &&
    all(values > 0) && all(values == as.integer(values))
  if (!valid) {
    stop(sprintf("Inventory %s must be positive integers", label), call. = FALSE)
  }
  invisible(TRUE)
}

validate_scatter_safe_slugs <- function(slugs) {
  canonical_slug_pattern <- "^[a-z0-9]+(_[a-z0-9]+)*$"
  if (!all(grepl(canonical_slug_pattern, slugs))) {
    stop("Inventory slugs must be safe filename tokens", call. = FALSE)
  }
  invisible(TRUE)
}

validate_output_prefix_token <- function(output_prefix) {
  safe_prefix_pattern <- "\\A[A-Za-z0-9][A-Za-z0-9._-]*\\z"
  valid <- is.character(output_prefix) && length(output_prefix) == 1L &&
    !is.na(output_prefix) &&
    grepl(safe_prefix_pattern, output_prefix, perl = TRUE)
  if (!valid) {
    stop(
      paste0(
        "output_prefix must be one safe basename token that starts with an ",
        "ASCII letter or number and contains only letters, numbers, dots, ",
        "underscores, or hyphens"
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

validate_scatter_inventory <- function(inventory) {
  if (!inherits(inventory, "data.frame") ||
      !identical(names(inventory), scatter_inventory_columns)) {
    stop("Inventory must contain the exact columns required for scatter", call. = FALSE)
  }
  if (nrow(inventory) == 0L) {
    stop("Inventory must contain at least one cell type", call. = FALSE)
  }

  validate_scatter_nonempty_unique(inventory$cell_group, "cell types")
  validate_scatter_nonempty_unique(inventory$slug, "slugs")
  validate_scatter_safe_slugs(inventory$slug)
  validate_scatter_nonempty_unique(inventory$path, "paths")
  validate_scatter_nonempty_unique(inventory$logical_name, "logical names")
  if (!is.character(inventory$sha256) || anyNA(inventory$sha256) ||
      any(!grepl("^[0-9a-f]{64}$", inventory$sha256))) {
    stop("Inventory sha256 values must be lowercase SHA-256 checksums", call. = FALSE)
  }
  if (!is.character(inventory$scale) || anyNA(inventory$scale) ||
      any(inventory$scale != "cpm")) {
    stop("Inventory scale must be cpm", call. = FALSE)
  }
  validate_scatter_positive_dimension(inventory$n_genes, "n_genes")
  validate_scatter_positive_dimension(inventory$n_samples, "n_samples")
  if (any(grepl("[/\\\\]", inventory$path)) ||
      any(basename(inventory$path) != inventory$path)) {
    stop("Inventory paths must be stable basenames", call. = FALSE)
  }
  invisible(TRUE)
}

validate_scatter_bed_paths <- function(bed_paths, expected_paths) {
  if (!is.character(bed_paths) || anyNA(bed_paths) || any(!nzchar(bed_paths))) {
    stop("BED paths must be non-empty character paths", call. = FALSE)
  }
  if (length(bed_paths) != length(expected_paths)) {
    stop("BED path count must match inventory row count", call. = FALSE)
  }
  bed_basenames <- basename(bed_paths)
  if (anyDuplicated(bed_basenames) > 0L ||
      !setequal(bed_basenames, expected_paths)) {
    stop("BED basenames must match inventory paths", call. = FALSE)
  }
  bed_paths[match(expected_paths, bed_basenames)]
}

prepare_scatter_contract <- function(inventory, bed_paths, output_prefix) {
  validate_scatter_inventory(inventory)
  bed_paths <- validate_scatter_bed_paths(bed_paths, inventory$path)
  validate_output_prefix_token(output_prefix)

  tibble::tibble(
    cell_type = inventory$cell_group,
    cell_type_slug = inventory$slug,
    expression_bed = bed_paths,
    output_prefix = paste(output_prefix, inventory$slug, sep = ".")
  )
}

build_cell_type_qtl_manifest <- function(
    cell_types,
    cell_type_slugs,
    int_beds,
    scaled_beds,
    int_pcs,
    int_pcs_all,
    scaled_pcs,
    scaled_pcs_all,
    int_covariates,
    scaled_covariates,
    int_outliers,
    scaled_outliers,
    source_beds = NULL,
    source_bed_slugs = NULL,
    filtered_beds = NULL,
    filtered_bed_slugs = NULL,
    negative_summary = NULL,
    gene_comparison = NULL,
    filter_metrics = NULL) {
  file_lists <- list(
    int_bed = int_beds,
    scaled_bed = scaled_beds,
    int_phenotype_pcs = int_pcs,
    int_phenotype_pcs_all = int_pcs_all,
    scaled_phenotype_pcs = scaled_pcs,
    scaled_phenotype_pcs_all = scaled_pcs_all,
    int_merged_covariates = int_covariates,
    scaled_merged_covariates = scaled_covariates,
    int_connectivity_outliers = int_outliers,
    scaled_connectivity_outliers = scaled_outliers
  )
  expected_length <- length(cell_types)
  if (expected_length == 0L || length(cell_type_slugs) != expected_length ||
      any(lengths(file_lists) != expected_length)) {
    stop(
      "Cell-type manifest arrays must have one equal nonzero length",
      call. = FALSE
    )
  }
  if (!is.character(cell_types) || anyNA(cell_types) || any(!nzchar(cell_types)) ||
      anyDuplicated(cell_types) > 0L || !is.character(cell_type_slugs) ||
      anyNA(cell_type_slugs) || any(!nzchar(cell_type_slugs)) ||
      anyDuplicated(cell_type_slugs) > 0L) {
    stop("Cell types and slugs must be nonempty and unique", call. = FALSE)
  }
  if (any(grepl("[[:cntrl:]]", c(cell_types, cell_type_slugs)))) {
    stop("Cell identities must not contain control characters", call. = FALSE)
  }
  if (any(!grepl("^[a-z0-9]+(_[a-z0-9]+)*$", cell_type_slugs))) {
    stop("Cell-type IDs must use safe lowercase slugs", call. = FALSE)
  }
  paths <- unlist(file_lists, use.names = FALSE)
  if (!all(vapply(file_lists, is.character, logical(1))) ||
      anyNA(paths) || any(!nzchar(paths)) || any(grepl("[[:cntrl:]]", paths))) {
    stop("Manifest paths must be nonempty strings without control characters", call. = FALSE)
  }
  # Paths are metadata supplied by completed upstream calls. Do not open them:
  # cloud objects and local runner host paths are not localized into this task.
  cloud_paths <- grepl("^gs://[^/[:space:]]+/[^[:cntrl:]]+$", paths) &
    !endsWith(paths, "/")
  local_paths <- startsWith(paths, "/")
  if (any(!cloud_paths & !local_paths)) {
    stop("Manifest paths must be full gs:// object URLs or absolute local paths", call. = FALSE)
  }
  duplicate_paths <- vapply(
    file_lists,
    function(paths) anyDuplicated(paths) > 0L,
    logical(1)
  )
  if (any(duplicate_paths)) {
    stop(
      "Manifest paths must be unique within each output category",
      call. = FALSE
    )
  }

  filter_metadata <- list(
    source_beds = source_beds,
    source_bed_slugs = source_bed_slugs,
    filtered_beds = filtered_beds,
    filtered_bed_slugs = filtered_bed_slugs,
    negative_summary = negative_summary,
    gene_comparison = gene_comparison,
    filter_metrics = filter_metrics
  )
  metadata_provided <- !vapply(filter_metadata, is.null, logical(1))
  if (any(metadata_provided) && !all(metadata_provided)) {
    stop("Filter manifest metadata must be supplied together", call. = FALSE)
  }
  if (all(metadata_provided)) {
    for (name in c("source_bed_slugs", "filtered_bed_slugs")) {
      slugs <- filter_metadata[[name]]
      if (!is.character(slugs) || length(slugs) != expected_length || anyNA(slugs) ||
          any(!nzchar(slugs)) || anyDuplicated(slugs) > 0L ||
          !setequal(slugs, cell_type_slugs)) {
        stop("Filter BED slugs must uniquely match manifest cell-type slugs", call. = FALSE)
      }
    }
    source_beds <- source_beds[match(cell_type_slugs, source_bed_slugs)]
    filtered_beds <- filtered_beds[match(cell_type_slugs, filtered_bed_slugs)]
    new_paths <- c(source_beds, filtered_beds, negative_summary, gene_comparison, filter_metrics)
    new_cloud_paths <- grepl("^gs://[^/[:space:]]+/[^[:cntrl:]]+$", new_paths) &
      !endsWith(new_paths, "/")
    new_local_paths <- startsWith(new_paths, "/")
    if (length(source_beds) != expected_length || length(filtered_beds) != expected_length ||
        anyNA(new_paths) || any(!nzchar(new_paths)) ||
        any(!new_cloud_paths & !new_local_paths)) {
      stop("Filter manifest paths must be full gs:// object URLs or absolute local paths",
           call. = FALSE)
    }
  } else {
    source_beds <- rep(NA_character_, expected_length)
    filtered_beds <- rep(NA_character_, expected_length)
    negative_summary <- NA_character_
    gene_comparison <- NA_character_
    filter_metrics <- NA_character_
  }

  file_columns <- tibble::as_tibble(file_lists)
  tibble::tibble(
    `entity:cell_type_id` = cell_type_slugs,
    cell_type = cell_types,
    cell_type_slug = cell_type_slugs
  ) |>
    dplyr::bind_cols(file_columns) |>
    dplyr::mutate(
      source_cpm_bed = source_beds,
      filtered_cpm_bed = filtered_beds,
      negative_expression_summary = negative_summary,
      reference_gene_comparison = gene_comparison,
      reference_filter_metrics = filter_metrics
    )
}

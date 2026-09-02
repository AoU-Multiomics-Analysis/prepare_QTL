scatter_inventory_columns <- c(
  "logical_name", "path", "n_genes", "n_samples", "scale", "cell_group", "slug"
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
  if (!is.character(inventory$scale) || anyNA(inventory$scale) ||
      any(inventory$scale != "log2_cpm")) {
    stop("Inventory scale must be log2_cpm", call. = FALSE)
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
  if (!identical(basename(bed_paths), expected_paths)) {
    stop("BED basenames must match inventory paths in order", call. = FALSE)
  }
  invisible(TRUE)
}

prepare_scatter_contract <- function(inventory, bed_paths, output_prefix) {
  validate_scatter_inventory(inventory)
  validate_scatter_bed_paths(bed_paths, inventory$path)
  if (!is.character(output_prefix) || length(output_prefix) != 1L ||
      is.na(output_prefix) || !nzchar(output_prefix)) {
    stop("output_prefix must be one non-empty string", call. = FALSE)
  }

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
    scaled_outliers) {
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
  if (!all(vapply(file_lists, is.character, logical(1))) ||
      any(!file.exists(unlist(file_lists, use.names = FALSE)))) {
    stop("Every cell-type manifest file must exist", call. = FALSE)
  }
  duplicate_basenames <- vapply(
    file_lists,
    function(paths) anyDuplicated(basename(paths)) > 0L,
    logical(1)
  )
  if (any(duplicate_basenames)) {
    stop(
      "Manifest file basenames must be unique within each output category",
      call. = FALSE
    )
  }

  file_columns <- file_lists |>
    purrr::map(basename) |>
    tibble::as_tibble()
  tibble::tibble(
    cell_type = cell_types,
    cell_type_slug = cell_type_slugs
  ) |>
    dplyr::bind_cols(file_columns)
}

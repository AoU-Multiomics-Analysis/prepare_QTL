read_cell_type_mapping <- function(path = NULL) {
  if (is.null(path) || !nzchar(path)) return(NULL)
  if (grepl("^(gs|s3|https?)://", path) || !file.exists(path) || file.access(path, 4) != 0) {
    stop("Cell-type mapping localization error: expected a readable local file", call. = FALSE)
  }
  x <- readr::read_tsv(path, col_types = readr::cols(.default = readr::col_character()))
  if (!identical(names(x), c("source", "target")) || nrow(x) < 2L ||
      anyNA(x) || any(!nzchar(trimws(x$source))) || any(!nzchar(trimws(x$target))) ||
      anyDuplicated(x$source) || any(trimws(x$source) != x$source) ||
      any(trimws(x$target) != x$target)) {
    stop("Mapping requires source and target columns, unique sources, and nonempty labels", call. = FALSE)
  }
  x
}

validate_mapping_columns <- function(mapping, columns) {
  if (!setequal(mapping$source, columns)) {
    stop("Mapping sources must exactly match the reference/proportion columns", call. = FALSE)
  }
  invisible(TRUE)
}

combine_mapped_proportions <- function(proportions, mapping) {
  validate_combined_proportions(proportions)
  validate_mapping_columns(mapping, colnames(proportions))
  if (any(abs(rowSums(proportions) - 1) > 1e-6)) {
    stop("Proportion rows must sum to one", call. = FALSE)
  }
  targets <- unique(mapping$target)
  combined <- purrr::map(targets, function(target) {
    rowSums(proportions[, mapping$source[mapping$target == target], drop = FALSE])
  }) |> purrr::set_names(targets) |> tibble::as_tibble() |> as.matrix()
  rownames(combined) <- rownames(proportions)
  combined
}

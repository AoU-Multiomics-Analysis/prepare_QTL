reference_population_map <- c(
  NveB = "B cells", MemB = "B cells", CD4T = "CD4 T cells",
  CD8T = "CD8 T cells", NK = "NK cells", Mono = "Monocyte/myeloid",
  MonoNonClassical = "Monocyte/myeloid", Neut = "Neutrophils",
  Eo = "Eosinophils", myDC = "Dendritic cells",
  myDC123 = "Dendritic cells", pDC = "Dendritic cells"
)

reference_required_populations <- split(names(reference_population_map), reference_population_map)

reference_mapping_caveats <- c(
  `B cells` = "The model group includes plasma cells, which are absent from the reference.",
  `Monocyte/myeloid` = "The model group includes macrophage components; the reference contains monocytes.",
  `Dendritic cells` = "The comparison pools distinct sorted dendritic populations."
)

reference_gene_key <- function(gene_ids, unique = FALSE) {
  if (!is.character(gene_ids) || anyNA(gene_ids) || any(!nzchar(gene_ids))) {
    stop("Gene IDs must be non-empty", call. = FALSE)
  }
  keys <- sub("^(ENS[A-Z]*G[0-9]+)[.][0-9]+$", "\\1", gene_ids)
  if (unique && anyDuplicated(keys) > 0L) {
    stop("Gene IDs have a duplicate matching key after version removal", call. = FALSE)
  }
  keys
}

map_reference_samples <- function(sample_names) {
  if (!is.character(sample_names) || anyNA(sample_names) ||
      any(!nzchar(sample_names)) || anyDuplicated(sample_names) > 0L) {
    stop("Reference sample names must be non-empty and unique", call. = FALSE)
  }
  population <- sub("[.][0-9]+$", "", sample_names)
  unknown <- setdiff(unique(population), names(reference_population_map))
  if (length(unknown) > 0L) {
    stop(sprintf("Unknown reference sample prefix: %s", paste(unknown, collapse = ", ")), call. = FALSE)
  }
  present_lineages <- unique(unname(reference_population_map[population]))
  for (cell_type in present_lineages) {
    missing <- setdiff(reference_required_populations[[cell_type]], population)
    if (length(missing) > 0L) {
      stop(sprintf("Reference lineage '%s' is missing required population(s): %s",
                   cell_type, paste(missing, collapse = ", ")), call. = FALSE)
    }
  }
  tibble::tibble(sample_id = sample_names, population = population,
                 cell_type = unname(reference_population_map[population]))
}

prepare_reference_matrix <- function(counts) {
  if (!is.matrix(counts) || !is.numeric(counts) || is.null(rownames(counts)) ||
      is.null(colnames(counts))) stop("counts must be a named numeric matrix", call. = FALSE)
  reference_gene_key(rownames(counts), unique = TRUE)
  sample_map <- map_reference_samples(colnames(counts))
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) {
    stop("Reference counts must be finite, nonnegative integers", call. = FALSE)
  }
  if (any(colSums(counts) <= 0)) stop("Reference library totals must be positive", call. = FALSE)
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("The edgeR package is required", call. = FALSE)
  dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts))
  list(cpm = edgeR::cpm(dge, log = FALSE), sample_map = sample_map,
       samples = dplyr::mutate(sample_map,
         library_size = unname(dge$samples$lib.size),
         norm_factor = unname(dge$samples$norm.factors)))
}

summarize_reference_expression <- function(cpm, sample_map) {
  purrr::map_dfr(unique(sample_map$cell_type), function(cell_type) {
    columns <- sample_map$sample_id[sample_map$cell_type == cell_type]
    values <- cpm[, columns, drop = FALSE]
    logged <- log2(values + 1)
    tibble::tibble(
      gene_id = rownames(cpm), cell_type = cell_type, n_samples = ncol(values),
      mean_log2_cpm1 = unname(rowMeans(logged)),
      median_log2_cpm1 = unname(apply(logged, 1L, stats::median))
    )
  })
}

read_reference_counts <- function(path) {
  table <- readr::read_tsv(path, col_types = readr::cols(.default = readr::col_character()),
                           name_repair = "minimal", progress = FALSE, show_col_types = FALSE)
  if (ncol(table) < 2L) stop("Reference counts must contain genes and samples", call. = FALSE)
  first <- names(table)[[1L]]
  if (!(first %in% c("", "V1", "gene_id", "...1"))) {
    stop("Reference first column must be empty, V1, or gene_id", call. = FALSE)
  }
  names(table)[[1L]] <- "gene_id"
  if (anyDuplicated(names(table)) > 0L) stop("Reference columns must be unique", call. = FALSE)
  ids <- table$gene_id
  reference_gene_key(ids, unique = TRUE)
  sample_names <- names(table)[-1L]
  values <- lapply(table[-1L], function(column) {
    parsed <- suppressWarnings(readr::parse_double(column))
    if (anyNA(parsed)) stop("Reference counts contain missing or nonnumeric values", call. = FALSE)
    parsed
  })
  counts <- do.call(cbind, values)
  dimnames(counts) <- list(ids, sample_names)
  counts
}

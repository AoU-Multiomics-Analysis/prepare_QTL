validate_log2_pseudocount <- function(x) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x)) {
    stop("log2_pseudocount must be one finite numeric value", call. = FALSE)
  }
  if (x < 0) {
    stop("log2_pseudocount must be non-negative", call. = FALSE)
  }
  as.numeric(x)
}

validate_cpm_matrix <- function(cpm) {
  if (!is.matrix(cpm) || !is.numeric(cpm)) {
    stop("cpm must be a numeric matrix", call. = FALSE)
  }
  if (nrow(cpm) == 0L || ncol(cpm) == 0L) {
    stop("cpm must contain at least one gene and one sample", call. = FALSE)
  }
  if (is.null(rownames(cpm)) || anyNA(rownames(cpm))) {
    stop("CPM gene identifiers must be non-missing", call. = FALSE)
  }

  gene_ids <- trimws(rownames(cpm))
  if (any(!nzchar(gene_ids))) {
    stop("CPM gene identifiers must be non-empty", call. = FALSE)
  }
  if (anyDuplicated(gene_ids) > 0L) {
    stop("CPM gene identifiers must be unique after trimming", call. = FALSE)
  }
  if (is.null(colnames(cpm)) || anyNA(colnames(cpm)) ||
      any(!nzchar(colnames(cpm)))) {
    stop("CPM sample identifiers must be non-missing and non-empty", call. = FALSE)
  }
  if (anyDuplicated(colnames(cpm)) > 0L) {
    stop("CPM sample identifiers must be unique", call. = FALSE)
  }
  if (any(!is.finite(cpm))) {
    stop("CPM values must be finite", call. = FALSE)
  }
  if (any(cpm < 0)) {
    stop("CPM values must be non-negative", call. = FALSE)
  }

  rownames(cpm) <- gene_ids
  cpm
}

extract_gtf_attribute <- function(attributes, key) {
  if (!is.character(key) || length(key) != 1L || is.na(key) || !nzchar(key)) {
    stop("key must be one non-missing, non-empty character value", call. = FALSE)
  }

  pattern <- paste0("(?:^|;)[[:space:]]*", key,
    "[[:space:]]+\\\"([^\\\"]+)\\\"")
  stringr::str_match(attributes, pattern)[, 2L]
}

validate_gtf_gene_annotation <- function(annotation, require_gene_name = TRUE) {
  if (!inherits(annotation, "data.frame")) {
    stop("GTF annotation must be a data frame", call. = FALSE)
  }
  if (!is.logical(require_gene_name) || length(require_gene_name) != 1L ||
      is.na(require_gene_name)) {
    stop("require_gene_name must be true or false", call. = FALSE)
  }

  required_columns <- c("gene_id", "gene_name", "gene_type")
  missing_columns <- setdiff(required_columns, names(annotation))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "GTF annotation is missing required columns: %s",
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  annotation <- tibble::as_tibble(annotation) |>
    dplyr::transmute(
      gene_id = trimws(as.character(.data$gene_id)),
      gene_name = trimws(as.character(.data$gene_name)),
      gene_type = trimws(as.character(.data$gene_type))
    )
  if (nrow(annotation) == 0L) {
    stop("GTF annotation must contain at least one gene record", call. = FALSE)
  }
  if (anyNA(annotation$gene_id) || any(!nzchar(annotation$gene_id))) {
    stop("GTF annotation gene_id values must be non-missing and non-empty", call. = FALSE)
  }
  if (anyDuplicated(annotation$gene_id) > 0L) {
    stop("GTF annotation gene_id values must be unique", call. = FALSE)
  }
  if (require_gene_name &&
      !any(!is.na(annotation$gene_name) & nzchar(annotation$gene_name))) {
    stop("GTF annotation must contain at least one usable gene_name", call. = FALSE)
  }

  annotation
}

validate_gene_types <- function(gene_types) {
  if (!is.character(gene_types) || length(gene_types) == 0L) {
    stop("gene_type must contain at least one value", call. = FALSE)
  }
  if (anyNA(gene_types)) {
    stop("gene_type values must be non-missing", call. = FALSE)
  }

  gene_types <- trimws(gene_types)
  if (any(!nzchar(gene_types))) {
    stop("gene_type values must be non-empty", call. = FALSE)
  }
  if (anyDuplicated(gene_types) > 0L) {
    stop("gene_type values must be unique", call. = FALSE)
  }
  gene_types
}

parse_gene_types <- function(value) {
  if (!is.character(value) || length(value) != 1L || is.na(value)) {
    stop("gene_types must be one comma-separated character value", call. = FALSE)
  }
  values <- strsplit(value, ",", fixed = TRUE)[[1L]]
  if (endsWith(value, ",")) {
    values <- c(values, "")
  }
  validate_gene_types(values)
}

filter_expression_by_gene_types <- function(expression, annotation, gene_types) {
  if (!is.list(expression) || !all(c("coordinates", "cpm") %in% names(expression))) {
    stop("expression must contain coordinates and cpm", call. = FALSE)
  }

  gene_types <- validate_gene_types(gene_types)
  annotation <- validate_gtf_gene_annotation(annotation, require_gene_name = FALSE)
  cpm <- validate_cpm_matrix(expression$cpm)
  coordinates <- tibble::as_tibble(expression$coordinates)
  if (!identical(names(coordinates), expression_bed_columns()) ||
      !identical(rownames(cpm), coordinates$gene_id)) {
    stop("Expression genes and BED coordinates must match exactly", call. = FALSE)
  }

  annotation_index <- match(coordinates$gene_id, annotation$gene_id)
  gene_names <- annotation$gene_name[annotation_index]
  annotated_gene_types <- annotation$gene_type[annotation_index]
  has_gene_type <- !is.na(annotated_gene_types) & nzchar(annotated_gene_types)
  retained <- !is.na(annotation_index) & has_gene_type &
    annotated_gene_types %in% gene_types
  report <- tibble::tibble(
    gene_id = coordinates$gene_id,
    gene_name = gene_names,
    gene_type = annotated_gene_types,
    retained = retained,
    filter_reason = dplyr::case_when(
      is.na(annotation_index) ~ "missing_gtf_gene_id",
      !has_gene_type ~ "missing_gene_type",
      !retained ~ "gene_type_not_selected",
      TRUE ~ "retained"
    )
  )
  if (!any(retained)) {
    stop("No expression genes remain after gene-type filtering", call. = FALSE)
  }

  filtered_expression <- list(
    coordinates = coordinates[retained, , drop = FALSE],
    cpm = cpm[retained, , drop = FALSE]
  )
  if ("log2_pseudocount" %in% names(expression)) {
    filtered_expression$log2_pseudocount <- expression$log2_pseudocount
  }
  list(expression = filtered_expression, report = report)
}

read_gtf_gene_annotation <- function(
    path,
    chunk_size = 100000L,
    require_gene_name = TRUE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("GTF path must be one non-missing, non-empty character value", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop(sprintf("GTF file does not exist: %s", path), call. = FALSE)
  }
  if (!is.numeric(chunk_size) || length(chunk_size) != 1L || is.na(chunk_size) ||
      !is.finite(chunk_size) || chunk_size < 1L || chunk_size != as.integer(chunk_size)) {
    stop("chunk_size must be a positive integer", call. = FALSE)
  }

  connection <- if (grepl("[.]gz$", path, ignore.case = TRUE)) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }
  on.exit(close(connection), add = TRUE)

  gene_rows <- list()
  chunk_index <- 0L
  repeat {
    lines <- readLines(connection, n = as.integer(chunk_size), warn = FALSE)
    if (length(lines) == 0L) {
      break
    }

    records <- lines[!startsWith(trimws(lines), "#")]
    if (length(records) == 0L) {
      next
    }
    fields_by_record <- strsplit(records, "\t", fixed = TRUE)
    field_counts <- lengths(fields_by_record)
    if (any(field_counts != 9L)) {
      stop("Each non-comment GTF record must contain nine tab-separated fields", call. = FALSE)
    }
    fields <- do.call(rbind, fields_by_record)
    is_gene <- fields[, 3L] == "gene"
    if (!any(is_gene)) {
      next
    }

    gene_fields <- fields[is_gene, , drop = FALSE]
    chunk_index <- chunk_index + 1L
    gene_rows[[chunk_index]] <- tibble::tibble(
      gene_id = extract_gtf_attribute(gene_fields[, 9L], "gene_id"),
      gene_name = extract_gtf_attribute(gene_fields[, 9L], "gene_name"),
      gene_type = dplyr::coalesce(
        extract_gtf_attribute(gene_fields[, 9L], "gene_type"),
        extract_gtf_attribute(gene_fields[, 9L], "gene_biotype")
      )
    )
  }

  if (length(gene_rows) == 0L) {
    stop("GTF annotation must contain at least one gene record", call. = FALSE)
  }
  validate_gtf_gene_annotation(
    dplyr::bind_rows(gene_rows),
    require_gene_name = require_gene_name
  )
}

collapse_cpm_to_gene_names <- function(cpm, annotation) {
  cpm <- validate_cpm_matrix(cpm)
  annotation <- validate_gtf_gene_annotation(annotation)
  annotation_index <- match(rownames(cpm), annotation$gene_id)
  mapped_gene_names <- annotation$gene_name[annotation_index]
  usable <- !is.na(mapped_gene_names) & nzchar(mapped_gene_names)
  if (!any(usable)) {
    return(matrix(
      numeric(),
      nrow = 0L,
      ncol = ncol(cpm),
      dimnames = list(character(), colnames(cpm))
    ))
  }

  mapped_gene_names <- mapped_gene_names[usable]
  gene_name_order <- unique(mapped_gene_names)
  group_index <- match(mapped_gene_names, gene_name_order)
  output <- rowsum(
    cpm[usable, , drop = FALSE],
    group = group_index,
    reorder = FALSE
  )
  storage.mode(output) <- "double"
  rownames(output) <- gene_name_order
  colnames(output) <- colnames(cpm)
  output
}

make_cpm_mapping_report <- function(cpm, annotation) {
  matching_indices <- match(rownames(cpm), annotation$gene_id)
  gene_names <- annotation$gene_name[matching_indices]
  usable_gene_names <- !is.na(gene_names) & nzchar(gene_names)
  duplicate_gene_names <- usable_gene_names &
    (duplicated(gene_names) | duplicated(gene_names, fromLast = TRUE))

  tibble::tibble(
    gene_id = rownames(cpm),
    gene_name = gene_names,
    mapping_action = dplyr::case_when(
      is.na(matching_indices) ~ "missing_gtf_gene_id",
      !usable_gene_names ~ "missing_gene_name",
      duplicate_gene_names ~ "duplicate_gene_name_aggregated_for_dtangle",
      TRUE ~ "mapped"
    )
  )
}

make_tca_expression <- function(expression, log2_pseudocount = 0) {
  if (!is.list(expression) || !all(c("coordinates", "cpm") %in% names(expression))) {
    stop("expression must contain coordinates and cpm", call. = FALSE)
  }
  log2_pseudocount <- validate_log2_pseudocount(log2_pseudocount)
  cpm <- validate_cpm_matrix(expression$cpm)
  if (!identical(rownames(cpm), expression$coordinates$gene_id)) {
    stop("Expression genes and BED coordinates must match exactly", call. = FALSE)
  }
  if (log2_pseudocount == 0 && any(cpm == 0)) {
    stop(
      "strictly positive CPM values are required when log2_pseudocount is zero",
      call. = FALSE
    )
  }
  log2(cpm + log2_pseudocount)
}

make_dtangle_expression <- function(expression, annotation, log2_pseudocount = 0) {
  if (!is.list(expression) || !all(c("coordinates", "cpm") %in% names(expression))) {
    stop("expression must contain coordinates and cpm", call. = FALSE)
  }
  log2_pseudocount <- validate_log2_pseudocount(log2_pseudocount)
  cpm <- validate_cpm_matrix(expression$cpm)
  annotation <- validate_gtf_gene_annotation(annotation)
  mapping_report <- make_cpm_mapping_report(cpm, annotation)
  symbol_cpm <- collapse_cpm_to_gene_names(cpm, annotation)
  if (log2_pseudocount == 0 && any(cpm == 0)) {
    stop(
      "strictly positive CPM values are required when log2_pseudocount is zero",
      call. = FALSE
    )
  }
  if (nrow(symbol_cpm) == 0L) {
    stop("No CPM genes map to a usable gene symbol", call. = FALSE)
  }
  list(
    log_expression = log2(symbol_cpm + log2_pseudocount),
    mapping_report = mapping_report,
    log2_pseudocount = log2_pseudocount
  )
}

make_excluded_genes <- function(mapping_report) {
  mapping_report |>
    dplyr::filter(!(.data$mapping_action %in% c(
      "mapped", "duplicate_gene_name_aggregated_for_dtangle"
    ))) |>
    dplyr::transmute(
      gene_id = .data$gene_id,
      gene_name = .data$gene_name,
      reason = .data$mapping_action
    )
}

# Summarize exported TCA values as stored, including zeros and negative estimates.
# Only a small block of expression rows is held in memory at a time.
summarize_cell_type_beds <- function(inventory, bed_paths, output_path, chunk_size = 256L) {
  validate_scatter_inventory(inventory)
  bed_paths <- validate_scatter_bed_paths(bed_paths, inventory$path)
  chunk_size <- validate_tensor_positive_integer(chunk_size, "chunk_size")
  if (!all(file.exists(bed_paths))) {
    stop("Every cell-type BED must exist", call. = FALSE)
  }
  if (normalizePath(output_path, mustWork = FALSE) %in% normalizePath(bed_paths)) {
    stop("Summary output must not overwrite an input BED", call. = FALSE)
  }
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  output_connection <- gzfile(output_path, open = "wb")
  on.exit(close(output_connection), add = TRUE)
  write_header <- TRUE
  total_rows <- 0L
  expected_samples <- NULL
  expected_coordinates <- NULL

  purrr::walk(seq_len(nrow(inventory)), function(index) {
    cell_type <- inventory$cell_group[[index]]
    bed_path <- bed_paths[[index]]
    header <- names(readr::read_tsv(
      bed_path, n_max = 0, col_types = readr::cols(.default = readr::col_character()),
      name_repair = "minimal", progress = FALSE, show_col_types = FALSE
    ))
    coordinate_columns <- expression_bed_columns()
    if (length(header) <= 4L || !identical(header[1:4], coordinate_columns)) {
      stop("BED must start with #chr, start, end, and gene_id", call. = FALSE)
    }
    if (anyNA(header) || any(!nzchar(header)) || anyDuplicated(header) > 0L) {
      stop("BED column and sample names must be non-empty and unique", call. = FALSE)
    }
    sample_ids <- header[-(1:4)]
    n_samples <- length(sample_ids)
    if (n_samples != inventory$n_samples[[index]]) {
      stop("BED sample count does not match inventory", call. = FALSE)
    }
    if (is.null(expected_samples)) expected_samples <<- sample_ids
    assert_identical_ids(expected_samples, sample_ids, "BED sample")
    coordinates <- list()
    seen_genes <- character()
    gene_count <- 0L
    message(sprintf(
      "stage=summarize_cell_type_beds cell_type=%s samples=%d start_time=%s",
      cell_type, n_samples, tensor_utc_time()
    ))

    callback <- readr::SideEffectChunkCallback$new(function(chunk, position) {
      if (nrow(chunk) == 0L) return(invisible(NULL))
      if (nrow(readr::problems(chunk)) > 0L) {
        stop("BED parse error: invalid or missing fields", call. = FALSE)
      }
      coords <- chunk |> dplyr::select(dplyr::all_of(coordinate_columns))
      if (anyNA(coords) || any(!nzchar(trimws(coords[["#chr"]]))) ||
          any(!nzchar(trimws(coords$gene_id))) ||
          any(coords$start < 0L) || any(coords$end <= coords$start)) {
        stop("BED coordinates and gene IDs must be valid and non-missing", call. = FALSE)
      }
      if (anyDuplicated(coords$gene_id) > 0L || any(coords$gene_id %in% seen_genes)) {
        stop("BED gene IDs must be unique", call. = FALSE)
      }
      seen_genes <<- c(seen_genes, coords$gene_id)
      coordinates[[length(coordinates) + 1L]] <<- coords
      values <- chunk |> dplyr::select(dplyr::all_of(sample_ids)) |> as.matrix()
      if (any(!is.finite(values))) {
        stop("BED expression values must be finite and non-missing", call. = FALSE)
      }
      statistics <- purrr::map_dfr(seq_len(nrow(values)), function(row_index) {
        row_values <- values[row_index, ]
        quartiles <- stats::quantile(row_values, c(0.25, 0.5, 0.75), type = 7, names = FALSE)
        standard_deviation <- stats::sd(row_values)
        tibble::tibble(
          mean_cpm = mean(row_values), median_cpm = quartiles[[2]],
          sd_cpm = standard_deviation, se_mean_cpm = standard_deviation / sqrt(n_samples),
          q1_cpm = quartiles[[1]], q3_cpm = quartiles[[3]],
          iqr_cpm = quartiles[[3]] - quartiles[[1]]
        )
      })
      result <- coords |>
        dplyr::mutate(cell_type = cell_type, n_samples = n_samples, scale = "cpm", .before = 1) |>
        dplyr::bind_cols(statistics)
      readr::write_tsv(result, output_connection, col_names = write_header, na = "NA")
      write_header <<- FALSE
      gene_count <<- gene_count + nrow(chunk)
      total_rows <<- total_rows + nrow(chunk)
      invisible(NULL)
    })
    readr::read_tsv_chunked(
      bed_path, callback = callback, chunk_size = chunk_size,
      col_types = readr::cols(
        `#chr` = readr::col_character(), start = readr::col_integer(),
        end = readr::col_integer(), gene_id = readr::col_character(),
        .default = readr::col_double()
      ),
      progress = FALSE
    )
    if (gene_count != inventory$n_genes[[index]]) {
      stop("BED gene count does not match inventory", call. = FALSE)
    }
    bed_coordinates <- dplyr::bind_rows(coordinates)
    if (is.null(expected_coordinates)) expected_coordinates <<- bed_coordinates
    if (!identical(bed_coordinates, expected_coordinates)) {
      stop("Cell-type BED gene coordinates and order must match", call. = FALSE)
    }
    message(sprintf(
      "stage=summarize_cell_type_beds cell_type=%s genes=%d completion_time=%s",
      cell_type, gene_count, tensor_utc_time()
    ))
  })
  message(sprintf(
    "stage=summarize_cell_type_beds dimensions=rows:%d,cell_types:%d outputs=%s completion_time=%s",
    total_rows, nrow(inventory), output_path, tensor_utc_time()
  ))
  invisible(output_path)
}

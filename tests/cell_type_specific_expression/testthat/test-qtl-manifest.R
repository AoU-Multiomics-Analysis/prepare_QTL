qtl_manifest_columns <- c(
  "cell_type", "cell_type_slug", "int_bed", "scaled_bed",
  "int_phenotype_pcs", "int_phenotype_pcs_all",
  "scaled_phenotype_pcs", "scaled_phenotype_pcs_all",
  "int_merged_covariates", "scaled_merged_covariates",
  "int_connectivity_outliers", "scaled_connectivity_outliers"
)

make_qtl_manifest_inputs <- function(root = tempfile("qtl-manifest-inputs-")) {
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  cell_types <- c("Monocytes", "CD4 T cells")
  cell_type_slugs <- c("monocytes", "cd4_t_cells")
  argument_names <- c(
    "int_beds", "scaled_beds", "int_pcs", "int_pcs_all",
    "scaled_pcs", "scaled_pcs_all", "int_covariates",
    "scaled_covariates", "int_outliers", "scaled_outliers"
  )
  extensions <- c(
    "int.bed.gz", "scaled.bed.gz", "int.pcs.tsv", "int.pcs.all.tsv",
    "scaled.pcs.tsv", "scaled.pcs.all.tsv", "int.covariates.tsv",
    "scaled.covariates.tsv", "int.outliers.tsv", "scaled.outliers.tsv"
  )

  file_lists <- purrr::map2(argument_names, extensions, function(argument, extension) {
    category_directory <- file.path(root, argument)
    dir.create(category_directory, recursive = TRUE, showWarnings = FALSE)
    paths <- file.path(category_directory, paste(cell_type_slugs, extension, sep = "."))
    purrr::walk(paths, ~ writeLines("fixture", .x))
    paths
  })
  names(file_lists) <- argument_names

  c(
    list(cell_types = cell_types, cell_type_slugs = cell_type_slugs),
    file_lists
  )
}

testthat::test_that("QTL manifest has the exact schema and preserves scatter order", {
  inputs <- make_qtl_manifest_inputs()

  manifest <- do.call(build_cell_type_qtl_manifest, inputs)

  testthat::expect_s3_class(manifest, "tbl_df")
  testthat::expect_identical(names(manifest), qtl_manifest_columns)
  testthat::expect_identical(manifest$cell_type, c("Monocytes", "CD4 T cells"))
  testthat::expect_identical(manifest$cell_type_slug, c("monocytes", "cd4_t_cells"))
  testthat::expect_identical(
    manifest$int_bed,
    c("monocytes.int.bed.gz", "cd4_t_cells.int.bed.gz")
  )
  file_values <- unlist(manifest[-c(1L, 2L)], use.names = FALSE)
  testthat::expect_false(any(grepl("[/\\\\]", file_values)))
})

testthat::test_that("QTL manifest requires equal nonzero array lengths", {
  inputs <- make_qtl_manifest_inputs()
  inputs$scaled_outliers <- inputs$scaled_outliers[-1L]
  testthat::expect_error(
    do.call(build_cell_type_qtl_manifest, inputs),
    "one equal nonzero length"
  )

  empty_inputs <- purrr::map(inputs, ~ .x[0])
  testthat::expect_error(
    do.call(build_cell_type_qtl_manifest, empty_inputs),
    "one equal nonzero length"
  )
})

testthat::test_that("QTL manifest requires nonempty unique names and slugs", {
  inputs <- make_qtl_manifest_inputs()
  empty_name_inputs <- inputs
  empty_name_inputs$cell_types[[1L]] <- ""
  testthat::expect_error(
    do.call(build_cell_type_qtl_manifest, empty_name_inputs),
    "nonempty and unique"
  )

  duplicate_slug_inputs <- inputs
  duplicate_slug_inputs$cell_type_slugs[[2L]] <-
    duplicate_slug_inputs$cell_type_slugs[[1L]]
  testthat::expect_error(
    do.call(build_cell_type_qtl_manifest, duplicate_slug_inputs),
    "nonempty and unique"
  )

  empty_slug_inputs <- inputs
  empty_slug_inputs$cell_type_slugs[[1L]] <- ""
  testthat::expect_error(
    do.call(build_cell_type_qtl_manifest, empty_slug_inputs),
    "nonempty and unique"
  )

  duplicate_name_inputs <- inputs
  duplicate_name_inputs$cell_types[[2L]] <- duplicate_name_inputs$cell_types[[1L]]
  testthat::expect_error(
    do.call(build_cell_type_qtl_manifest, duplicate_name_inputs),
    "nonempty and unique"
  )
})

testthat::test_that("QTL manifest requires every input file to exist", {
  inputs <- make_qtl_manifest_inputs()
  unlink(inputs$int_pcs[[2L]])

  testthat::expect_error(
    do.call(build_cell_type_qtl_manifest, inputs),
    "must exist"
  )
})

testthat::test_that("QTL manifest rejects duplicate basenames within a category", {
  inputs <- make_qtl_manifest_inputs()
  duplicate_directory <- file.path(dirname(inputs$int_beds[[1L]]), "localized")
  dir.create(duplicate_directory)
  duplicate_path <- file.path(duplicate_directory, basename(inputs$int_beds[[1L]]))
  writeLines("fixture", duplicate_path)
  inputs$int_beds[[2L]] <- duplicate_path

  testthat::expect_error(
    do.call(build_cell_type_qtl_manifest, inputs),
    "basenames must be unique"
  )
})

testthat::test_that("QTL manifest CLI writes stable basenames in scatter order", {
  root <- testthat::test_path("..", "..", "..")
  script <- file.path(
    root,
    "scripts",
    "cell_type_specific_expression",
    "build_qtl_manifest.R"
  )
  temporary_directory <- withr::local_tempdir()
  inputs <- make_qtl_manifest_inputs(file.path(temporary_directory, "localized inputs"))
  argument_files <- purrr::imap(inputs, function(values, argument_name) {
    path <- file.path(temporary_directory, paste0(argument_name, ".txt"))
    writeLines(values, path)
    path
  })
  output_path <- file.path(temporary_directory, "outputs", "cell_type_qtl_manifest.tsv")
  cli_names <- c(
    cell_types = "--cell-types",
    cell_type_slugs = "--cell-type-slugs",
    int_beds = "--int-beds",
    scaled_beds = "--scaled-beds",
    int_pcs = "--int-pcs",
    int_pcs_all = "--int-pcs-all",
    scaled_pcs = "--scaled-pcs",
    scaled_pcs_all = "--scaled-pcs-all",
    int_covariates = "--int-covariates",
    scaled_covariates = "--scaled-covariates",
    int_outliers = "--int-outliers",
    scaled_outliers = "--scaled-outliers"
  )
  cli_arguments <- unlist(purrr::imap(
    cli_names,
    ~ c(.x, argument_files[[.y]])
  ), use.names = FALSE)

  status <- suppressWarnings(system2(
    "Rscript",
    c(script, cli_arguments, "--output", output_path)
  ))

  testthat::expect_identical(status, 0L)
  if (status != 0L) {
    return(invisible())
  }
  manifest <- readr::read_tsv(
    output_path,
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  testthat::expect_identical(names(manifest), qtl_manifest_columns)
  testthat::expect_identical(manifest$cell_type, inputs$cell_types)
  testthat::expect_identical(
    manifest$int_bed,
    basename(inputs$int_beds)
  )
  testthat::expect_false(any(grepl("[/\\\\]", unlist(manifest[-c(1L, 2L)]))))
})

testthat::test_that("QTL manifest WDL task accepts and forwards all aligned arrays", {
  path <- testthat::test_path(
    "..", "..", "..", "workflows", "cell_type_specific_expression",
    "tasks", "integration.wdl"
  )
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  testthat::expect_match(text, "task BuildQtlManifest", fixed = TRUE)
  expected_inputs <- c(
    "Array[String] cell_types", "Array[String] cell_type_slugs",
    "Array[File] int_beds", "Array[File] scaled_beds",
    "Array[File] int_phenotype_pcs", "Array[File] int_phenotype_pcs_all",
    "Array[File] scaled_phenotype_pcs", "Array[File] scaled_phenotype_pcs_all",
    "Array[File] int_merged_covariates", "Array[File] scaled_merged_covariates",
    "Array[File] int_connectivity_outliers",
    "Array[File] scaled_connectivity_outliers"
  )
  purrr::walk(expected_inputs, ~ testthat::expect_match(text, .x, fixed = TRUE))
  testthat::expect_match(text, "validated_cell_count", fixed = TRUE)
  testthat::expect_match(text, "manifest_path", fixed = TRUE)
  testthat::expect_match(text, "completion_time", fixed = TRUE)
  testthat::expect_match(
    text,
    'File manifest = "outputs/cell_type_qtl_manifest.tsv"',
    fixed = TRUE
  )
  argument_files <- c(
    "inputs/cell_types.txt", "inputs/cell_type_slugs.txt", "inputs/int_beds.txt",
    "inputs/scaled_beds.txt", "inputs/int_phenotype_pcs.txt",
    "inputs/int_phenotype_pcs_all.txt", "inputs/scaled_phenotype_pcs.txt",
    "inputs/scaled_phenotype_pcs_all.txt", "inputs/int_merged_covariates.txt",
    "inputs/scaled_merged_covariates.txt", "inputs/int_connectivity_outliers.txt",
    "inputs/scaled_connectivity_outliers.txt"
  )
  purrr::walk(argument_files, function(argument_file) {
    testthat::expect_match(text, paste0("cat > ", argument_file), fixed = TRUE)
  })
  manifest_task <- stringr::str_extract(
    text,
    "task BuildQtlManifest \\{[\\s\\S]*\\z"
  )
  testthat::expect_identical(
    stringr::str_count(manifest_task, stringr::fixed("~{sep='\\n' ")),
    length(argument_files)
  )
})

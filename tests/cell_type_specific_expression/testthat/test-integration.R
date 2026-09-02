make_scatter_inventory <- function() {
  tibble::tibble(
    logical_name = c("cd4_t_cells_expression", "monocytes_expression"),
    path = c("cd4_t_cells.bed.gz", "monocytes.bed.gz"),
    n_genes = c(120L, 120L),
    n_samples = c(80L, 80L),
    scale = c("log2_cpm", "log2_cpm"),
    cell_group = c("CD4 T cells", "Monocytes"),
    slug = c("cd4_t_cells", "monocytes")
  )
}

testthat::test_that("scatter contract preserves inventory order and display names", {
  inventory <- make_scatter_inventory()

  contract <- prepare_scatter_contract(
    inventory = inventory,
    bed_paths = c("/localized/cd4_t_cells.bed.gz", "/localized/monocytes.bed.gz"),
    output_prefix = "cohort"
  )

  testthat::expect_s3_class(contract, "tbl_df")
  testthat::expect_identical(
    names(contract),
    c("cell_type", "cell_type_slug", "expression_bed", "output_prefix")
  )
  testthat::expect_equal(contract$cell_type, c("CD4 T cells", "Monocytes"))
  testthat::expect_equal(contract$cell_type_slug, c("cd4_t_cells", "monocytes"))
  testthat::expect_equal(
    contract$expression_bed,
    c("/localized/cd4_t_cells.bed.gz", "/localized/monocytes.bed.gz")
  )
  testthat::expect_equal(
    contract$output_prefix,
    c("cohort.cd4_t_cells", "cohort.monocytes")
  )
})

testthat::test_that("scatter contract rejects unsafe output-prefix tokens", {
  inventory <- make_scatter_inventory()
  unsafe_prefixes <- c(
    "../outside", "cohort name", "cohort\nname", "cohort\r", "cohort'name",
    'cohort"name', "$(touch_x)", "cohort/name", "cohort\\name",
    "cohort;name", ".", ".."
  )

  purrr::walk(unsafe_prefixes, function(unsafe_prefix) {
    testthat::expect_error(
      prepare_scatter_contract(inventory, inventory$path, unsafe_prefix),
      "safe basename token",
      info = encodeString(unsafe_prefix)
    )
  })
  testthat::expect_silent(
    prepare_scatter_contract(inventory, inventory$path, "cohort.v2-test_1")
  )
})

testthat::test_that("scatter contract requires the exact inventory schema", {
  inventory <- make_scatter_inventory()

  testthat::expect_error(
    prepare_scatter_contract(inventory[-1L], inventory$path, "cohort"),
    "exact columns"
  )
  testthat::expect_error(
    prepare_scatter_contract(
      dplyr::mutate(inventory, unexpected = "value"),
      inventory$path,
      "cohort"
    ),
    "exact columns"
  )
})

testthat::test_that("scatter contract rejects duplicate public identifiers", {
  inventory <- make_scatter_inventory()

  testthat::expect_error(
    prepare_scatter_contract(
      dplyr::mutate(inventory, cell_group = c("CD4 T cells", "CD4 T cells")),
      inventory$path,
      "cohort"
    ),
    "cell types must be unique"
  )
  testthat::expect_error(
    prepare_scatter_contract(
      dplyr::mutate(inventory, slug = c("cd4_t_cells", "cd4_t_cells")),
      inventory$path,
      "cohort"
    ),
    "slugs must be unique"
  )
  testthat::expect_error(
    prepare_scatter_contract(
      dplyr::mutate(inventory, path = c("shared.bed.gz", "shared.bed.gz")),
      inventory$path,
      "cohort"
    ),
    "paths must be unique"
  )
})

testthat::test_that("scatter contract requires canonical safe slug tokens", {
  inventory <- make_scatter_inventory()
  unsafe_slugs <- c(
    "cd4/t_cells", "../unsafe", "..", "cd4 t cells", "cd4$t_cells"
  )

  purrr::walk(unsafe_slugs, function(unsafe_slug) {
    testthat::expect_error(
      prepare_scatter_contract(
        dplyr::mutate(inventory, slug = c(unsafe_slug, "monocytes")),
        inventory$path,
        "cohort"
      ),
      "safe filename tokens"
    )
  })
})

testthat::test_that("scatter contract validates scale and dimensions", {
  inventory <- make_scatter_inventory()

  testthat::expect_error(
    prepare_scatter_contract(
      dplyr::mutate(inventory, scale = c("cpm", "log2_cpm")),
      inventory$path,
      "cohort"
    ),
    "log2_cpm"
  )
  testthat::expect_error(
    prepare_scatter_contract(
      dplyr::mutate(inventory, n_genes = c(0L, 120L)),
      inventory$path,
      "cohort"
    ),
    "n_genes must be positive"
  )
  testthat::expect_error(
    prepare_scatter_contract(
      dplyr::mutate(inventory, n_samples = c(80L, -1L)),
      inventory$path,
      "cohort"
    ),
    "n_samples must be positive"
  )
})

testthat::test_that("scatter contract requires order-aligned BED basenames", {
  inventory <- make_scatter_inventory()

  testthat::expect_error(
    prepare_scatter_contract(
      inventory,
      c("/localized/monocytes.bed.gz", "/localized/cd4_t_cells.bed.gz"),
      "cohort"
    ),
    "basenames must match"
  )
  testthat::expect_error(
    prepare_scatter_contract(inventory, inventory$path[-1L], "cohort"),
    "count must match"
  )
})

testthat::test_that("scatter contract rejects task-local public paths", {
  inventory <- make_scatter_inventory()

  testthat::expect_error(
    prepare_scatter_contract(
      dplyr::mutate(inventory, path = c("scatter/cd4_t_cells.bed.gz", "monocytes.bed.gz")),
      inventory$path,
      "cohort"
    ),
    "stable basenames"
  )
})

testthat::test_that("scatter-input CLI writes aligned metadata files", {
  root <- testthat::test_path("..", "..", "..")
  script <- file.path(
    root,
    "scripts",
    "cell_type_specific_expression",
    "prepare_scatter_inputs.R"
  )
  temporary_directory <- withr::local_tempdir()
  inventory_path <- file.path(temporary_directory, "inventory.tsv")
  bed_paths_path <- file.path(temporary_directory, "bed_paths.txt")
  output_prefix_path <- file.path(temporary_directory, "output_prefix.txt")
  output_directory <- file.path(temporary_directory, "scatter")
  readr::write_tsv(make_scatter_inventory(), inventory_path)
  writeLines(
    c("/localized/cd4_t_cells.bed.gz", "/localized/monocytes.bed.gz"),
    bed_paths_path
  )
  writeLines("cohort", output_prefix_path)

  status <- system2(
    "Rscript",
    c(
      script,
      "--inventory", inventory_path,
      "--bed-paths", bed_paths_path,
      "--output-prefix-file", output_prefix_path,
      "--output-dir", output_directory
    )
  )

  testthat::expect_identical(status, 0L)
  if (status != 0L) {
    return(invisible())
  }
  testthat::expect_equal(
    readLines(file.path(output_directory, "cell_types.txt")),
    c("CD4 T cells", "Monocytes")
  )
  testthat::expect_equal(
    readLines(file.path(output_directory, "cell_type_slugs.txt")),
    c("cd4_t_cells", "monocytes")
  )
  testthat::expect_equal(
    readLines(file.path(output_directory, "expression_beds.txt")),
    c("/localized/cd4_t_cells.bed.gz", "/localized/monocytes.bed.gz")
  )
  testthat::expect_equal(
    readLines(file.path(output_directory, "output_prefixes.txt")),
    c("cohort.cd4_t_cells", "cohort.monocytes")
  )
})

testthat::test_that("scatter-input CLI rejects a multiline output-prefix file", {
  root <- testthat::test_path("..", "..", "..")
  script <- file.path(
    root,
    "scripts",
    "cell_type_specific_expression",
    "prepare_scatter_inputs.R"
  )
  temporary_directory <- withr::local_tempdir()
  inventory_path <- file.path(temporary_directory, "inventory.tsv")
  bed_paths_path <- file.path(temporary_directory, "bed_paths.txt")
  output_prefix_path <- file.path(temporary_directory, "output_prefix.txt")
  readr::write_tsv(make_scatter_inventory(), inventory_path)
  writeLines(make_scatter_inventory()$path, bed_paths_path)
  writeLines(c("cohort", "injected"), output_prefix_path)

  status <- suppressWarnings(system2(
    "Rscript",
    c(
      script,
      "--inventory", inventory_path,
      "--bed-paths", bed_paths_path,
      "--output-prefix-file", output_prefix_path,
      "--output-dir", file.path(temporary_directory, "scatter")
    ),
    stdout = FALSE,
    stderr = FALSE
  ))

  testthat::expect_identical(status, 1L)
})

testthat::test_that("scatter-input CLI preserves and rejects carriage returns", {
  root <- testthat::test_path("..", "..", "..")
  script <- file.path(
    root,
    "scripts",
    "cell_type_specific_expression",
    "prepare_scatter_inputs.R"
  )
  temporary_directory <- withr::local_tempdir()
  inventory_path <- file.path(temporary_directory, "inventory.tsv")
  bed_paths_path <- file.path(temporary_directory, "bed_paths.txt")
  output_prefix_path <- file.path(temporary_directory, "output_prefix.txt")
  readr::write_tsv(make_scatter_inventory(), inventory_path)
  writeLines(make_scatter_inventory()$path, bed_paths_path)
  writeBin(charToRaw("cohort\r\n"), output_prefix_path)

  status <- suppressWarnings(system2(
    "Rscript",
    c(
      script,
      "--inventory", inventory_path,
      "--bed-paths", bed_paths_path,
      "--output-prefix-file", output_prefix_path,
      "--output-dir", file.path(temporary_directory, "scatter")
    ),
    stdout = FALSE,
    stderr = FALSE
  ))

  testthat::expect_identical(status, 1L)
})

testthat::test_that("scatter WDL returns aligned metadata without copying BEDs", {
  path <- testthat::test_path(
    "..",
    "..",
    "..",
    "workflows",
    "cell_type_specific_expression",
    "tasks",
    "integration.wdl"
  )
  testthat::expect_true(file.exists(path), info = path)
  if (!file.exists(path)) {
    return(invisible())
  }
  text <- paste(readLines(path, warn = FALSE), collapse = "\n")

  testthat::expect_identical(readLines(path, n = 1L, warn = FALSE), "version 1.0")
  testthat::expect_match(text, "task PrepareScatterInputs", fixed = TRUE)
  testthat::expect_match(text, "File cell_type_bed_inventory", fixed = TRUE)
  testthat::expect_match(text, "Array[File] cell_type_beds", fixed = TRUE)
  testthat::expect_match(
    text,
    "Array[String] cell_types = read_lines(\"scatter/cell_types.txt\")",
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    "Array[String] cell_type_slugs = read_lines(\"scatter/cell_type_slugs.txt\")",
    fixed = TRUE
  )
  testthat::expect_match(text, "Array[File] expression_beds = cell_type_beds", fixed = TRUE)
  testthat::expect_match(
    text,
    "Array[String] output_prefixes = read_lines(\"scatter/output_prefixes.txt\")",
    fixed = TRUE
  )
  testthat::expect_match(text, "validated_cell_count", fixed = TRUE)
  testthat::expect_match(text, "output_paths", fixed = TRUE)
  testthat::expect_match(
    text,
    "File output_prefix_file = write_lines([output_prefix])",
    fixed = TRUE
  )
  testthat::expect_match(text, "--output-prefix-file", fixed = TRUE)
  testthat::expect_false(grepl("--output-prefix '~{output_prefix}'", text, fixed = TRUE))
  testthat::expect_false(grepl("\\bcp\\b", text))
})

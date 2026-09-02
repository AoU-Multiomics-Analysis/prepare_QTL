wdl_path <- function(...) {
  testthat::test_path("..", "..", "..", "workflows", "cell_type_specific_expression", ...)
}

wdl_text <- function(...) {
  path <- wdl_path(...)
  testthat::expect_true(file.exists(path), info = path)
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

wdl_text_path <- function(path) {
  parts <- strsplit(path, "/", fixed = TRUE)[[1L]]
  do.call(wdl_text, as.list(parts))
}

testthat::test_that("standalone deconvolution sources use WDL 1.0", {
  paths <- c(
    wdl_path("deconvolution.wdl"),
    wdl_path("tasks", "dtangle.wdl"),
    wdl_path("tasks", "proportions.wdl"),
    wdl_path("tasks", "tca.wdl"),
    wdl_path("tasks", "qc.wdl")
  )

  purrr::walk(paths, function(path) {
    testthat::expect_true(file.exists(path), info = path)
    testthat::expect_identical(readLines(path, n = 1L, warn = FALSE), "version 1.0")
  })
})

testthat::test_that("workflow records and forwards the log2 pseudocount", {
  text <- wdl_text("deconvolution.wdl")

  testthat::expect_match(text, "workflow CellTypeDeconvolution", fixed = TRUE)
  testthat::expect_match(text, "Float log2_pseudocount = 0.0", fixed = TRUE)
  testthat::expect_match(text, "log2_pseudocount: log2_pseudocount", fixed = TRUE)

  calls <- c("RunDtangle", "FitTca", "ExportTcaBeds")
  purrr::walk(calls, function(call_name) {
    call_block <- stringr::str_extract(
      text,
      paste0("call [^{]+", call_name, " \\{[\\s\\S]*?\\n  \\}")
    )
    testthat::expect_false(is.na(call_block), info = call_name)
    testthat::expect_match(
      call_block,
      "log2_pseudocount = log2_pseudocount",
      fixed = TRUE,
      info = call_name
    )
  })
})

testthat::test_that("expression tasks forward the log2 pseudocount to scripts", {
  task_paths <- c("dtangle.wdl", "tca.wdl")
  purrr::walk(task_paths, function(path) {
    text <- wdl_text("tasks", path)
    testthat::expect_match(text, "Float log2_pseudocount", fixed = TRUE)
    testthat::expect_match(
      text,
      "--log2-pseudocount '~{log2_pseudocount}'",
      fixed = TRUE
    )
  })
})

testthat::test_that("migrated WDL scripts use the prepare_qtl runtime path", {
  paths <- c(
    "deconvolution.wdl", "tasks/dtangle.wdl", "tasks/proportions.wdl",
    "tasks/tca.wdl", "tasks/qc.wdl"
  )
  texts <- purrr::map(paths, wdl_text_path)
  text <- paste(texts, collapse = "\n")

  testthat::expect_match(
    text,
    "/opt/prepare_qtl/scripts/cell_type_specific_expression/",
    fixed = TRUE
  )
  testthat::expect_false(grepl("/opt/celltype", text, fixed = TRUE))
})

testthat::test_that("every container R command target is copied into the image", {
  root <- testthat::test_path("..", "..", "..")
  dockerfile <- paste(readLines(
    file.path(root, "envs", "CellTypeSpecificExpression", "Dockerfile"),
    warn = FALSE
  ), collapse = "\n")
  wdl_paths <- list.files(
    file.path(root, "workflows", "cell_type_specific_expression"),
    pattern = "[.]wdl$",
    recursive = TRUE,
    full.names = TRUE
  )
  wdl <- paste(
    unlist(purrr::map(wdl_paths, readLines, warn = FALSE), use.names = FALSE),
    collapse = "\n"
  )
  runtime_prefix <- "/opt/prepare_qtl/scripts/cell_type_specific_expression/"
  targets <- stringr::str_extract_all(
    wdl,
    paste0(stringr::fixed(runtime_prefix), "[A-Za-z0-9_./-]+[.]R")
  )[[1L]] |>
    unique()
  source_targets <- file.path(
    root,
    "scripts",
    "cell_type_specific_expression",
    sub(runtime_prefix, "", targets, fixed = TRUE)
  )

  testthat::expect_true(length(targets) > 0L)
  testthat::expect_match(
    dockerfile,
    paste(
      "scripts/cell_type_specific_expression",
      "/opt/prepare_qtl/scripts/cell_type_specific_expression"
    ),
    fixed = TRUE
  )
  testthat::expect_true(
    all(file.exists(source_targets)),
    info = paste(targets[!file.exists(source_targets)], collapse = ", ")
  )
  testthat::expect_true(any(grepl(
    "build_deconvolution_manifest[.]R$",
    targets
  )))
})

testthat::test_that("TCA exports its BED array from an ordered path file", {
  text <- wdl_text("tasks", "tca.wdl")

  testthat::expect_match(
    text,
    'Array[File] cell_type_beds = read_lines("outputs/cell_type_bed_paths.txt")',
    fixed = TRUE
  )
  testthat::expect_match(
    text,
    'File cell_type_bed_paths = "outputs/cell_type_bed_paths.txt"',
    fixed = TRUE
  )
  testthat::expect_false(grepl('glob("outputs/*.bed.gz")', text, fixed = TRUE))
})

testthat::test_that("standalone environment artifacts are present", {
  root <- testthat::test_path("..", "..", "..")
  dockerfile <- file.path(root, "envs", "CellTypeSpecificExpression", "Dockerfile")
  environment <- file.path(root, "envs", "CellTypeSpecificExpression", "environment.yml")

  testthat::expect_true(file.exists(dockerfile), info = dockerfile)
  testthat::expect_true(file.exists(environment), info = environment)
})

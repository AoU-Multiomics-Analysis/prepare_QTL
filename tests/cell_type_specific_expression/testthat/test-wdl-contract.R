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
    wdl_path("tasks", "hspe.wdl"),
    wdl_path("tasks", "expression.wdl"),
    wdl_path("tasks", "proportions.wdl"),
    wdl_path("tasks", "tca.wdl"),
    wdl_path("tasks", "qc.wdl")
  )

  purrr::walk(paths, function(path) {
    testthat::expect_true(file.exists(path), info = path)
    if (!file.exists(path)) {
      return(invisible(NULL))
    }
    testthat::expect_identical(readLines(path, n = 1L, warn = FALSE), "version 1.0")
  })
})

testthat::test_that("gene types filter expression before deconvolution", {
  standalone <- wdl_text("deconvolution.wdl")
  integrated <- wdl_text("prepare_cell_type_eQTL.wdl")
  task_path <- wdl_path("tasks", "expression.wdl")
  testthat::expect_true(file.exists(task_path), info = task_path)
  if (!file.exists(task_path)) {
    return(invisible(NULL))
  }
  task <- wdl_text("tasks", "expression.wdl")

  purrr::walk(list(standalone, integrated), function(text) {
    testthat::expect_match(
      text,
      'Array[String] gene_type = ["protein_coding", "lncRNA"]',
      fixed = TRUE
    )
  })
  testthat::expect_match(task, "task FilterExpressionGenes", fixed = TRUE)
  testthat::expect_match(task, "--gene-types '~{sep=\",\" gene_type}'", fixed = TRUE)
  testthat::expect_match(task, "File expression = \"outputs/filtered_expression.bed.gz\"", fixed = TRUE)
  testthat::expect_match(task, "File report = \"outputs/gene_type_filter_report.tsv\"", fixed = TRUE)

  testthat::expect_equal(
    stringr::str_count(
      standalone,
      stringr::fixed("expression = FilterExpressionGenes.expression")
    ),
    4L
  )
  testthat::expect_match(standalone, "gene_type = gene_type", fixed = TRUE)
  testthat::expect_match(integrated, "gene_type = gene_type", fixed = TRUE)
})

testthat::test_that("LM22 is required and precomputed proportions control estimation", {
  workflow_paths <- c("deconvolution.wdl", "prepare_cell_type_eQTL.wdl")
  purrr::walk(workflow_paths, function(path) {
    text <- wdl_text(path)
    testthat::expect_match(text, "File lm22", fixed = TRUE, info = path)
    testthat::expect_false(grepl("File? lm22", text, fixed = TRUE), info = path)
    testthat::expect_match(
      text,
      "File? precomputed_proportions",
      fixed = TRUE,
      info = path
    )
  })

  task_text <- wdl_text("tasks", "proportions.wdl")
  testthat::expect_false(grepl("File? lm22", task_text, fixed = TRUE))
  testthat::expect_false(grepl("--lm22-defined", task_text, fixed = TRUE))
  testthat::expect_match(
    task_text,
    "--precomputed-defined '~{defined(precomputed_proportions)}'",
    fixed = TRUE
  )
})

testthat::test_that("TCA task memory defaults are 256 GB", {
  text <- wdl_text("tasks", "tca.wdl")
  task_names <- c("FitTca", "ExportTcaBeds")

  purrr::walk(task_names, function(task_name) {
    input_block <- stringr::str_extract(
      text,
      paste0("task ", task_name, " \\{\\n  input \\{[\\s\\S]*?\\n  \\}")
    )
    testthat::expect_false(is.na(input_block), info = task_name)
    testthat::expect_match(
      input_block,
      'String memory = "256 GB"',
      fixed = TRUE,
      info = task_name
    )
  })
})

testthat::test_that("TCA export rebuilds expression without a persisted matrix", {
  task <- wdl_text("tasks", "tca.wdl")
  standalone <- wdl_text("deconvolution.wdl")
  integrated <- wdl_text("prepare_cell_type_eQTL.wdl")

  purrr::walk(c(task, standalone, integrated), function(text) {
    testthat::expect_false(grepl("tca_expression", text, fixed = TRUE))
  })
  testthat::expect_match(
    task,
    "--expression '~{expression}'",
    fixed = TRUE
  )
})

testthat::test_that("TCA parallel execution is explicit and disabled by default", {
  task_text <- wdl_text("tasks", "tca.wdl")
  workflow_text <- wdl_text("deconvolution.wdl")
  integrated_text <- wdl_text("prepare_cell_type_eQTL.wdl")

  testthat::expect_equal(
    stringr::str_count(task_text, stringr::fixed("Boolean parallel = false")),
    2L
  )
  testthat::expect_equal(
    stringr::str_count(task_text, stringr::fixed('if parallel then "--parallel" else ""')),
    2L
  )
  testthat::expect_equal(
    stringr::str_count(task_text, stringr::fixed("~{parallel_argument}")),
    2L
  )
  testthat::expect_match(
    workflow_text,
    "Boolean tca_parallel = false",
    fixed = TRUE
  )
  testthat::expect_equal(
    stringr::str_count(workflow_text, stringr::fixed("parallel = tca_parallel")),
    3L
  )
  testthat::expect_match(
    integrated_text,
    "Boolean tca_parallel = false",
    fixed = TRUE
  )
  testthat::expect_match(
    integrated_text,
    "tca_parallel = tca_parallel",
    fixed = TRUE
  )
})

testthat::test_that("TCA smoke fixtures retain 8 GB overrides", {
  fixture_names <- c("hspe-e2e.inputs.json", "precomputed-e2e.inputs.json")

  purrr::walk(fixture_names, function(fixture_name) {
    inputs <- jsonlite::read_json(
      testthat::test_path("..", "fixtures", fixture_name),
      simplifyVector = TRUE
    )
    testthat::expect_identical(
      inputs[["PrepareCellTypeEqtlWorkflow.fit_memory"]],
      "8 GB",
      info = fixture_name
    )
    testthat::expect_identical(
      inputs[["PrepareCellTypeEqtlWorkflow.export_memory"]],
      "8 GB",
      info = fixture_name
    )
  })
})

testthat::test_that("workflow applies the log2 pseudocount only to hspe", {
  text <- wdl_text("deconvolution.wdl")

  testthat::expect_match(text, "workflow CellTypeDeconvolution", fixed = TRUE)
  testthat::expect_match(text, "Float log2_pseudocount = 0.0", fixed = TRUE)
  testthat::expect_match(text, "log2_pseudocount = log2_pseudocount", fixed = TRUE)

  hspe_call <- stringr::str_extract(
    text,
    "call [^{]+RunHspe \\{[\\s\\S]*?\\n  \\}"
  )
  testthat::expect_match(
    hspe_call,
    "log2_pseudocount = log2_pseudocount",
    fixed = TRUE
  )
  purrr::walk(c("FitTca", "ExportTcaBeds"), function(call_name) {
    call_block <- stringr::str_extract(
      text,
      paste0("call [^{]+", call_name, " \\{[\\s\\S]*?\\n  \\}")
    )
    testthat::expect_false(
      grepl("log2_pseudocount", call_block, fixed = TRUE),
      info = call_name
    )
  })
})

testthat::test_that("effective parameters are materialized inside BuildManifest", {
  workflow_text <- wdl_text("deconvolution.wdl")
  task_text <- wdl_text("tasks", "qc.wdl")

  testthat::expect_false(grepl("write_json", workflow_text, fixed = TRUE))
  testthat::expect_match(
    workflow_text,
    "File effective_parameters_file = BuildManifest.effective_parameters_file",
    fixed = TRUE
  )
  testthat::expect_match(
    task_text,
    "--effective-parameters-output outputs/effective_parameters.json",
    fixed = TRUE
  )
  testthat::expect_match(
    task_text,
    'File effective_parameters_file = "outputs/effective_parameters.json"',
    fixed = TRUE
  )
})

testthat::test_that("filter validation and hspe receive the log2 pseudocount", {
  hspe_text <- wdl_text("tasks", "hspe.wdl")
  expression_text <- wdl_text("tasks", "expression.wdl")
  tca_text <- wdl_text("tasks", "tca.wdl")
  testthat::expect_match(hspe_text, "Float log2_pseudocount", fixed = TRUE)
  testthat::expect_match(
    hspe_text,
    "--log2-pseudocount '~{log2_pseudocount}'",
    fixed = TRUE
  )
  testthat::expect_match(expression_text, "Float log2_pseudocount", fixed = TRUE)
  testthat::expect_match(
    expression_text,
    "--log2-pseudocount '~{log2_pseudocount}'",
    fixed = TRUE
  )
  testthat::expect_false(grepl("log2_pseudocount", tca_text, fixed = TRUE))
})

testthat::test_that("migrated WDL scripts use the prepare_qtl runtime path", {
  paths <- c(
    "deconvolution.wdl", "tasks/hspe.wdl", "tasks/proportions.wdl",
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

testthat::test_that("TCA exposes generated BED files for cloud delocalization", {
  text <- wdl_text("tasks", "tca.wdl")

  testthat::expect_match(
    text,
    'Array[File] cell_type_beds = glob("outputs/*.bed.gz")',
    fixed = TRUE
  )
  testthat::expect_false(grepl("read_lines", text, fixed = TRUE))
})

testthat::test_that("BuildManifest does not relocalize cell-type BED files", {
  task_text <- wdl_text("tasks", "qc.wdl")
  workflow_text <- wdl_text("deconvolution.wdl")
  call_text <- stringr::str_match(
    workflow_text,
    "call qc_tasks[.]BuildManifest \\{[\\s\\S]*?\\n  \\}"
  )[[1L]]

  testthat::expect_false(grepl("Array[File] cell_type_beds", task_text, fixed = TRUE))
  testthat::expect_false(grepl("localized_bed_files", task_text, fixed = TRUE))
  testthat::expect_false(is.na(call_text))
  testthat::expect_false(grepl(
    "cell_type_beds = ExportTcaBeds.cell_type_beds",
    call_text,
    fixed = TRUE
  ))
})

testthat::test_that("standalone environment artifacts are present", {
  root <- testthat::test_path("..", "..", "..")
  dockerfile <- file.path(root, "envs", "CellTypeSpecificExpression", "Dockerfile")
  environment <- file.path(root, "envs", "CellTypeSpecificExpression", "environment.yml")

  testthat::expect_true(file.exists(dockerfile), info = dockerfile)
  testthat::expect_true(file.exists(environment), info = environment)
})

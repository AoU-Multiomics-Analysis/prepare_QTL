source(testthat::test_path("helper-load.R"), local = .GlobalEnv)

validation_path <- testthat::test_path(
  "..", "..", "..", "scripts", "cell_type_specific_expression", "R",
  "workflow_validation.R"
)
if (file.exists(validation_path)) {
  source(validation_path, local = .GlobalEnv)
}

testthat::test_that("proportion mode requires exactly one optional input", {
  testthat::expect_true(exists("validate_proportion_mode", mode = "function"))
  if (!exists("validate_proportion_mode", mode = "function")) {
    return(invisible(NULL))
  }

  dtangle <- validate_proportion_mode(TRUE, FALSE)
  precomputed <- validate_proportion_mode(FALSE, TRUE)
  testthat::expect_identical(
    dtangle,
    list(selected_mode = "dtangle", estimate_proportions = TRUE)
  )
  testthat::expect_identical(
    precomputed,
    list(selected_mode = "precomputed", estimate_proportions = FALSE)
  )
  testthat::expect_error(
    validate_proportion_mode(TRUE, TRUE),
    "exactly one.*lm22.*precomputed_proportions"
  )
  testthat::expect_error(
    validate_proportion_mode(FALSE, FALSE),
    "exactly one.*lm22.*precomputed_proportions"
  )
})

run_validation_cli <- function(lm22_defined, precomputed_defined) {
  script <- testthat::test_path(
    "..", "..", "..", "scripts", "cell_type_specific_expression",
    "validate_proportion_mode.R"
  )
  if (!file.exists(script)) {
    return(list(status = 127L, output = "validation script is absent"))
  }
  output_directory <- tempfile("proportion-mode-")
  dir.create(output_directory)
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      script,
      "--lm22-defined", tolower(as.character(lm22_defined)),
      "--precomputed-defined", tolower(as.character(precomputed_defined)),
      "--output-dir", output_directory
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  list(status = as.integer(status), output = paste(output, collapse = "\n"))
}

testthat::test_that("mode validation CLI rejects both optional inputs without Docker", {
  result <- run_validation_cli(TRUE, TRUE)

  testthat::expect_false(identical(result$status, 0L))
  testthat::expect_match(
    result$output,
    "exactly one.*lm22.*precomputed_proportions"
  )
})

testthat::test_that("mode validation CLI rejects neither optional input without Docker", {
  result <- run_validation_cli(FALSE, FALSE)

  testthat::expect_false(identical(result$status, 0L))
  testthat::expect_match(
    result$output,
    "exactly one.*lm22.*precomputed_proportions"
  )
})

fixture_generator_path <- function() {
  testthat::test_path(
    "..", "..", "..", "scripts", "cell_type_specific_expression",
    "generate_synthetic_fixture.R"
  )
}

run_fixture_generator <- function(output_directory) {
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    shQuote(c(fixture_generator_path(), output_directory)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  list(status = as.integer(status), output = paste(output, collapse = "\n"))
}

integrated_fixture_files <- c(
  "synthetic_expression.bed",
  "synthetic_expression_with_zero.bed",
  "synthetic.gtf",
  "synthetic_signature.tsv",
  "precomputed_proportions.tsv",
  "samples.tsv",
  "additional_covariates.tsv",
  "expected_groups.txt",
  "dtangle-e2e.inputs.json",
  "precomputed-e2e.inputs.json"
)

testthat::test_that("integrated fixtures are deterministic and checked in", {
  first_directory <- withr::local_tempdir(pattern = "fixture-first-")
  second_directory <- withr::local_tempdir(pattern = "fixture-second-")
  first_result <- run_fixture_generator(first_directory)
  second_result <- run_fixture_generator(second_directory)

  testthat::expect_identical(first_result$status, 0L, info = first_result$output)
  testthat::expect_identical(second_result$status, 0L, info = second_result$output)
  generated_paths <- file.path(first_directory, integrated_fixture_files)
  testthat::expect_true(
    all(file.exists(generated_paths)),
    info = paste(
      "Missing generated fixtures:",
      paste(integrated_fixture_files[!file.exists(generated_paths)], collapse = ", ")
    )
  )
  if (!all(file.exists(generated_paths))) {
    return(invisible(NULL))
  }

  difference <- suppressWarnings(system2(
    "diff",
    c("-rq", shQuote(first_directory), shQuote(second_directory)),
    stdout = TRUE,
    stderr = TRUE
  ))
  difference_status <- attr(difference, "status")
  if (is.null(difference_status)) {
    difference_status <- 0L
  }
  testthat::expect_identical(as.integer(difference_status), 0L)
  testthat::expect_length(difference, 0L)

  checked_in_directory <- testthat::test_path("..", "fixtures")
  checked_in_paths <- file.path(checked_in_directory, integrated_fixture_files)
  testthat::expect_true(
    all(file.exists(checked_in_paths)),
    info = paste(
      "Missing checked-in fixtures:",
      paste(integrated_fixture_files[!file.exists(checked_in_paths)], collapse = ", ")
    )
  )
  if (!all(file.exists(checked_in_paths))) {
    return(invisible(NULL))
  }

  checked_in_difference <- suppressWarnings(system2(
    "diff",
    c("-rq", shQuote(first_directory), shQuote(checked_in_directory)),
    stdout = TRUE,
    stderr = TRUE
  ))
  checked_in_status <- attr(checked_in_difference, "status")
  if (is.null(checked_in_status)) {
    checked_in_status <- 0L
  }
  testthat::expect_identical(as.integer(checked_in_status), 0L)
  testthat::expect_length(checked_in_difference, 0L)
})

testthat::test_that("end-to-end fixtures cover both proportion and pseudocount modes", {
  fixture_directory <- withr::local_tempdir(pattern = "fixture-contract-")
  generation <- run_fixture_generator(fixture_directory)
  testthat::expect_identical(generation$status, 0L, info = generation$output)

  input_paths <- file.path(
    fixture_directory,
    c("dtangle-e2e.inputs.json", "precomputed-e2e.inputs.json")
  )
  testthat::expect_true(all(file.exists(input_paths)))
  if (!all(file.exists(input_paths))) {
    return(invisible(NULL))
  }

  dtangle_inputs <- jsonlite::read_json(input_paths[[1L]], simplifyVector = TRUE)
  precomputed_inputs <- jsonlite::read_json(input_paths[[2L]], simplifyVector = TRUE)
  prefix <- "PrepareCellTypeEqtlWorkflow."
  required_shared_inputs <- c(
    "expression", "gtf", "SampleList", "AdditionalCovariates", "OutputPrefix",
    "deconvolution_docker_image", "qtl_docker_image", "log2_pseudocount"
  )
  purrr::walk(list(dtangle_inputs, precomputed_inputs), function(inputs) {
    testthat::expect_true(all(paste0(prefix, required_shared_inputs) %in% names(inputs)))
    testthat::expect_identical(
      inputs[[paste0(prefix, "deconvolution_docker_image")]],
      "cell-type-specific-expression:test"
    )
    testthat::expect_identical(
      inputs[[paste0(prefix, "qtl_docker_image")]],
      "prepare-qtl:test"
    )
    testthat::expect_false(any(grepl("residual", names(inputs), ignore.case = TRUE)))
  })

  testthat::expect_true(paste0(prefix, "lm22") %in% names(dtangle_inputs))
  testthat::expect_false(
    paste0(prefix, "precomputed_proportions") %in% names(dtangle_inputs)
  )
  testthat::expect_equal(
    dtangle_inputs[[paste0(prefix, "log2_pseudocount")]],
    0
  )
  testthat::expect_true(
    paste0(prefix, "precomputed_proportions") %in% names(precomputed_inputs)
  )
  testthat::expect_false(paste0(prefix, "lm22") %in% names(precomputed_inputs))
  testthat::expect_gt(
    precomputed_inputs[[paste0(prefix, "log2_pseudocount")]],
    0
  )

  read_expression_values <- function(path) {
    readr::read_tsv(
      file.path(fixture_directory, basename(path)),
      show_col_types = FALSE,
      progress = FALSE
    ) |>
      dplyr::select(-dplyr::all_of(c("#chr", "start", "end", "gene_id"))) |>
      as.matrix()
  }
  dtangle_values <- read_expression_values(
    dtangle_inputs[[paste0(prefix, "expression")]]
  )
  precomputed_values <- read_expression_values(
    precomputed_inputs[[paste0(prefix, "expression")]]
  )
  testthat::expect_true(all(dtangle_values > 0))
  testthat::expect_true(any(precomputed_values == 0))
  testthat::expect_false(any(precomputed_values < 0))
})

source(testthat::test_path("helper-load.R"), local = .GlobalEnv)

validation_path <- testthat::test_path(
  "..", "..", "..", "scripts", "cell_type_specific_expression", "R",
  "workflow_validation.R"
)
if (file.exists(validation_path)) {
  source(validation_path, local = .GlobalEnv)
}

assertion_helpers_path <- testthat::test_path(
  "..", "smoke", "deconvolution_assertion_helpers.R"
)
if (file.exists(assertion_helpers_path)) {
  source(assertion_helpers_path, local = .GlobalEnv)
}

testthat::test_that("precomputed proportions select whether dtangle runs", {
  testthat::expect_true(exists("validate_proportion_mode", mode = "function"))
  if (!exists("validate_proportion_mode", mode = "function")) {
    return(invisible(NULL))
  }

  dtangle <- validate_proportion_mode(FALSE)
  precomputed <- validate_proportion_mode(TRUE)
  testthat::expect_identical(
    dtangle,
    list(selected_mode = "dtangle", estimate_proportions = TRUE)
  )
  testthat::expect_identical(
    precomputed,
    list(selected_mode = "precomputed", estimate_proportions = FALSE)
  )
  testthat::expect_error(validate_proportion_mode(NA), "true or false")
  testthat::expect_error(validate_proportion_mode(1), "true or false")
})

run_validation_cli <- function(precomputed_defined) {
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

testthat::test_that("mode validation CLI selects precomputed proportions", {
  result <- run_validation_cli(TRUE)

  testthat::expect_identical(result$status, 0L, info = result$output)
  testthat::expect_match(result$output, "selected_mode=precomputed")
  testthat::expect_match(result$output, "estimate_proportions=false")
})

testthat::test_that("mode validation CLI selects dtangle when proportions are absent", {
  result <- run_validation_cli(FALSE)

  testthat::expect_identical(result$status, 0L, info = result$output)
  testthat::expect_match(result$output, "selected_mode=dtangle")
  testthat::expect_match(result$output, "estimate_proportions=true")
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

testthat::test_that("generated numeric fixture text uses eight decimal places", {
  fixture_directory <- withr::local_tempdir(pattern = "fixture-precision-")
  generation <- run_fixture_generator(fixture_directory)
  testthat::expect_identical(generation$status, 0L, info = generation$output)

  fixed_decimal_pattern <- "^-?[0-9]+[.][0-9]{8}$"
  numeric_text <- list(
    expression = readr::read_tsv(
      file.path(fixture_directory, "synthetic_expression.bed"),
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE
    )[-(1:4)],
    expression_with_zero = readr::read_tsv(
      file.path(fixture_directory, "synthetic_expression_with_zero.bed"),
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE
    )[-(1:4)],
    signature = readr::read_tsv(
      file.path(fixture_directory, "synthetic_signature.tsv"),
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE
    )[-1],
    proportions = readr::read_tsv(
      file.path(fixture_directory, "precomputed_proportions.tsv"),
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE
    )[-1],
    covariates = readr::read_tsv(
      file.path(fixture_directory, "additional_covariates.tsv"),
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE
    )[-1]
  )
  purrr::iwalk(numeric_text, function(table, label) {
    values <- unlist(table, use.names = FALSE)
    testthat::expect_true(
      all(grepl(fixed_decimal_pattern, values)),
      info = paste(label, "must use exactly eight digits after the decimal")
    )
  })
})

testthat::test_that("rounded input and grouped proportions use workflow tolerance", {
  tolerance_name <- "workflow_input_proportion_row_sum_tolerance"
  testthat::expect_true(exists(tolerance_name, inherits = TRUE))
  testthat::expect_true(exists("require_row_sums_within_tolerance", mode = "function"))
  if (!exists(tolerance_name, inherits = TRUE) ||
      !exists("require_row_sums_within_tolerance", mode = "function")) {
    return(invisible(NULL))
  }

  row_sum_tolerance <- get(tolerance_name, inherits = TRUE)
  fixture_path <- testthat::test_path("..", "fixtures", "precomputed_proportions.tsv")
  proportions <- readr::read_tsv(
    fixture_path,
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    tibble::column_to_rownames("sample_id") |>
    as.matrix()
  testthat::expect_identical(row_sum_tolerance, 1e-6)
  testthat::expect_lte(
    max(abs(rowSums(proportions) - 1)),
    row_sum_tolerance
  )
  testthat::expect_invisible(require_row_sums_within_tolerance(
    proportions,
    row_sum_tolerance,
    "LM22 proportion"
  ))

  invalid_proportions <- proportions
  invalid_proportions[1L, 1L] <- invalid_proportions[1L, 1L] +
    (1 + 2 * row_sum_tolerance - sum(invalid_proportions[1L, ]))
  testthat::expect_error(
    require_row_sums_within_tolerance(
      invalid_proportions,
      row_sum_tolerance,
      "LM22 proportion"
    ),
    "LM22 proportion rows must sum to one within 1e-06"
  )

  fixture_inputs <- jsonlite::read_json(
    testthat::test_path("..", "fixtures", "precomputed-e2e.inputs.json"),
    simplifyVector = TRUE
  )
  workflow_prefix <- "PrepareCellTypeEqtlWorkflow."
  expected <- derive_expected_proportion_outputs(
    proportions,
    fixture_inputs[[paste0(workflow_prefix, "group_mean_threshold")]],
    fixture_inputs[[paste0(workflow_prefix, "zero_floor")]]
  )
  combined_error <- max(abs(rowSums(expected$combined) - 1))
  testthat::expect_gt(combined_error, 1e-8)
  testthat::expect_lte(combined_error, row_sum_tolerance)
  testthat::expect_invisible(require_row_sums_within_tolerance(
    expected$combined,
    row_sum_tolerance,
    "combined proportion"
  ))

  invalid_combined <- expected$combined
  invalid_combined[1L, 1L] <- invalid_combined[1L, 1L] +
    (1 + 2 * row_sum_tolerance - sum(invalid_combined[1L, ]))
  testthat::expect_error(
    require_row_sums_within_tolerance(
      invalid_combined,
      row_sum_tolerance,
      "combined proportion"
    ),
    "combined proportion rows must sum to one within 1e-06"
  )
  testthat::expect_lt(
    max(abs(rowSums(expected$tca_weights) - 1)),
    1e-8
  )

  assertion_text <- readLines(
    testthat::test_path("..", "smoke", "assert_deconvolution_outputs.R"),
    warn = FALSE
  ) |>
    paste(collapse = "\n")
  row_sum_helper_calls <- gregexpr(
    "require_row_sums_within_tolerance\\(",
    assertion_text
  )[[1L]]
  row_sum_helper_call_count <- if (identical(row_sum_helper_calls, -1L)) {
    0L
  } else {
    length(row_sum_helper_calls)
  }
  testthat::expect_identical(row_sum_helper_call_count, 2L)
})

testthat::test_that("smoke expectations derive canonical combined groups and weights", {
  required_helpers <- c(
    "canonical_lm22_group_map", "derive_expected_proportion_outputs",
    "require_matrix_equal"
  )
  purrr::walk(
    required_helpers,
    ~ testthat::expect_true(exists(.x, mode = "function"), info = .x)
  )
  if (!all(vapply(required_helpers, exists, logical(1), mode = "function"))) {
    return(invisible(NULL))
  }

  cell_types <- unlist(canonical_lm22_group_map(), use.names = FALSE)
  proportions <- matrix(
    1 / 22,
    nrow = 2L,
    ncol = 22L,
    dimnames = list(c("S1", "S2"), cell_types)
  )
  proportions["S2", ] <- 0
  proportions["S2", "B cells naive"] <- 0.4
  proportions["S2", "T cells CD4 naive"] <- 0.6

  expected <- derive_expected_proportion_outputs(
    proportions,
    mean_threshold = 0.08,
    zero_floor = 0.1
  )
  group_counts <- c(3, 5, 1, 1, 2, 4, 1, 1, 2, 2)
  expected_combined <- rbind(
    S1 = group_counts / 22,
    S2 = c(0.4, 0.6, rep(0, 8L))
  )
  colnames(expected_combined) <- names(canonical_lm22_group_map())
  expected_weights <- rbind(
    S1 = c(3 / 12, 5 / 12, 4 / 12),
    S2 = c(0.4 / 1.1, 0.6 / 1.1, 0.1 / 1.1)
  )
  colnames(expected_weights) <- c("B cells", "CD4 T cells", "Monocyte/myeloid")

  testthat::expect_identical(colnames(expected$combined), colnames(expected_combined))
  testthat::expect_equal(expected$combined, expected_combined, tolerance = 1e-12)
  testthat::expect_identical(
    expected$retained_groups,
    c("B cells", "CD4 T cells", "Monocyte/myeloid")
  )
  testthat::expect_equal(expected$tca_weights, expected_weights, tolerance = 1e-12)
  testthat::expect_invisible(
    require_matrix_equal(expected$combined, expected_combined, 1e-12, "combined")
  )
  changed <- expected_combined
  changed[[1L]] <- changed[[1L]] + 1e-4
  testthat::expect_error(
    require_matrix_equal(expected$combined, changed, 1e-8, "combined"),
    "tolerance"
  )
  reordered <- expected_combined[, rev(seq_len(ncol(expected_combined))), drop = FALSE]
  testthat::expect_error(
    require_matrix_equal(expected$combined, reordered, 1e-8, "combined"),
    "row or column order"
  )
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
    "expression", "gtf", "lm22", "SampleList", "AdditionalCovariates", "OutputPrefix",
    "deconvolution_docker_image", "qtl_docker_image", "log2_pseudocount",
    "gene_type"
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
    testthat::expect_identical(
      inputs[[paste0(prefix, "gene_type")]],
      c("protein_coding", "lncRNA")
    )
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
  testthat::expect_true(paste0(prefix, "lm22") %in% names(precomputed_inputs))
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

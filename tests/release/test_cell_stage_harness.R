#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
definitions <- new.env(parent = globalenv())
purrr::walk(parse("tests/release/test_cell_stage.R"), function(expression) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.call(expression[[3L]]) && identical(expression[[3L]][[1L]], as.name("function"))) {
    eval(expression, definitions)
  }
})

testthat::test_that("stage workspace reads bundled code, not mounted candidate code", {
  testthat::expect_true(exists("prepare_stage_workspace", definitions, inherits = FALSE))
  if (!exists("prepare_stage_workspace", definitions, inherits = FALSE)) return(invisible(NULL))
  directory <- tempfile("stage-harness-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  source <- file.path(directory, "source")
  bundled <- file.path(directory, "bundled")
  work <- file.path(directory, "work")
  purrr::walk(c("tests", "workflows", "scripts/cell_type_specific_expression"),
             ~ dir.create(file.path(source, .x), recursive = TRUE))
  dir.create(bundled)
  writeLines("bundled", file.path(bundled, "bootstrap.R"))
  writeLines("candidate", file.path(source, "scripts/cell_type_specific_expression/bootstrap.R"))
  writeLines("fixture", file.path(source, "tests/fixture.txt"))
  writeLines("version 1.0", file.path(source, "workflows/main.wdl"))
  definitions$prepare_stage_workspace(source, bundled, work)
  testthat::expect_identical(readLines(file.path(work, "scripts/cell_type_specific_expression/bootstrap.R")), "bundled")
  testthat::expect_identical(readLines(file.path(work, "tests/fixture.txt")), "fixture")
  testthat::expect_identical(readLines(file.path(source, "scripts/cell_type_specific_expression/bootstrap.R")), "candidate")
  testthat::expect_error(definitions$prepare_stage_workspace(source, file.path(directory, "missing"),
                                                          file.path(directory, "bad")), "Missing bundled")
})

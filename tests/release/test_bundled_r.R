#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(RNOmni))
main <- function() {
family <- commandArgs(trailingOnly = TRUE)[[1L]]
stopifnot(family %in% c("proteomics", "splicing"))
script_root <- Sys.getenv("BUNDLED_SCRIPT_ROOT", "/tmp")
script <- file.path(script_root, if (family == "proteomics") "PrepareProteomics.R" else "PrepareSpliceData.R")

# Load the actual bundled function definitions without running data downloads.
# The complete proteomics workflow still needs the external BioMart service.
definitions <- parse(script)
functions <- new.env(parent = globalenv())
purrr::walk(definitions, function(expression) {
  if (is.call(expression) && identical(expression[[1L]], as.name("<-")) &&
      is.call(expression[[3L]]) && identical(expression[[3L]][[1L]], as.name("function"))) {
    eval(expression, envir = functions)
  }
})
values <- c(2, 4, 7, 9, 13)
stopifnot(isTRUE(all.equal(functions$transform_phenotype(values, FALSE), as.numeric(scale(values)))))
stopifnot(isTRUE(all.equal(functions$transform_phenotype(values, TRUE), RNOmni::RankNorm(values))))

if (family == "proteomics") {
  data <- tibble(SampleID = c("bad1", "bad2"), ResearchID = c("s1", "s2"))
  stopifnot(identical(functions$select_sample_id_column(data, c("s1", "s2")), "ResearchID"))
  stopifnot(inherits(try(functions$select_sample_id_column(data, "missing"), silent = TRUE), "try-error"))
  stopifnot(system2("Rscript", c(shQuote(script), "--help")) == 0L)
  stopifnot(system2("Rscript", c(shQuote(file.path(script_root, "NormalizeProteomics.R")), "--help")) == 0L)
} else {
  directory <- tempfile("splicing-smoke-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  samples <- paste0("s", 1:5)
  bed <- tibble(`#chr` = "chr1", start = c(100, 200, 300), end = c(101, 201, 301),
                ID = paste0("junction_ENSG", 1:3)) |>
    bind_cols(as_tibble(matrix(c(1, 3, 9, 2, 5, 8, 4, 8, 6, 7, 6, 4, 9, 7, 2), nrow = 3,
                              dimnames = list(NULL, samples))))
  input <- file.path(directory, "input.bed")
  sample_file <- file.path(directory, "samples.tsv")
  prefix <- file.path(directory, "test")
  write_tsv(bed, input)
  write_tsv(tibble(ID = samples), sample_file)
  stopifnot(system2("Rscript", c(shQuote(script), "--SpliceData", shQuote(input), "--SampleList",
                               shQuote(sample_file), "--OutputPrefix", shQuote(prefix))) == 0L)
  for (suffix in c("INT", "scaled", "raw")) {
    output <- read_tsv(paste0(prefix, ".splicing.", suffix, ".bed.gz"), show_col_types = FALSE)
    stopifnot(nrow(output) == 3L, identical(names(output)[-(1:4)], samples))
    if (suffix == "scaled") stopifnot(max(abs(rowMeans(as.matrix(output[, samples])))) < 1e-10)
  }
}
message("Bundled ", family, " smoke checks passed")
}
main()

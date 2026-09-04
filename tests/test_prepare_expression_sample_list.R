#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1]]
test_path <- gsub("~+~", " ", sub("^--file=", "", file_arg), fixed = TRUE)
root <- normalizePath(file.path(dirname(test_path), ".."))
script <- Sys.getenv("PREPARE_EXPRESSION_SCRIPT",
  unset = file.path(root, "scripts/expression/PrepareExpression.R"))

check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

run_tests <- function() {
  directory <- tempfile("sample-list-headers-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE), add = TRUE)
  samples_file <- file.path(directory, "samples.tsv")
  expected_samples <- c("0007", "S3", "S2")
  writeLines(c("research_id", expected_samples), samples_file)
  bed <- tibble::tibble(`#chr` = c("chr1", "chr2"), start = c(0L, 10L),
    end = c(1L, 11L), gene_id = c("g1", "g2"),
    S2 = c(1, 3), `0007` = c(4, 2), S3 = c(7, 6))
  bed_path <- file.path(directory, "expression.bed.gz")
  readr::write_tsv(bed, bed_path)
  gct <- file.path(directory, "counts.gct")
  writeLines(c("#1.2", "2\t3", "Name\tDescription\tS2\t0007\tS3",
    "g1\tone\t10\t20\t40", "g2\ttwo\t40\t30\t20"), gct)
  gtf <- file.path(directory, "genes.gtf")
  writeLines(c(
    'chr1\ttest\tgene\t1\t10\t.\t+\t.\tgene_id "g1";',
    'chr2\ttest\tgene\t11\t20\t.\t+\t.\tgene_id "g2";'), gtf)

  # Real CLI calls catch header handling that is not wired into an input mode.
  for (mode in c("CpmBed", "Log2CpmBed", "CountGCT")) {
    prefix <- file.path(directory, mode)
    args <- c(script, paste0("--", mode), if (mode == "CountGCT") gct else bed_path,
      "--SampleList", samples_file, "--OutputPrefix", prefix)
    if (mode == "CountGCT") args <- c(args, "--AnnotationGTF", gtf)
    logs <- suppressWarnings(system2("Rscript", shQuote(args), stdout = TRUE, stderr = TRUE))
    check(is.null(attr(logs, "status")), paste(c(mode, tail(logs, 10)), collapse = "\n"))
    check(any(grepl("Number of sample in sample list:3", logs, fixed = TRUE)),
      "The header was included in the sample count")
    for (suffix in c("raw", "INT", "scaled")) {
      output <- readr::read_tsv(paste0(prefix, ".expression.", suffix, ".bed.gz"),
        show_col_types = FALSE)
      check(identical(names(output)[-(1:4)], expected_samples),
        paste(mode, suffix, "changed sample order or leading zeros"))
      if (suffix == "raw" && mode != "CountGCT") {
        check(isTRUE(all.equal(as.matrix(output[expected_samples]),
          as.matrix(bed[expected_samples]), tolerance = 0)), "Supplied expression values changed")
      }
    }
  }

  # Load the real reader without executing the CLI or its unrelated packages.
  definition <- Filter(function(expr) {
    is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      identical(expr[[2]], as.name("read_expression_sample_list"))
  }, as.list(parse(script)))
  check(length(definition) == 1L, "The sample-list reader is missing")
  eval(definition[[1L]])
  read_samples <- function(lines, available = expected_samples) {
    writeLines(lines, samples_file)
    read_expression_sample_list(samples_file, available, "CPM BED")
  }
  for (header in c("research_id", "RESEARCH_ID", "sample_id", "Sample_ID", "ID", "id")) {
    check(identical(read_samples(c(header, expected_samples)), expected_samples),
      paste("Header was not recognized:", header))
  }
  check(identical(read_samples(expected_samples), expected_samples), "Headerless IDs changed")
  check(identical(read_samples(c('"research_id"', '"0007"', '"S3"', '"S2"')),
    expected_samples), "Quoted IDs changed")
  for (real_id in c("ID", "research_id", "sample_id")) {
    real_samples <- c(real_id, expected_samples)
    check(identical(read_samples(real_samples, real_samples), real_samples),
      "A real first sample was mistaken for a header")
  }
  invalid <- list(
    list(c("research_id", "0007", "0007"), "duplicate"),
    list(c("research_id", "0007", "", "S2"), "blank"),
    list(c("research_id", "0007", "   ", "S2"), "blank"),
    list(character(), "empty"),
    list("research_id", "no sample IDs"),
    list(c("missing_sample", expected_samples), "missing_sample"),
    list(c("0007", "research_id", "S2"), "research_id"),
    list(c("sample\tgroup", "0007\tA", "S2\tB"), "one column")
  )
  for (case in invalid) {
    error <- tryCatch({ read_samples(case[[1]]); "NO ERROR" }, error = conditionMessage)
    check(grepl(case[[2]], error, fixed = TRUE), paste("Wrong validation result:", error))
  }
  message("Sample-list header tests passed for CPM, log2 CPM, and raw counts")
}

run_tests()

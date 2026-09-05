invisible(lapply(file.path(script_root, "R", c(
  "bed_outputs.R", "tca_stage.R", "expression_bed.R", "reference_expression.R",
  "reference_filter.R", "reference_filter_plots.R",
  "reference_filter_io.R"
)), source, local = .GlobalEnv))

testthat::test_that("reference sample mapping pools samples and rejects incomplete lineages", {
  mapped <- map_reference_samples(c("NveB.1", "NveB.2", "MemB.1", "CD4T.1"))
  testthat::expect_equal(mapped$population, c("NveB", "NveB", "MemB", "CD4T"))
  testthat::expect_equal(mapped$cell_type, c("B cells", "B cells", "B cells", "CD4 T cells"))
  testthat::expect_error(map_reference_samples(c("NveB.1", "CD4T.1")), "missing.*MemB")
  testthat::expect_error(map_reference_samples("mystery.1"), "Unknown")
})

testthat::test_that("gene matching removes numeric versions and rejects collisions", {
  testthat::expect_equal(reference_gene_key(c("ENSG1.2", "ABC.2")), c("ENSG1", "ABC.2"))
  testthat::expect_error(reference_gene_key(c("ENSG1.1", "ENSG1.2"), unique = TRUE), "duplicate")
})

testthat::test_that("reference preparation matches full-matrix edgeR normalization", {
  testthat::skip_if_not_installed("edgeR")
  counts <- matrix(c(100, 10, 50, 5, 4, 80, 3, 70, 40, 20, 30, 10), nrow = 3, byrow = TRUE,
                   dimnames = list(c("ENSG1.1", "ENSG2", "ENSG3"),
                                   c("NveB.1", "NveB.2", "MemB.1", "CD4T.1")))
  expected_dge <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts))
  expected <- edgeR::cpm(expected_dge, log = FALSE)
  result <- prepare_reference_matrix(counts)
  testthat::expect_equal(result$cpm, expected, tolerance = 1e-12)
  summary <- summarize_reference_expression(result$cpm, result$sample_map)
  b_gene <- dplyr::filter(summary, gene_id == "ENSG1.1", cell_type == "B cells")
  testthat::expect_equal(b_gene$mean_log2_cpm1, mean(log2(expected[1, 1:3] + 1)))
  testthat::expect_equal(b_gene$n_samples, 3)
})

testthat::test_that("negative and log statistics use sample-level values", {
  x <- matrix(c(0, 3, 15, -1e-12, 2, 4), nrow = 2, byrow = TRUE,
              dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
  result <- summarize_bed_values(x)
  testthat::expect_equal(result$mean_log2_cpm1[[1]], 2)
  testthat::expect_equal(result$negative_count, c(0L, 1L))
  testthat::expect_true(is.na(result$mean_negative_cpm[[1]]))
  testthat::expect_equal(result$mean_negative_cpm[[2]], -1e-12)
  testthat::expect_false(result$nonnegative[[2]])
  testthat::expect_error(summarize_bed_values(matrix(c(1, Inf), nrow = 1)), "finite")
})

testthat::test_that("reference decisions use strict thresholds and explicit unavailable states", {
  bed <- tibble::tibble(gene_id = c("a", "b", "c"), nonnegative = TRUE,
                        mean_log2_cpm1 = c(0.01, 0.02, 1))
  ref <- tibble::tibble(gene_id = c("a", "b"), cell_type = "B cells",
                        mean_log2_cpm1 = c(1, 0.01))
  result <- apply_reference_rules(bed, "B cells", ref, 0.01)
  testthat::expect_equal(result$retained, c(FALSE, FALSE, FALSE))
  testthat::expect_equal(result$comparison_status, c("compared", "compared", "reference_gene_unmatched"))
  unsupported <- apply_reference_rules(bed, "Mast cells", ref, 0.01)
  testthat::expect_true(all(unsupported$retained))
  testthat::expect_true(all(unsupported$comparison_status == "no_reference_cell_type"))
  absent <- apply_reference_rules(bed, "B cells", NULL, 0.01)
  testthat::expect_true(all(absent$retained))
  testthat::expect_true(all(absent$comparison_status == "reference_not_provided"))
  multi_lineage <- dplyr::bind_rows(ref,
    dplyr::mutate(ref, cell_type = "CD4 T cells", mean_log2_cpm1 = 10))
  testthat::expect_equal(
    apply_reference_rules(bed, "B cells", multi_lineage, 0.01)$reference_mean_log2_cpm1[1:2],
    c(1, 0.01)
  )
})

testthat::test_that("baseline regression uses rstandard and one pass cutoff", {
  comparison <- tibble::tibble(
    gene_id = paste0("g", 1:6),
    reference_mean_log2_cpm1 = 1:6,
    deconvolution_mean_log2_cpm1 = c(1, 2, 3, 4, 5, 12)
  )
  fit <- fit_reference_regression(comparison)
  expected <- stats::lm(deconvolution_mean_log2_cpm1 ~ reference_mean_log2_cpm1,
                        data = comparison)
  testthat::expect_equal(fit$genes$standardized_residual, stats::rstandard(expected))
  filtered <- apply_residual_cutoff(fit$genes, 1.5)
  testthat::expect_equal(filtered$residual_excluded,
                         abs(stats::rstandard(expected)) > 1.5)
  testthat::expect_error(fit_reference_regression(comparison[1:2, ]), "three")
  perfect <- fit_reference_regression(dplyr::mutate(comparison, deconvolution_mean_log2_cpm1 = 2 * reference_mean_log2_cpm1))
  testthat::expect_true(all(is.na(perfect$genes$standardized_residual)))
  testthat::expect_error(apply_residual_cutoff(perfect$genes, 2), "undefined")
})

testthat::test_that("filter CLI preserves retained BED rows across chunks and paths with spaces", {
  tmp <- tempfile("reference filter ")
  dir.create(tmp)
  bed <- file.path(tmp, "B cells input.bed.gz")
  input <- tibble::tibble(`#chr` = c("1", "1", "2"), start = c(0L, 10L, 20L),
                          end = c(1L, 11L, 21L), gene_id = c("g1", "g2", "g3"),
                          s1 = c(0, -1e-12, 4), s2 = c(3, 2, 4), s3 = c(15, 4, 4))
  readr::write_tsv(input, bed)
  inventory <- tibble::tibble(
    logical_name = "cell_type_bed", path = basename(bed),
    sha256 = digest::digest(file = bed, algo = "sha256", serialize = FALSE),
    n_genes = 3L, n_samples = 3L, scale = "cpm", cell_group = "B cells", slug = "b_cells")
  inventory_path <- file.path(tmp, "inventory.tsv")
  readr::write_tsv(inventory, inventory_path)
  bed_list <- file.path(tmp, "bed paths with spaces.txt")
  writeLines(bed, bed_list)
  output <- file.path(tmp, "output with spaces")
  status <- system2("Rscript", c(shQuote(file.path(script_root, "filter_cell_type_beds.R")),
                                  "--inventory", shQuote(inventory_path),
                                  "--bed-list", shQuote(bed_list), "--chunk-size", "2",
                                  "--output-dir", shQuote(output)))
  testthat::expect_equal(status, 0L)
  filtered <- readr::read_tsv(file.path(output, "beds", "b_cells.filtered.bed.gz"), show_col_types = FALSE)
  testthat::expect_equal(filtered$gene_id, c("g1", "g3"))
  testthat::expect_equal(as.character(filtered[["#chr"]]), c("1", "2"))
  testthat::expect_equal(dplyr::select(filtered, start, end, s1, s2, s3),
                         dplyr::select(input[c(1, 3), ], start, end, s1, s2, s3))
  testthat::expect_true(file.exists(file.path(output, "plots", "negative_overview.pdf")))
  new_inventory <- readr::read_tsv(file.path(output, "filtered_inventory.tsv"), show_col_types = FALSE)
  testthat::expect_identical(names(new_inventory), scatter_inventory_columns)
  testthat::expect_equal(new_inventory$n_genes, 2)
  testthat::expect_equal(new_inventory$sha256,
    digest::digest(file = file.path(output, "beds", "b_cells.filtered.bed.gz"), algo = "sha256", serialize = FALSE))
  testthat::expect_identical(readLines(file.path(output, "filtered_beds.txt")),
                             file.path(output, "beds", "b_cells.filtered.bed.gz"))
})

testthat::test_that("reference preparation CLI writes normalized provenance outputs", {
  testthat::skip_if_not_installed("edgeR")
  tmp <- tempfile("reference preparation ")
  dir.create(tmp)
  counts_path <- file.path(tmp, "raw counts.tsv.gz")
  counts <- tibble::tibble(gene_id = c("ENSG1.1", "ENSG2"),
    NveB.1 = c(10L, 4L), NveB.2 = c(20L, 2L), MemB.1 = c(3L, 30L), CD4T.1 = c(5L, 8L))
  readr::write_tsv(counts, counts_path)
  output <- file.path(tmp, "reference output")
  status <- system2("Rscript", c(shQuote(file.path(script_root, "prepare_haemopedia.R")),
                                  shQuote(counts_path), shQuote(output)))
  testthat::expect_equal(status, 0L)
  summary <- readr::read_tsv(file.path(output, "reference_summary.tsv.gz"), show_col_types = FALSE)
  testthat::expect_identical(names(summary), c("gene_id", "cell_type", "n_samples",
    "mean_log2_cpm1", "median_log2_cpm1"))
  samples <- readr::read_tsv(file.path(output, "reference_samples.tsv"), show_col_types = FALSE)
  testthat::expect_true(all(c("library_size", "norm_factor") %in% names(samples)))
  metadata <- jsonlite::read_json(file.path(output, "reference_metadata.json"), simplifyVector = TRUE)
  testthat::expect_match(metadata$normalization, "full raw count matrix")
  testthat::expect_equal(metadata$input$n_genes, 2L)
})

testthat::test_that("filtered BED writing preserves the original text exactly", {
  directory <- tempfile("exact bed ")
  dir.create(directory)
  input <- file.path(directory, "input.bed")
  output <- file.path(directory, "output.bed.gz")
  lines <- c(
    "#chr\tstart\tend\tgene_id\ts1\ts2",
    "1\t0\t1\tg1\t0.12345678901234567\t10000000000000003",
    "1\t1\t2\tg2\t2.0000000000000000\t3e-09",
    "2\t2\t3\tg3\t4.50\t6.700000000000001"
  )
  writeLines(lines, input)
  write_filtered_bed(input, output, c("g1", "g3"), chunk_size = 2L)
  testthat::expect_identical(readLines(gzfile(output)), lines[c(1L, 2L, 4L)])
})

testthat::test_that("filter metrics preserve unavailable reference states and sample counts", {
  testthat::expect_error(validate_residual_cutoff(0), "positive")
  testthat::expect_error(validate_residual_cutoff(Inf), "positive")
  testthat::expect_null(validate_residual_cutoff(NULL))
  metric <- make_filter_metric(
    tibble::tibble(n_genes = 0L, pearson_r = NA_real_, spearman_rho = NA_real_,
      r_squared = NA_real_, intercept = NA_real_, slope = NA_real_),
    cell_type = "B cells", slug = "b_cells",
    comparison_status = "reference_cell_type_unavailable",
    n_original = 4L, n_negative_excluded = 1L, n_reference_excluded = 0L,
    n_residual_excluded = 0L, n_retained = 3L,
    deconvolution_n_samples = 10L, reference_n_samples = NA_integer_,
    metric_set = "baseline"
  )
  testthat::expect_equal(metric$comparison_status, "reference_cell_type_unavailable")
  testthat::expect_equal(metric$deconvolution_n_samples, 10L)
  testthat::expect_true(is.na(metric$reference_n_samples))
  other_lineage <- tibble::tibble(gene_id = "g1", cell_type = "CD4 T cells",
    n_samples = 4L, mean_log2_cpm1 = 1, median_log2_cpm1 = 1)
  testthat::expect_identical(reference_comparison_status(other_lineage, "B cells"),
                             "reference_cell_type_unavailable")
})

testthat::test_that("a supplied reference must contain each supported inventory lineage", {
  directory <- tempfile("lineage validation ")
  dir.create(directory)
  bed <- file.path(directory, "cell.bed.gz")
  values <- tibble::tibble(`#chr` = "1", start = 0:2, end = 1:3,
    gene_id = paste0("g", 1:3), s1 = c(1, 2, 3), s2 = c(2, 3, 5))
  readr::write_tsv(values, bed)
  inventory <- tibble::tibble(
    logical_name = "cell_type_bed", path = basename(bed),
    sha256 = digest::digest(file = bed, algo = "sha256", serialize = FALSE),
    n_genes = 3L, n_samples = 2L, scale = "cpm",
    cell_group = "B cells", slug = "b_cells"
  )
  cd4_reference <- tibble::tibble(
    gene_id = paste0("g", 1:3), cell_type = "CD4 T cells", n_samples = 2L,
    mean_log2_cpm1 = c(1, 2, 3), median_log2_cpm1 = c(1, 2, 3)
  )
  testthat::expect_error(
    filter_cell_type_beds(inventory, bed, file.path(directory, "missing-no-cutoff"),
                          reference_summary = cd4_reference),
    "missing.*B cells"
  )
  testthat::expect_error(
    filter_cell_type_beds(inventory, bed, file.path(directory, "missing-with-cutoff"),
                          reference_summary = cd4_reference, residual_cutoff = 2),
    "missing.*B cells"
  )

  cd4_inventory <- dplyr::mutate(inventory, cell_group = "CD4 T cells", slug = "cd4_t_cells")
  testthat::expect_no_error(
    filter_cell_type_beds(cd4_inventory, bed, file.path(directory, "valid-subset"),
                          reference_summary = cd4_reference)
  )
})

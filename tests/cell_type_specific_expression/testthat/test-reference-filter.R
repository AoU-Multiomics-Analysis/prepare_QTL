invisible(lapply(file.path(script_root, "R", c(
  "reference_expression.R", "reference_filter.R", "reference_filter_plots.R",
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
  config <- list(inventory = inventory_path, bed_paths = list(bed), reference_summary = NULL,
                 min_mean_log2_cpm1 = 0.01, residual_cutoff = NULL, chunk_size = 2L)
  config_path <- file.path(tmp, "config with spaces.json")
  jsonlite::write_json(config, config_path, auto_unbox = TRUE, null = "null")
  output <- file.path(tmp, "output with spaces")
  status <- system2("Rscript", c(shQuote(file.path(script_root, "filter_cell_type_beds.R")),
                                  shQuote(config_path), shQuote(output)))
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

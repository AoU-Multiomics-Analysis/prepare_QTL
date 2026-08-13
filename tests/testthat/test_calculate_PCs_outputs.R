source(file.path("..", "..", "scripts", "common", "PCOutputUtils.R"))

testthat::test_that("full output retains every rotated PC and selected is a prefix", {
    rotated <- data.frame(
        PC1 = c(1, 2), PC2 = c(3, 4), PC3 = c(5, 6),
        row.names = c("Xsample_a", "Xsample_b")
    )
    outputs <- format_pca_outputs(rotated, n_pcs = 2)

    testthat::expect_identical(names(outputs$all), c("ID", "PC1", "PC2", "PC3"))
    testthat::expect_identical(outputs$all$ID, c("sample_a", "sample_b"))
    testthat::expect_identical(names(outputs$selected), c("ID", "PC1", "PC2"))
    testthat::expect_equal(outputs$selected[, -1], outputs$all[, c("PC1", "PC2")])
})

testthat::test_that("writer creates selected and all files", {
    rotated <- data.frame(PC1 = c(1, 2), PC2 = c(3, 4), row.names = c("Xone", "Xtwo"))
    selected_path <- tempfile(fileext = ".tsv")
    all_path <- tempfile(fileext = ".all.tsv")
    write_pca_outputs(rotated, n_pcs = 1, selected_path, all_path)

    testthat::expect_true(file.exists(selected_path))
    testthat::expect_true(file.exists(all_path))
    testthat::expect_equal(ncol(readr::read_tsv(selected_path, show_col_types = FALSE)), 2)
    testthat::expect_equal(ncol(readr::read_tsv(all_path, show_col_types = FALSE)), 3)
})

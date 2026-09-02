repo_path <- function(...) {
    testthat::test_path("..", "..", ...)
}

read_document <- function(path) {
    if (!file.exists(path)) {
        return("")
    }
    paste(readLines(path, warn = FALSE), collapse = "\n")
}

expect_documented <- function(text, patterns) {
    normalized_text <- gsub("\\s+", " ", text)
    for (pattern in patterns) {
        testthat::expect_match(normalized_text, pattern, perl = TRUE)
    }
}

testthat::test_that("Dockstore publishes both cell-type workflow entry points", {
    dockstore <- read_document(repo_path(".dockstore.yml"))
    entries <- c(
        cell_type_deconvolution =
            "/workflows/cell_type_specific_expression/deconvolution.wdl",
        prepare_cell_type_eQTL =
            "/workflows/cell_type_specific_expression/prepare_cell_type_eQTL.wdl"
    )

    for (workflow_name in names(entries)) {
        entry_pattern <- paste0(
            "primaryDescriptorPath:\\s*", entries[[workflow_name]],
            "\\s*\\n\\s*name:\\s*", workflow_name, "(?:\\s|$)"
        )
        testthat::expect_match(dockstore, entry_pattern, perl = TRUE)
    }
})

testthat::test_that("repository catalogs link to the cell-type workflow guide", {
    guide_path <- repo_path("docs", "cell-type-specific-expression.md")
    testthat::expect_true(file.exists(guide_path))

    for (path in c(repo_path("README.md"), repo_path("workflows", "README.md"))) {
        contents <- read_document(path)
        testthat::expect_match(
            contents,
            "cell-type-specific-expression\\.md",
            perl = TRUE,
            info = paste("Missing cell-type guide link in", path)
        )
    }
})

testthat::test_that("guide states the expression and proportion input contracts", {
    guide <- read_document(
        repo_path("docs", "cell-type-specific-expression.md")
    )

    expect_documented(guide, c(
        "linear CPM",
        "finite.*nonnegative",
        "log2_pseudocount.*defaults to.*0",
        "log2_pseudocount.*0.*strictly positive",
        "[Zz]ero.*only.*log2_pseudocount.*greater than zero",
        "[Nn]egative.*always.*invalid",
        "log2\\(CPM \\+ log2_pseudocount\\).*exactly once",
        "LM22.*dtangle",
        "precomputed LM22 proportions",
        "exactly one.*lm22.*precomputed_proportions"
    ))
})

testthat::test_that("guide states the scattered QTL preparation contract", {
    guide <- read_document(
        repo_path("docs", "cell-type-specific-expression.md")
    )

    expect_documented(guide, c(
        "Log2CpmBed",
        "does not apply.*second log2",
        "AdditionalCovariates.*required",
        "independent.*connectivity filtering",
        "cell type.*INT and scaled branch",
        "stable basenames",
        "Array\\[File\\].*authoritative",
        "raw or residualized BED"
    ))
})

testthat::test_that("guide publishes the exact QTL manifest schema", {
    guide_lines <- readLines(
        repo_path("docs", "cell-type-specific-expression.md"),
        warn = FALSE
    )
    manifest_header <- "| Column | Meaning |"
    header_index <- match(manifest_header, guide_lines)
    testthat::expect_false(is.na(header_index))

    table_lines <- guide_lines[(header_index + 2L):length(guide_lines)]
    table_lines <- table_lines[grepl("^\\| `[^`]+` \\|", table_lines)]
    documented_columns <- sub("^\\| `([^`]+)` \\|.*$", "\\1", table_lines)
    expected_columns <- c(
        "cell_type",
        "cell_type_slug",
        "int_bed",
        "scaled_bed",
        "int_phenotype_pcs",
        "int_phenotype_pcs_all",
        "scaled_phenotype_pcs",
        "scaled_phenotype_pcs_all",
        "int_merged_covariates",
        "scaled_merged_covariates",
        "int_connectivity_outliers",
        "scaled_connectivity_outliers"
    )

    testthat::expect_identical(documented_columns, expected_columns)
})

testthat::test_that("guide states canonical ownership and the two-branch merge order", {
    guide <- read_document(
        repo_path("docs", "cell-type-specific-expression.md")
    )

    expect_documented(guide, c(
        "prepare_QTL.*only maintained.*canonical",
        "feat/log2-cpm-bed-input.*first",
        "codex/cell-type-specific-expression.*after",
        "old CellTypeDeconvolution repository.*deprecated.*after.*prepare_QTL.*main",
        "[Dd]o not delete"
    ))
})

repo_path <- function(...) {
    file.path("..", "..", ...)
}

read_file <- function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
}

expect_contains_all <- function(path, patterns) {
    contents <- read_file(path)
    for (pattern in patterns) {
        testthat::expect_match(
            contents,
            pattern,
            perl = TRUE,
            info = paste("Missing pattern in", path, ":", pattern)
        )
    }
}

testthat::test_that("common phenotype PC workflow exposes selected and all-PC outputs", {
    expect_contains_all(
        repo_path("workflows", "common", "calculate_phenotypePCs.wdl"),
        c(
            "File\\s+PhenotypePCsAllTSV\\s*=\\s*\"~\\{OutputPrefix\\}_phenotype_PCs~\\{OutputSuffix\\}\\.all\\.tsv\"",
            "File\\s+OutPhenotypePCsAll\\s*=\\s*ComputePCs\\.PhenotypePCsAllTSV",
            "File\\s+OutPhenotypePCs\\s*=\\s*ComputePCs\\.PhenotypePCsTSV"
        )
    )
})

testthat::test_that("public modality workflows pass through all-PC phenotype outputs", {
    workflow_expectations <- list(
        expression = list(
            path = repo_path("workflows", "expression", "prepare_eQTL.wdl"),
            patterns = c(
                "File\\s+IntPhenotypePCsAllOut\\s*=\\s*IntPhenotypePCs\\.OutPhenotypePCsAll",
                "File\\s+ScaledPhenotypePCsAllOut\\s*=\\s*ScaledPhenotypePCs\\.OutPhenotypePCsAll"
            )
        ),
        proteomics = list(
            path = repo_path("workflows", "proteomics", "prepare_pQTL.wdl"),
            patterns = c(
                "File\\s+IntPhenotypePCsAllOut\\s*=\\s*IntPhenotypePCs\\.OutPhenotypePCsAll",
                "File\\s+ScaledPhenotypePCsAllOut\\s*=\\s*ScaledPhenotypePCs\\.OutPhenotypePCsAll"
            )
        ),
        splicing = list(
            path = repo_path("workflows", "splicing", "prepare_sQTL.wdl"),
            patterns = c(
                "File\\s+IntPhenotypePCsAllOut\\s*=\\s*IntPhenotypePCs\\.OutPhenotypePCsAll",
                "File\\s+ScaledPhenotypePCsAllOut\\s*=\\s*ScaledPhenotypePCs\\.OutPhenotypePCsAll"
            )
        ),
        methylation_qtl_covariates = list(
            path = repo_path("workflows", "methylation", "qtl_covariates.wdl"),
            patterns = c(
                "File\\s+IntPhenotypePCsAllOut\\s*=\\s*IntPhenotypePCs\\.OutPhenotypePCsAll"
            )
        ),
        methylation_aggregate_arrays = list(
            path = repo_path("workflows", "methylation", "AggregateMethylationCohortArrays.wdl"),
            patterns = c(
                "File\\s+IntPhenotypePCsAllOut\\s*=\\s*IntPhenotypePCs\\.OutPhenotypePCsAll"
            )
        ),
        methylation_merge = list(
            path = repo_path("workflows", "methylation", "merge_methylation.wdl"),
            patterns = c(
                "File\\s+IntPhenotypePCsAllOut\\s*=\\s*CohortMerge\\.IntPhenotypePCsAllOut"
            )
        ),
        methylation_aggregate = list(
            path = repo_path("workflows", "methylation", "AggregateMethylationCohort.wdl"),
            patterns = c(
                "File\\s+IntPhenotypePCsAllOut\\s*=\\s*PrepareQtlCovariates\\.IntPhenotypePCsAllOut"
            )
        )
    )

    for (workflow in workflow_expectations) {
        expect_contains_all(workflow$path, workflow$patterns)
    }
})

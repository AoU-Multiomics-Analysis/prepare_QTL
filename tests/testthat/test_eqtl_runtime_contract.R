repo_path <- function(...) {
    testthat::test_path("..", "..", ...)
}

read_file <- function(path) {
    paste(readLines(path, warn = FALSE), collapse = "\n")
}

expect_runtime_defaults <- function(path, expected_count) {
    contents <- read_file(path)
    defaults <- c(
        'String\\s+DockerImage\\s*=\\s*"ghcr\\.io/aou-multiomics-analysis/prepare_qtl:main"',
        "Int\\s+preemptible_attempts\\s*=\\s*2",
        "Int\\s+max_retries\\s*=\\s*2"
    )

    for (pattern in defaults) {
        testthat::expect_equal(
            length(gregexpr(pattern, contents, perl = TRUE)[[1L]][
                gregexpr(pattern, contents, perl = TRUE)[[1L]] > 0L
            ]),
            expected_count,
            info = paste("Expected runtime default in", path, ":", pattern)
        )
    }
}

expect_task_runtime_and_logging <- function(path) {
    contents <- read_file(path)
    patterns <- c(
        "stage=",
        "start_time=",
        "completion_time=",
        "dimensions=",
        "outputs=",
        "docker:\\s*DockerImage",
        "preemptible:\\s*preemptible_attempts",
        "maxRetries:\\s*max_retries"
    )

    for (pattern in patterns) {
        testthat::expect_match(
            contents,
            pattern,
            perl = TRUE,
            info = paste("Missing runtime or logging contract in", path, ":", pattern)
        )
    }
}

testthat::test_that("eQTL runtime policy defaults are available to every public wrapper", {
    expect_runtime_defaults(
        repo_path("workflows", "expression", "prepare_eQTL.wdl"),
        expected_count = 2L
    )
    expect_runtime_defaults(
        repo_path("workflows", "common", "calculate_phenotypePCs.wdl"),
        expected_count = 2L
    )
    expect_runtime_defaults(
        repo_path("workflows", "common", "MergeCovariates.wdl"),
        expected_count = 2L
    )
    expect_runtime_defaults(
        repo_path("workflows", "common", "ResidualizePhenotypes.wdl"),
        expected_count = 1L
    )
})

testthat::test_that("eQTL tasks use the runtime policy and emit complete lifecycle logs", {
    workflow_paths <- c(
        repo_path("workflows", "expression", "prepare_eQTL.wdl"),
        repo_path("workflows", "common", "calculate_phenotypePCs.wdl"),
        repo_path("workflows", "common", "MergeCovariates.wdl"),
        repo_path("workflows", "common", "ResidualizePhenotypes.wdl")
    )

    for (path in workflow_paths) {
        expect_task_runtime_and_logging(path)
    }
})

testthat::test_that("eQTL preparation forwards its runtime policy through all calls", {
    contents <- read_file(repo_path("workflows", "expression", "prepare_eQTL.wdl"))
    for (input_name in c("DockerImage", "preemptible_attempts", "max_retries")) {
        pattern <- paste0(input_name, "\\s*=\\s*", input_name)
        testthat::expect_equal(
            length(gregexpr(pattern, contents, perl = TRUE)[[1L]][
                gregexpr(pattern, contents, perl = TRUE)[[1L]] > 0L
            ]),
            7L,
            info = paste("Expected runtime policy forwarding for", input_name)
        )
    }
})

testthat::test_that("common workflow wrappers forward their runtime policy to tasks", {
    workflow_paths <- c(
        repo_path("workflows", "common", "calculate_phenotypePCs.wdl"),
        repo_path("workflows", "common", "MergeCovariates.wdl")
    )

    for (path in workflow_paths) {
        contents <- read_file(path)
        for (input_name in c("DockerImage", "preemptible_attempts", "max_retries")) {
            testthat::expect_match(
                contents,
                paste0(input_name, "\\s*=\\s*", input_name),
                perl = TRUE,
                info = paste("Expected runtime policy forwarding in", path, "for", input_name)
            )
        }
    }
})

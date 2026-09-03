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

extract_table_lines <- function(lines, heading, header) {
    heading_index <- match(heading, lines)
    testthat::expect_false(is.na(heading_index))
    if (is.na(heading_index)) {
        return(character())
    }

    candidate_indices <- which(seq_along(lines) > heading_index & lines == header)
    testthat::expect_true(length(candidate_indices) > 0L)
    if (length(candidate_indices) == 0L) {
        return(character())
    }

    first_row <- candidate_indices[[1L]] + 2L
    remaining_lines <- lines[first_row:length(lines)]
    row_count <- match(FALSE, grepl("^\\|", remaining_lines), nomatch = 0L) - 1L
    if (row_count < 0L) {
        row_count <- length(remaining_lines)
    }
    remaining_lines[seq_len(row_count)]
}

parse_documented_input_contract <- function(lines, heading) {
    table_lines <- extract_table_lines(
        lines,
        heading,
        "| Input | Type | Default | Contract |"
    )
    if (length(table_lines) == 0L) {
        return(data.frame(
            name = character(),
            type = character(),
            default = character(),
            stringsAsFactors = FALSE
        ))
    }
    cells <- lapply(table_lines, function(line) {
        values <- trimws(strsplit(line, "|", fixed = TRUE)[[1L]])
        values[nzchar(values)][seq_len(3L)]
    })
    contract <- as.data.frame(
        do.call(rbind, cells),
        stringsAsFactors = FALSE
    )
    names(contract) <- c("name", "type", "default")
    contract[] <- lapply(contract, function(values) gsub("`", "", values, fixed = TRUE))
    contract
}

parse_expected_contract <- function(rows) {
    fields <- strsplit(rows, "|", fixed = TRUE)
    contract <- as.data.frame(do.call(rbind, fields), stringsAsFactors = FALSE)
    names(contract) <- c("name", "type", "default")
    contract
}

parse_wdl_input_contract <- function(path) {
    lines <- readLines(path, warn = FALSE)
    workflow_index <- grep("^workflow ", lines)[[1L]]
    input_index <- which(seq_along(lines) > workflow_index & lines == "  input {")[[1L]]
    end_index <- which(seq_along(lines) > input_index & lines == "  }")[[1L]]
    declarations <- lines[(input_index + 1L):(end_index - 1L)]
    declaration_pattern <- paste0(
        "^    ([A-Za-z]+(?:\\[[A-Za-z]+\\])?\\??) ",
        "([A-Za-z][A-Za-z0-9_]*)(?: = (.+))?$"
    )
    matches <- regexec(declaration_pattern, declarations, perl = TRUE)
    fields <- regmatches(declarations, matches)
    fields <- fields[lengths(fields) > 0L]

    contract <- do.call(rbind, lapply(fields, function(values) {
        type <- values[[2L]]
        default <- if (length(values) == 4L && nzchar(values[[4L]])) {
            values[[4L]]
        } else if (grepl("\\?$", type)) {
            "None"
        } else {
            "Required"
        }
        c(name = values[[3L]], type = type, default = default)
    }))
    as.data.frame(contract, stringsAsFactors = FALSE)
}

standalone_input_contract <- parse_expected_contract(c(
    "expression|File|Required",
    "gtf|File|Required",
    "lm22|File|Required",
    "precomputed_proportions|File?|None",
    "covariates|File?|None",
    "deconvolution_docker_image|String|\"ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression:main\"",
    "preemptible_attempts|Int|2",
    "max_retries|Int|2",
    "min_lm22_overlap|Float|0.80",
    "hspe_marker_fraction|Float|0.10",
    "hspe_quantile_normalize|Boolean|false",
    "hspe_batch_size|Int|100",
    "hspe_batch_memory|String|\"4 GB\"",
    "hspe_batch_disk_gb|Int|10",
    "group_mean_threshold|Float|0.0001",
    "zero_floor|Float|0.000001",
    "tca_max_iters|Int|10",
    "tca_parallel|Boolean|false",
    "random_seed|Int|20260901",
    "log2_pseudocount|Float|0.0",
    "gene_type|Array[String]|[\"protein_coding\", \"lncRNA\"]",
    "hspe_cpu|Int|4",
    "hspe_memory|String|\"32 GB\"",
    "hspe_disk_gb|Int|100",
    "proportions_cpu|Int|2",
    "proportions_memory|String|\"16 GB\"",
    "proportions_disk_gb|Int|50",
    "fit_cpu|Int|16",
    "fit_memory|String|\"256 GB\"",
    "fit_disk_gb|Int|750",
    "export_cpu|Int|8",
    "export_memory|String|\"256 GB\"",
    "export_disk_gb|Int|500",
    "manifest_cpu|Int|4",
    "manifest_memory|String|\"32 GB\"",
    "manifest_disk_gb|Int|100"
))

integrated_input_contract <- parse_expected_contract(c(
    "expression|File|Required",
    "gtf|File|Required",
    "lm22|File|Required",
    "precomputed_proportions|File?|None",
    "deconvolution_covariates|File?|None",
    "SampleList|File|Required",
    "AdditionalCovariates|File|Required",
    "OutputPrefix|String|Required",
    "deconvolution_docker_image|String|\"ghcr.io/aou-multiomics-analysis/prepare_qtl-cell-type-specific-expression:main\"",
    "qtl_docker_image|String|\"ghcr.io/aou-multiomics-analysis/prepare_qtl:main\"",
    "preemptible_attempts|Int|2",
    "max_retries|Int|2",
    "min_lm22_overlap|Float|0.80",
    "hspe_marker_fraction|Float|0.10",
    "hspe_quantile_normalize|Boolean|false",
    "hspe_batch_size|Int|100",
    "hspe_batch_memory|String|\"4 GB\"",
    "hspe_batch_disk_gb|Int|10",
    "group_mean_threshold|Float|0.0001",
    "zero_floor|Float|0.000001",
    "tca_max_iters|Int|10",
    "tca_parallel|Boolean|false",
    "random_seed|Int|20260901",
    "log2_pseudocount|Float|0.0",
    "gene_type|Array[String]|[\"protein_coding\", \"lncRNA\"]",
    "hspe_cpu|Int|4",
    "hspe_memory|String|\"32 GB\"",
    "hspe_disk_gb|Int|100",
    "proportions_cpu|Int|2",
    "proportions_memory|String|\"16 GB\"",
    "proportions_disk_gb|Int|50",
    "fit_cpu|Int|16",
    "fit_memory|String|\"256 GB\"",
    "fit_disk_gb|Int|750",
    "export_cpu|Int|8",
    "export_memory|String|\"256 GB\"",
    "export_disk_gb|Int|500",
    "manifest_cpu|Int|4",
    "manifest_memory|String|\"32 GB\"",
    "manifest_disk_gb|Int|100",
    "scatter_cpu|Int|1",
    "scatter_memory|String|\"4 GB\"",
    "scatter_disk_gb|Int|20",
    "eqtl_cpu|Int|8",
    "eqtl_memory|Int|64",
    "eqtl_disk_gb|Int|200"
))

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

testthat::test_that("guide and WDL publish the exact standalone input contract", {
    guide_lines <- readLines(
        repo_path("docs", "cell-type-specific-expression.md"),
        warn = FALSE
    )
    documented <- parse_documented_input_contract(
        guide_lines,
        "### Standalone deconvolution inputs"
    )
    implemented <- parse_wdl_input_contract(repo_path(
        "workflows", "cell_type_specific_expression", "deconvolution.wdl"
    ))

    testthat::expect_identical(documented, standalone_input_contract)
    testthat::expect_identical(implemented, standalone_input_contract)
})

testthat::test_that("guide and WDL publish the exact integrated input contract", {
    guide_lines <- readLines(
        repo_path("docs", "cell-type-specific-expression.md"),
        warn = FALSE
    )
    documented <- parse_documented_input_contract(
        guide_lines,
        "### Integrated cell-type eQTL inputs"
    )
    implemented <- parse_wdl_input_contract(repo_path(
        "workflows", "cell_type_specific_expression", "prepare_cell_type_eQTL.wdl"
    ))

    testthat::expect_identical(documented, integrated_input_contract)
    testthat::expect_identical(implemented, integrated_input_contract)
})

testthat::test_that("guide states the expression and proportion input contracts", {
    guide <- read_document(
        repo_path("docs", "cell-type-specific-expression.md")
    )

    expect_documented(guide, c(
        "linear CPM",
        "first four columns.*#chr.*start.*end.*gene_id",
        "at least one sample column",
        "modeled.*nonconstant gene.*preserves.*coordinate.*order",
        "finite.*nonnegative",
        "log2_pseudocount.*defaults to.*0",
        "[Zz]ero bulk CPM.*positive pseudocount",
        "TCA accepts.*zero CPM",
        "[Nn]egative.*always.*invalid",
        "hspe uses.*log2\\(CPM \\+ log2_pseudocount\\)",
        "log2\\(LM22 \\+ log2_pseudocount\\)",
        "TCA uses.*linear CPM",
        "gene_type.*protein_coding.*lncRNA",
        "gene-type filtering.*before.*hspe.*TCA",
        "LM22.*hspe",
        "precomputed LM22 proportions",
        "LM22.*required",
        "precomputed_proportions.*provided.*skip.*hspe",
        "exactly the 22 standard LM22 columns",
        "sample_id.*first column",
        "sample IDs.*unique.*non-empty",
        "values.*finite.*nonnegative",
        "rows.*sum to one.*1e-6",
        "sample row order.*expression BED sample-column order",
        "does not reorder samples",
        "fewer than two major groups.*group_mean_threshold.*fails",
        "zero_floor.*exact zero values",
        "preemptible_attempts.*2",
        "max_retries.*2"
    ))

    lm22_cell_types <- c(
        "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
        "T cells CD4 naive", "T cells CD4 memory resting",
        "T cells CD4 memory activated", "T cells follicular helper",
        "T cells regulatory (Tregs)", "T cells gamma delta",
        "NK cells resting", "NK cells activated", "Monocytes",
        "Macrophages M0", "Macrophages M1", "Macrophages M2",
        "Dendritic cells resting", "Dendritic cells activated",
        "Mast cells resting", "Mast cells activated", "Eosinophils",
        "Neutrophils"
    )
    for (cell_type in lm22_cell_types) {
        testthat::expect_match(guide, cell_type, fixed = TRUE)
    }
})

testthat::test_that("guide states the scattered QTL preparation contract", {
    guide <- read_document(
        repo_path("docs", "cell-type-specific-expression.md")
    )

    expect_documented(guide, c(
        "Log2CpmBed",
        "OutputPrefix.*start.*ASCII letter or number.*only letters, numbers, dots, underscores, and hyphens",
        "skips count.*normalization.*transformation",
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

    table_lines <- extract_table_lines(
        guide_lines,
        "## QTL manifest",
        manifest_header
    )
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

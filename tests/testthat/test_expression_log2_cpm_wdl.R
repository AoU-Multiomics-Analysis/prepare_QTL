testthat::test_that("the eQTL workflow declares optional raw-count and log2 CPM inputs", {
    miniwdl <- Sys.which("miniwdl")
    testthat::skip_if(!nzchar(miniwdl), "miniwdl is not installed")
    workflow <- testthat::test_path(
        "..", "..", "workflows", "expression", "prepare_eQTL.wdl"
    )
    miniwdl_python <- sub("^#!", "", readLines(miniwdl, n = 1L, warn = FALSE))
    inspector_path <- tempfile(fileext = ".py")
    template_path <- tempfile(fileext = ".json")
    writeLines(c(
        "import json",
        "import sys",
        "import WDL",
        "document = WDL.load(sys.argv[1])",
        "inputs = {binding.name: str(binding.value.type) for binding in document.workflow.available_inputs}",
        "print(json.dumps(inputs, sort_keys=True))"
    ), inspector_path)
    status <- system2(
        miniwdl_python,
        c(inspector_path, workflow),
        stdout = template_path,
        stderr = FALSE
    )
    testthat::expect_identical(status, 0L)

    template <- jsonlite::read_json(template_path, simplifyVector = TRUE)
    testthat::expect_identical(template[["Log2CpmBed"]], "File?")
    testthat::expect_identical(template[["CountGCT"]], "File?")
    testthat::expect_identical(template[["AnnotationGTF"]], "File?")
})

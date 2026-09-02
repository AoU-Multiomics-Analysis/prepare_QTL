pipeline_root <- normalizePath(file.path(testthat::test_path(), "../../.."))
script_root <- file.path(
  pipeline_root,
  "scripts",
  "cell_type_specific_expression"
)
Sys.setenv(CELL_TYPE_SPECIFIC_EXPRESSION_ROOT = script_root)
r_files <- file.path(script_root, "R", c("constants.R", "io.R", "integration.R"))
invisible(lapply(r_files[file.exists(r_files)], source, local = .GlobalEnv))

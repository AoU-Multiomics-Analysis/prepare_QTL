script_root <- Sys.getenv(
  "CELL_TYPE_SPECIFIC_EXPRESSION_ROOT",
  unset = "/opt/prepare_qtl/scripts/cell_type_specific_expression"
)
# Modules load from scripts/cell_type_specific_expression/R.
module_root <- file.path(script_root, "R")
r_files <- list.files(module_root, pattern = "[.]R$", full.names = TRUE)
invisible(lapply(sort(r_files), source))

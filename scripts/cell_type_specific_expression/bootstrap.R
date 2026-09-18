script_root <- Sys.getenv(
  "CELL_TYPE_SPECIFIC_EXPRESSION_ROOT",
  unset = "/opt/prepare_qtl/scripts/cell_type_specific_expression"
)
# Modules live in scripts/cell_type_specific_expression/R.
module_root <- file.path(script_root, "R")
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
entry_name <- if (length(file_arg)) {
  basename(sub("^--file=", "", file_arg[[1L]]))
} else {
  ""
}
module_manifest <- file.path(script_root, "modules", sub("[.]R$", ".txt", entry_name))
if (nzchar(entry_name) && file.exists(module_manifest)) {
  modules <- readLines(module_manifest, warn = FALSE)
  if (!length(modules) || any(!grepl("^[A-Za-z0-9_]+[.]R$", modules)) || anyDuplicated(modules)) {
    stop("Invalid task module manifest: ", module_manifest, call. = FALSE)
  }
  r_files <- file.path(module_root, modules)
  if (any(!file.exists(r_files))) stop("Missing task helper module", call. = FALSE)
} else {
  # Tests and development utilities may load the complete module library.
  entry_parent <- if (length(file_arg)) basename(dirname(sub("^--file=", "", file_arg[[1L]]))) else ""
  if (entry_parent %in% c("estimation", "fit", "export", "downstream")) {
    stop("Missing task module manifest: ", module_manifest, call. = FALSE)
  }
  r_files <- list.files(module_root, pattern = "[.]R$", full.names = TRUE)
}
invisible(lapply(sort(r_files), source))

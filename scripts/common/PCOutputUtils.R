format_pca_outputs <- function(rotated, n_pcs) {
    rotated <- as.data.frame(rotated, check.names = FALSE)
    pc_names <- names(rotated)
    all_output <- cbind(
        ID = sub("X", "", rownames(rotated), fixed = TRUE),
        rotated
    )
    selected_output <- all_output[, c("ID", pc_names[seq_len(n_pcs)]), drop = FALSE]
    list(selected = selected_output, all = all_output)
}

write_pca_outputs <- function(rotated, n_pcs, selected_path, all_path) {
    outputs <- format_pca_outputs(rotated, n_pcs)
    readr::write_tsv(outputs$selected, selected_path)
    readr::write_tsv(outputs$all, all_path)
    invisible(outputs)
}

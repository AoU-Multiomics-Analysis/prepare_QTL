args <- commandArgs(trailingOnly = TRUE)
if (length(args) %% 2 != 0) stop('Arguments must be named value pairs')
opts <- setNames(as.list(args[seq(2, length(args), 2)]), sub('^--', '', args[seq(1, length(args), 2)]))
required <- c('library', 'counts-dir', 'capture-plan', 'selected', 'feature-map', 'monaco', 'hpca', 'out-dir', 'threads')
if (!all(required %in% names(opts))) stop('Required named arguments are missing')
.libPaths(c(normalizePath(opts$library, mustWork = TRUE), .libPaths()))
suppressPackageStartupMessages({library(SingleR); library(Matrix); library(dplyr); library(readr); library(tibble)})
log_msg <- function(...) message(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), ' ', ...)
for (key in c('capture-plan', 'selected', 'feature-map', 'monaco', 'hpca')) {
  if (!file.exists(opts[[key]]) || file.access(opts[[key]], 4) != 0) stop('Unreadable input: ', key)
}
out <- opts[['out-dir']]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
threads <- as.integer(opts$threads)
stopifnot(threads > 0)
set.seed(20261002)
selected <- read_tsv(opts$selected, show_col_types = FALSE)
stopifnot(!anyDuplicated(selected$cell_id), nrow(selected) == 65191, all(selected$scDblFinder.class == 'singlet'))
features <- read_tsv(opts[['feature-map']], show_col_types = FALSE) |>
  mutate(exclude_sex = tolower(as.character(exclude_sex)) == 'true')
stopifnot(!anyDuplicated(features$feature_id), !anyNA(features$exclude_sex))
symbols <- sort(unique(features$feature_symbol[!features$exclude_sex & !is.na(features$feature_symbol) & nzchar(features$feature_symbol)]))
monaco <- readRDS(opts$monaco)
hpca <- readRDS(opts$hpca)
references <- list(Monaco_main = list(ref = monaco, labels = monaco$label.main),
                   Monaco_fine = list(ref = monaco, labels = monaco$label.fine),
                   HPCA_main = list(ref = hpca, labels = hpca$label.main))
trained <- list()
training_stats <- list()
for (model in names(references)) {
  log_msg('Training ', model)
  ref <- references[[model]]$ref
  overlap <- intersect(rownames(ref), symbols)
  stopifnot(length(overlap) > 5000)
  trained[[model]] <- trainSingleR(ref = ref, labels = references[[model]]$labels,
    restrict = overlap, de.method = 'classic', num.threads = threads)
  training_stats[[model]] <- tibble(model = model, shared_genes = length(overlap),
    reference_samples = ncol(ref), reference_labels = length(unique(references[[model]]$labels)))
}
write_tsv(bind_rows(training_stats), file.path(out, 'training_stats.tsv'))
write_tsv(as_tibble(as.data.frame(colData(monaco))) |> count(label.main, label.fine), file.path(out, 'Monaco_label_map.tsv'))
predictions <- setNames(lapply(names(trained), function(x) list()), names(trained))
capture_stats <- list()
marker_rows <- list()
panel <- c('CLC', 'HDC', 'FCER1A', 'GATA2', 'CPA3', 'MS4A2', 'CCR3', 'S100A8', 'S100A9', 'S100A12',
 'TRAC', 'CD3D', 'CD3E', 'CD3G', 'CD247', 'CD8A', 'CD8B', 'NKG7', 'GNLY', 'KLRD1', 'LYZ', 'LST1',
 'FCN1', 'TYROBP', 'FCER1G', 'CSF3R', 'FCGR3B', 'CEACAM8', 'MS4A1', 'CD79A', 'MZB1', 'JCHAIN', 'HBB', 'HBA1', 'PPBP', 'PF4')
captures <- read_tsv(opts[['capture-plan']], show_col_types = FALSE)
capture_col <- intersect(c('capture', 'capture_id', '10X_run'), names(captures))
if (length(capture_col) != 1) stop('Capture plan must have one capture ID field')
for (capture in captures[[capture_col]]) {
  log_msg('Loading ', capture)
  root <- file.path(opts[['counts-dir']], capture)
  paths <- file.path(root, c('counts.mtx.gz', 'features.tsv', 'cells.tsv'))
  if (!all(file.exists(paths)) || any(file.access(paths, 4) != 0)) stop('Unreadable capture inputs: ', capture)
  cells <- read_tsv(paths[3], show_col_types = FALSE)
  ids <- read_tsv(paths[2], show_col_types = FALSE)$feature_id
  con <- gzfile(paths[1], open = 'rt')
  counts <- as(readMM(con), 'CsparseMatrix')
  close(con)
  stopifnot(nrow(counts) == length(ids), ncol(counts) == nrow(cells), !anyDuplicated(cells$cell_id), !anyDuplicated(ids))
  chosen <- which(cells$cell_id %in% selected$cell_id)
  stopifnot(length(chosen) > 0)
  counts <- counts[, chosen, drop = FALSE]
  cell_ids <- cells$cell_id[chosen]
  gene_map <- features[match(ids, features$feature_id), ]
  stopifnot(!anyNA(gene_map$feature_id))
  valid <- which(!gene_map$exclude_sex & gene_map$feature_symbol %in% symbols)
  aggregation <- sparseMatrix(i = match(gene_map$feature_symbol[valid], symbols), j = valid,
    x = 1, dims = c(length(symbols), length(ids)))
  test <- aggregation %*% counts
  rownames(test) <- symbols
  colnames(test) <- cell_ids
  stopifnot(all(colSums(test) > 0), all(is.finite(test@x)))
  marker <- as_tibble(t(as.matrix(test[match(panel, symbols), , drop = FALSE])), .name_repair = 'minimal')
  names(marker) <- paste0(panel, '_UMI')
  marker_rows[[capture]] <- bind_cols(tibble(cell_id = cell_ids, nonsex_raw_total_UMI = as.numeric(colSums(test))), marker)
  for (model in names(trained)) {
    log_msg('Classifying ', length(cell_ids), ' cells: ', capture, ' / ', model)
    pred <- classifySingleR(test = test, trained = trained[[model]], prune = FALSE, num.threads = threads)
    stopifnot(identical(rownames(pred), cell_ids), all(is.finite(pred$scores)))
    predictions[[model]][[capture]] <- pred
  }
  capture_stats[[capture]] <- tibble(capture = capture, exported_cells = nrow(cells),
    reviewed_cells = length(cell_ids), raw_features = length(ids), excluded_sex_features = sum(gene_map$exclude_sex),
    summed_nonsex_UMI = sum(test), summed_sex_UMI = sum(counts[gene_map$exclude_sex, , drop = FALSE]))
  rm(counts, test, aggregation)
  gc()
}
write_tsv(bind_rows(capture_stats), file.path(out, 'capture_stats.tsv'))
markers <- bind_rows(marker_rows)
stopifnot(!anyDuplicated(markers$cell_id), setequal(markers$cell_id, selected$cell_id))
markers <- markers[match(selected$cell_id, markers$cell_id), ]
write_tsv(markers, file.path(out, 'raw_marker_counts.tsv.gz'))
all_calls <- selected
for (model in names(predictions)) {
  log_msg('Applying global score pruning: ', model)
  pred <- do.call(rbind, predictions[[model]])
  stopifnot(!anyDuplicated(rownames(pred)), setequal(rownames(pred), selected$cell_id))
  pred <- pred[match(selected$cell_id, rownames(pred)), ]
  discard <- pruneScores(pred)
  pred$pruned.labels <- pred$labels
  pred$pruned.labels[discard] <- NA_character_
  saveRDS(pred, file.path(out, paste0(model, '_predictions.rds')))
  calls <- tibble(label = pred$labels, pruned_label = pred$pruned.labels,
    delta_next = pred$delta.next, delta_median = getDeltaFromMedian(pred))
  names(calls) <- paste0(model, '_', names(calls))
  all_calls <- bind_cols(all_calls, calls)
  scores <- bind_cols(tibble(cell_id = selected$cell_id), as_tibble(pred$scores, .name_repair = 'minimal'))
  write_tsv(scores, file.path(out, paste0(model, '_scores.tsv.gz')))
}
write_tsv(all_calls, file.path(out, 'cell_predictions.tsv.gz'))
writeLines(capture.output(sessionInfo()), file.path(out, 'session_info.txt'))
log_msg('Review complete: ', nrow(all_calls), ' retained cells; current reference labels were not changed')

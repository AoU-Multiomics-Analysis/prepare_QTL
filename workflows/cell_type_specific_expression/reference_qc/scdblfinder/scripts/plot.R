suppressPackageStartupMessages(library(tidyverse))
args <- commandArgs(trailingOnly = TRUE)
data <- read_tsv(file.path(args[[1]], 'reference_cell_type_summary.tsv'), show_col_types = FALSE)
p <- ggplot(data, aes(reorder(matrix_header, doublet_fraction), doublet_fraction)) +
  geom_col(fill = '#426B90', width = 0.75) + coord_flip() +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = NULL, y = 'Predicted doublets') + theme_classic(base_size = 11)
ggsave(file.path(args[[1]], 'predicted_doublets_by_cell_type.png'), p, width = 7, height = 4.5, dpi = 160)
ggsave(file.path(args[[1]], 'predicted_doublets_by_cell_type.pdf'), p, width = 7, height = 4.5)

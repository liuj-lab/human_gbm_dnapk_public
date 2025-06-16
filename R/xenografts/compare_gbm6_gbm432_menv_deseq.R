# Compare the GBM6 and GBM43_v2 DESeq outputs for Figure 2 of the manuscript.
#
# Create for manuscript preparation on 5/5/25. Last updated 5/5/25.

library(Seurat)
library(dplyr)
library(SCpubr)
library(ggplot2)
library(reticulate)
library(EnhancedVolcano)
library(orthogene)
library(clusterProfiler)
library(tidyr)
library(ComplexHeatmap)
library(colorRamp2)

output_dir <- "output/xenografts/playground/compare_gbm6_gbm432_menv_deseq"

utils.deseq <- new.env()
source("R/utils/deseq.R", local = utils.deseq)

# =========================================================
# Import data
# =========================================================
# Import DESeq data. Specifically, we have a bunch of different samples here:
# GBM6:
#   Macrophages + microglia
#   Macrophages + microglia + NK
# GBM432:
#   Macrophages + microglia

gbm6_macrophages_microglia <- utils.deseq$extract_deseq_list(
  paste0(
    "output/xenografts/playground/gbm6_tumor_menv_clustering/",
    "gbm6_immune/macrophage_microglia/tables"
  ),
  delim = " ",
  header = TRUE
)
gbm6_macrophages_microglia_nk <- utils.deseq$extract_deseq_list(
  paste0(
    "output/xenografts/playground/gbm6_tumor_menv_clustering/",
    "gbm6_immune/macrophage_microglia_nk/tables"
  ),
  delim = "\t",
  header = TRUE
)
gbm432_macrophages_microglia <- utils.deseq$extract_deseq_list(
  paste0(
    "output/xenografts/playground/gbm432_tumor_menv_clustering/",
    "gbm432_5_immune/macrophages_microglias/tables"
  ),
  delim = "\t",
  header = TRUE
)
deseq_raw <- list(
  gbm6_macrophages_microglia = gbm6_macrophages_microglia,
  gbm6_macrophages_microglia_nk = gbm6_macrophages_microglia_nk,
  gbm432_macrophages_microglia = gbm432_macrophages_microglia
)
raw_deseq_dfs <- unlist(deseq_raw, recursive = FALSE)

# Get a list of interferon-related genes
mouse_interferon_genes <- read.gmt(
  "data/external_gene_lists/REACTOME_INTERFERON_SIGNALING.v2024.1.Mm.gmt"
)$gene
human_interferon_genes <- read.gmt(
  "data/external_gene_lists/REACTOME_INTERFERON_SIGNALING.v2023.2.Hs.gmt"
)$gene

# Custom gene set
human_interferon_subset <- c(
  "TRIM22", "IFITM1", "MX1", "ISG20", "RSAD2", "OAS2", "IFI6", "OASL",
  "MX2", "USP18", "ISG15", "STAT2", "SOCS3", "EGR1", "SP100", "OAS3",
  "IRF7", "STAT1", "TRIM5", "TRIM25", "ADAR", "OAS1", "SAMHD1", "EIF2AK2",
  "UBE2L6", "TRIM45", "TRIM21", "HLA-F", "HLA-C", "HLA-A", "STAT3", "XAF1",
  "IFIT3", "IFIT1", "HLA-B", "PRKCD", "VCAM1", "HLA-E", "GBP4", "HERC5",
  "IFI27", "TRIM14", "SPHK1", "TUBB4A", "PTPN6", "TUBA8", "SNCA", "MAPT",
  "PSMB8", "RNASEL", "GBP3", "MT2A", "CAMK2A", "BST2", "IFI35", "IFIT5",
  "IFITM3", "IFITM2", "SOCS1"
)

# Convert the human genes to their mouse orthologs
mouse_ortholog_human_interferon_genes <- orthogene::convert_orthologs(
  human_interferon_genes,
  input_species = "human", output_species = "mouse"
)
mouse_ortholog_human_interferon_genes <- rownames(
  mouse_ortholog_human_interferon_genes
)
mouse_ortholog_human_subset <- orthogene::convert_orthologs(
  human_interferon_subset,
  input_species = "human", output_species = "mouse"
)
mouse_ortholog_human_subset <- rownames(mouse_ortholog_human_subset)

# Filter the DESeq outputs to include only these genes
mouse_gene_filter <- function(df) {
  rownames(df) <- df$feature
  df <- df %>% filter(feature %in% mouse_interferon_genes)
}
mouse_ortholog_gene_filter <- function(df) {
  rownames(df) <- df$feature
  df <- df %>% filter(feature %in% mouse_ortholog_human_interferon_genes)
}
mouse_subset_gene_filter <- function(df) {
  rownames(df) <- df$feature
  df <- df %>% filter(feature %in% mouse_ortholog_human_subset)
}
mouse_deseq_dfs <- lapply(raw_deseq_dfs, mouse_gene_filter)
mouse_ortholog_deseq_dfs <- lapply(raw_deseq_dfs, mouse_ortholog_gene_filter)
mouse_subset_deseq_dfs <- lapply(raw_deseq_dfs, mouse_subset_gene_filter)

# =========================================================
# Plotting
# =========================================================
deseq_dfs <- mouse_subset_deseq_dfs

# Get LFC matrices
lfc_deseq_dfs <- lapply(
  names(deseq_dfs), function(name) {
    df <- deseq_dfs[[name]]
    df <- df %>% dplyr::select(feature, log_fc)
    colnames(df)[colnames(df) == "log_fc"] <- paste0("log_fc_", name)
    return(df)
  }
)
lfc_df <- Reduce(
  function(x, y) merge(x, y, by = "feature", all = TRUE),
  lfc_deseq_dfs
)
rownames(lfc_df) <- lfc_df$feature
lfc_df <- lfc_df[, -1]
lfc_df[is.na(lfc_df)] <- 0
colnames(lfc_df) <- gsub("_untreated_noRT", "", colnames(lfc_df))

ht <- Heatmap(
  as.matrix(lfc_df),
  col = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red")),
  name = "LFC across cell conditions",
  width = unit(10, "cm"),
  height = unit(28, "cm"),
  show_row_names = TRUE
)
pdf(
  file.path(
    output_dir,
    "heatmaps/mouse_subset_ortholog_lfc_heatmap.pdf"
  ),
  width = 8, height = 20
)
drawn_heatmap <- draw(ht)
dev.off()

# Make a bubble plot to show significance

# Get a padj matrix
padj_dfs <- lapply(names(deseq_dfs), function(name) {
  df <- deseq_dfs[[name]]
  df <- df %>% dplyr::select(feature, padj)
  colnames(df)[colnames(df) == "padj"] <- paste0("padj_", name)
  return(df)
})
padj_df <- Reduce(function(x, y) merge(x, y, by = "feature", all = TRUE), padj_dfs)
rownames(padj_df) <- padj_df$feature
padj_df[is.na(padj_df)] <- 1 # Assign NA padj values such that they effectively disappear
colnames(padj_df) <- gsub("_untreated_noRT", "", colnames(padj_df))

# Select out the interferon genes and columns that we want
lfc_df$feature <- rownames(lfc_df)
log_fc_long <- lfc_df %>%
  pivot_longer(cols = starts_with("log_fc"), names_to = "experiment", values_to = "log_fc")
log_fc_long$experiment <- gsub("log_fc_", "", log_fc_long$experiment)
padj_long <- padj_df %>%
  pivot_longer(cols = starts_with("padj"), names_to = "experiment", values_to = "padj")
padj_long$experiment <- gsub("padj_", "", padj_long$experiment)
plot_data <- left_join(log_fc_long, padj_long, by = c("feature", "experiment"))
plot_data <- plot_data %>%
  mutate(size = -log10(padj))

# Fix some things before plotting

# Sort the rows based on the ordering in the heatmap, which is clustered
gene_order <- rownames(lfc_df)[row_order(drawn_heatmap)]
plot_data$feature <- factor(plot_data$feature, levels = gene_order)
plot_data <- plot_data %>% arrange(feature)

# Sort the columns based on the column order of the heatmap
column_order <- c(
  "gbm6_macrophages_microglia.untreated_RT",
  "gbm6_macrophages_microglia.nedisertib_noRT",
  "gbm6_macrophages_microglia.nedisertib_RT",
  "gbm6_macrophages_microglia_nk.untreated_RT",
  "gbm6_macrophages_microglia_nk.nedisertib_noRT",
  "gbm6_macrophages_microglia_nk.nedisertib_RT",
  "gbm432_macrophages_microglia.untreated_RT",
  "gbm432_macrophages_microglia.nedisertib_noRT",
  "gbm432_macrophages_microglia.nedisertib_RT"
)
plot_data$experiment <- factor(plot_data$experiment, levels = column_order)
plot_data <- plot_data %>% arrange(experiment)

# Take all padj values and cap them at 40, which has the effect of making the
# Inf values show up
plot_data$size <- pmin(plot_data$size, 40)
plot_data$log_fc_clipped <- pmax(pmin(plot_data$log_fc, 2), -2)
bubble_plot <- ggplot(plot_data, aes(x = experiment)) +
  geom_point(aes(y = feature, size = size, color = log_fc_clipped)) +
  geom_point(
    data = subset(plot_data, padj < 0.05),
    aes(y = feature, size = size, shape = "padj < 0.05"),
    color = "black", fill = NA, alpha = 0.8
  ) +
  scale_size(range = c(1, 5)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
  ) +
  scale_shape_manual(
    name = "Significance",
    values = c(`padj < 0.05` = 21),
    guide = guide_legend(override.aes = list(color = "black", fill = NA))
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.key.height = unit(0.3, "cm"),
    axis.line = element_line(colour = "black"),
    legend.key = element_blank()
  ) +
  labs(
    title = "Mouse Interferon Subset Gene Set",
    x = "Experiment",
    y = "Feature",
    size = "-log10(padj), capped at 40",
    color = "log_fc"
  )
ggsave(file.path(output_dir, "mouse_subset_interferon_bubble_plot.pdf"),
  device = "pdf", height = 9, width = 4
)
ggsave(file.path(output_dir, "mouse_subset_interferon_bubble_plot.png"),
  device = "png", height = 9, width = 4
)

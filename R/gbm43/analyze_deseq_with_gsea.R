# Analyze DESeq outputs with gene set analysis. This includes enrichment
# and analysis of gene sets with various heatmaps
# This pipeline uses independent filtering, separates positive and negative
# signals, and generates leading edge gene set summaries.

# Note that this was generated on 3/14/25 and supersedes previous GSEA implementations.
# All other GSEA implementations are considered deprecated.

# =========================================================
# Library imports and constants
# =========================================================

library(fgsea)
library(msigdbr)
library(dplyr)
library(parallel)
library(patchwork)
library(purrr)
library(ggplot2)
library(colorRamp2)
library(ComplexHeatmap)
library(Seurat)
library(tidyverse)

utils.deseq <- new.env()
source("R/utils/deseq.R", local = utils.deseq)

utils.gsea <- new.env()
source("R/utils/gsea.R", local = utils.gsea)

OUTPUT_DIR <- "output/gbm43/gsea/noRTNormalized"
DESEQ_DIR <- "output/gbm43/deseq/noRTNormalized"
DATA_DIR <- "data/gbm43_perturb_seq"

# =========================================================
# Import DESeq data
# =========================================================
deseq_data <- utils.deseq$extract_deseq_list(
  DESEQ_DIR
)

# Process all to have gene symbols
deseq_data <- lapply(deseq_data, function(df) {
  df$feature <- gsub("GRCh38-", "", df$feature)
  return(df)
})

# Remove NAs adjusted p-values, which accounts for independent filtering
deseq_data <- lapply(deseq_data, function(df) {
  df <- df[!is.na(df$log_fc), ]
  df <- df[!is.na(df$padj), ]
  return(df)
})

significant_deseq_data <- lapply(deseq_data, function(df) {
  df <- df[df$pvalue < 0.05, ]
})

# =========================================================
# Generate average and median LFC matrices over known gene sets
# and plot them as heatmaps
# =========================================================
msigdbr_df <- msigdbr(species = "Homo sapiens", category = "H")
msigdbr_gene_set_lists <- split(msigdbr_df$gene_symbol, msigdbr_df$gs_name)

# Remove perturbations that have less than 100 significant measurements.
significant_deseq_data <- Filter(function(df) nrow(df) > 100, significant_deseq_data)

# Generate a blank matrix to store values in
generate_blank_matrix <- function() {
  gene_set_names <- names(msigdbr_gene_set_lists)
  perturb_names <- names(significant_deseq_data)
  result_matrix <- matrix(nrow = length(gene_set_names), ncol = length(perturb_names))
  rownames(result_matrix) <- gene_set_names
  colnames(result_matrix) <- perturb_names
  return(result_matrix)
}

# Average LFC over gene sets and plot heatmap. Plot only those perturbations
# that perturb more than 20% of a given gene set.
mean_matrix <- generate_blank_matrix()
count_matrix <- generate_blank_matrix()
for (i in seq_along(msigdbr_gene_set_lists)) {
  gene_set <- msigdbr_gene_set_lists[[i]]
  for (j in seq_along(significant_deseq_data)) {
    perturb_df <- significant_deseq_data[[j]]
    perturb_df <- perturb_df[perturb_df$feature %in% gene_set, ]
    count_matrix[i, j] <- nrow(perturb_df)
    if (nrow(perturb_df) > 0.2 * length(gene_set)) {
      avg_log_fc <- mean(perturb_df$log_fc, na.rm = TRUE)
      mean_matrix[i, j] <- avg_log_fc
    } else {
      mean_matrix[i, j] <- 0
    }
  }
}
heatmap <- Heatmap(
  mean_matrix,
  name = "Mean LFC for overrepresented gene sets",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  column_title = "Perturbs",
  row_title = "Pathways",
  show_column_names = TRUE,
  show_row_names = TRUE,
  cluster_rows = TRUE, # Enable clustering for rows
  cluster_columns = TRUE, # Enable clustering for columns
  row_names_gp = gpar(fontsize = 7), # Adjust row name font size
  column_names_gp = gpar(fontsize = 8), # Adjust column name font size
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)
pdf(file.path(OUTPUT_DIR, "mean_LFC_heatmap.pdf"), width = 20, height = 10)
draw(heatmap, heatmap_legend_side = "top")
dev.off()

# Median LFC over gene set
median_matrix <- generate_blank_matrix()
for (i in seq_along(msigdbr_gene_set_lists)) {
  gene_set <- msigdbr_gene_set_lists[[i]]
  for (j in seq_along(significant_deseq_data)) {
    perturb_df <- significant_deseq_data[[j]]
    perturb_df <- perturb_df[perturb_df$feature %in% gene_set, ]
    median_lfc <- median(perturb_df$log_fc, na.rm = TRUE)
    if (nrow(perturb_df) > 0.2 * length(gene_set)) {
      median_log_fc <- median(perturb_df$log_fc, na.rm = TRUE)
      median_matrix[i, j] <- median_log_fc
    } else {
      median_matrix[i, j] <- 0
    }
  }
}
heatmap <- Heatmap(
  median_matrix,
  name = "Median LFC for overrepresented gene sets",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  column_title = "Perturbs",
  row_title = "Pathways",
  show_column_names = TRUE,
  show_row_names = TRUE,
  cluster_rows = TRUE, # Enable clustering for rows
  cluster_columns = TRUE, # Enable clustering for columns
  row_names_gp = gpar(fontsize = 7), # Adjust row name font size
  column_names_gp = gpar(fontsize = 8), # Adjust column name font size
  heatmap_legend_param = list(
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)
pdf(file.path(OUTPUT_DIR, "median_LFC_heatmap.pdf"), width = 20, height = 10)
draw(heatmap, heatmap_legend_side = "top")
dev.off()

# =========================================================
# Generate GSEA results
# =========================================================
# Genreate ranked list
deseq_ranked_lists <- lapply(deseq_data, function(df) {
  set.seed(5220)
  df <- df %>% mutate(rank = rank(log_fc, ties.method = "random"))
  df <- df[order(-df$rank), ]
  lfc_list <- df$log_fc
  gene_names <- df$feature
  names(lfc_list) <- gene_names
  return(lfc_list)
})

# Generate GSEA
withWarnings <- function(expr) {
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(value = val, warnings = myWarnings)
}

# Run fgsea and store warnings
fgsea_output_list <- lapply(deseq_ranked_lists, function(ranked_list) {
  set.seed(42)
  output <- list()
  result_with_warnings <- withWarnings(
    fgsea(
      pathways = msigdbr_gene_set_lists,
      stats = ranked_list,
      maxSize = 500,
      eps = 0.0,
      nPermSimple = 10000,
      scoreType = "std"
    )
  )
  output[["Results"]] <- result_with_warnings$value
  output[["Warnings"]] <- result_with_warnings$warnings
  return(output)
})
saveRDS(fgsea_output_list, file.path(OUTPUT_DIR, "fgsea_std_output_list.rds"))
fgsea_results_list <- lapply(fgsea_output_list, function(output) {
  return(output$Results)
})
fgsea_warnings_list <- lapply(fgsea_output_list, function(output) {
  return(output$Warnings)
})

# =========================================================
# Prepare for GSEA plotting
# =========================================================
# Retrieve knockdown values. Average these values across the samples.
kd_data <- read.table(file.path(DATA_DIR, "gbm43_remainingRNATbl.txt"))
noRT_kds <- rowMeans(kd_data[, grepl("noRT", colnames(kd_data))])
RT_kds <- rowMeans(kd_data[, !grepl("noRT", colnames(kd_data))])
names(noRT_kds) <- paste(names(noRT_kds), "_noRT", sep = "")
names(RT_kds) <- paste(names(RT_kds), "_RT", sep = "")
kd_map <- c(noRT_kds, RT_kds)
saveRDS(kd_map, file.path(DATA_DIR, "gbm43_sample_averaged_kd_map.rds"))

# Calculate dropout values by calculating LFC against the same condition
# non-targeting control. These should be controlled per sample
raw_seurat_data <- readRDS("/raleighlab/data1/liuj/gbm_perturb/analysis/GBM43_1_malignant_only_annotated_20230817.Rds")

# For each FACS-sorted sample, calculate the ratio from the
# number of non-targeting guides for each. We could also calculate p-values
# for these based on idk...a t-test?
metadata_table <- raw_seurat_data@meta.data
dropout_table <- matrix(
  data = 0, nrow = length(unique(metadata_table$orig.ident)),
  ncol = length(unique(metadata_table$sgRNA))
)
rownames(dropout_table) <- unique(metadata_table$orig.ident)
colnames(dropout_table) <- unique(metadata_table$sgRNA)
for (sample in unique(raw_seurat_data$orig.ident)) {
  # Retrieve the portion of the dataframe with that sample
  df <- metadata_table %>% dplyr::filter(orig.ident == sample)
  guide_counts <- table(df$sgRNA)
  for (sgRNA in names(guide_counts)) {
    dropout_table[sample, sgRNA] <- guide_counts[sgRNA] / guide_counts["non-targeting"]
  }
}
write.table(dropout_table, file.path(DATA_DIR, "gbm43_dropout_ratios_per_sample.txt"))
dropout_table <- read.table(file.path(DATA_DIR, "gbm43_dropout_ratios_per_sample.txt"))

# Impute 0s with the minimum detection of the sample (conservative wrt effect),
# then find log fold changes and average across RT and noRT samples
dropout_table <- t(apply(dropout_table, 1, function(row) {
  row[row == 0] <- min(row[row != 0], na.rm = TRUE)
  return(row)
}))
log2_dropout_table <- log2(dropout_table)
RT_values <- log2_dropout_table[!grepl("noRT", rownames(log2_dropout_table)), ]
noRT_values <- log2_dropout_table[grepl("noRT", rownames(log2_dropout_table)), ]
colnames(RT_values) <- paste(colnames(RT_values), "_RT", sep = "")
colnames(noRT_values) <- paste(colnames(noRT_values), "_noRT", sep = "")
dropout_map <- c(colMeans(RT_values), colMeans(noRT_values))
saveRDS(dropout_map, file.path(DATA_DIR, "gbm43_sample_averaged_log2FC_dropout_map.rds"))

# Calculate a KD map purely based on the sgRNACond counts over all the samples.
sgRNACond_counts <- table(metadata_table$sgRNACond)
noRT_perturbs <- names(sgRNACond_counts)[grepl("noRT", names(sgRNACond_counts))]
RT_perturbs <- names(sgRNACond_counts)[!grepl("noRT", names(sgRNACond_counts))]
noRT_ratios <- log2(sgRNACond_counts[noRT_perturbs] / sgRNACond_counts["non-targeting_noRT"])
RT_ratios <- log2(sgRNACond_counts[RT_perturbs] / sgRNACond_counts["non-targeting_RT"])
global_dropout_map <- c(noRT_ratios, RT_ratios)
saveRDS(global_dropout_map, file.path(DATA_DIR, "gbm43_global_average_log2FC_dropout_map.rds"))

# =========================================================
# Unified plotting
# =========================================================
# Plot the basic GSEA plot, with overlaid knockdown and dropout values
bubble_plot_output <- utils.gsea$make_bubble_plot(
  fgsea_results_list,
  color_overlay = FALSE,
  normalization_scheme = "noRTNormalized",
  color_label = "NES",
)
ggsave(file.path(OUTPUT_DIR, "basic_gsea_plot_with_NES.pdf"), bubble_plot_output$plot,
  height = 13,
  width = 6.7 + length(fgsea_results_list) / 27 * (11 - 6.7) / 0.8,
  device = "pdf"
)

# Overlay LFCs if the significant deseq data contains 20% of the gene set
median_matrix <- generate_blank_matrix()
for (i in seq_along(msigdbr_gene_set_lists)) {
  gene_set <- msigdbr_gene_set_lists[[i]]
  for (j in seq_along(significant_deseq_data)) {
    perturb_df <- significant_deseq_data[[j]]
    perturb_df <- perturb_df[perturb_df$feature %in% gene_set, ]
    median_lfc <- median(perturb_df$log_fc, na.rm = TRUE)
    if (nrow(perturb_df) > 0.2 * length(gene_set)) {
      median_log_fc <- median(perturb_df$log_fc, na.rm = TRUE)
      median_matrix[i, j] <- median_log_fc
    } else {
      median_matrix[i, j] <- 0
    }
  }
}
rownames(median_matrix) <- gsub("HALLMARK_", "", rownames(median_matrix))
colnames(median_matrix) <- gsub("_non-targeting_noRT", "", colnames(median_matrix))
bubble_plot_output <- utils.gsea$make_bubble_plot(
  fgsea_output_list,
  color_overlay = TRUE,
  color_lookup_matrix = median_matrix,
  normalization_scheme = "noRTNormalized",
  color_label = "LFC",
  color_limits = c(-0.75, 0.75)
)
ggsave(file.path(OUTPUT_DIR, "gsesa_plot_with_LFC_overlay.pdf"), bubble_plot_output$plot,
  height = 13,
  width = 6.7 + length(bubble_plot_output$col_order) / 27 * (11 - 6.7) / 0.8,
  device = "pdf"
)

# Overlay LFCs from leading edge genes only
# Pull leading edge genes out of the fgsea output
leading_edge_gene_sets <- lapply(fgsea_results_list, function(df) {
  df <- df %>%
    select(pathway, leadingEdge) %>%
    column_to_rownames("pathway")
  return(df)
})

# Generate a lookup table that includes the leading edge genes
msigdbr_df <- msigdbr(species = "Homo sapiens", category = "H")
msigdbr_gene_set_lists <- split(msigdbr_df$gene_symbol, msigdbr_df$gs_name)

# Remove perturbations that have less than 100 significant measurements.
significant_deseq_data <- Filter(function(df) nrow(df) > 100, significant_deseq_data)

# Generate a blank matrix to store values in
generate_blank_matrix <- function() {
  gene_set_names <- names(msigdbr_gene_set_lists)
  perturb_names <- names(deseq_data)
  result_matrix <- matrix(nrow = length(gene_set_names), ncol = length(perturb_names))
  rownames(result_matrix) <- gene_set_names
  colnames(result_matrix) <- perturb_names
  return(result_matrix)
}

# Average LFC over gene sets and plot heatmap.
leading_edge_mean_LFC_matrix <- generate_blank_matrix()
for (pathway in names(msigdbr_gene_set_lists)) {
  for (perturb in names(leading_edge_gene_sets)) {
    leading_edge_gene_set <- leading_edge_gene_sets[[perturb]][pathway, "leadingEdge"][[1]]
    perturb_df <- deseq_data[[perturb]]
    perturb_df <- perturb_df[perturb_df$feature %in% leading_edge_gene_set, ]
    mean_lfc <- median(perturb_df$log_fc, na.rm = TRUE)
    leading_edge_mean_LFC_matrix[pathway, perturb] <- mean_lfc
  }
}
rownames(leading_edge_mean_LFC_matrix) <- gsub("HALLMARK_", "", rownames(leading_edge_mean_LFC_matrix))
colnames(leading_edge_mean_LFC_matrix) <- gsub("_non-targeting_noRT", "", colnames(leading_edge_mean_LFC_matrix))
bubble_plot_output <- utils.gsea$make_bubble_plot(
  fgsea_results_list,
  color_overlay = TRUE,
  color_lookup_matrix = leading_edge_mean_LFC_matrix,
  normalization_scheme = "noRTNormalized",
  color_label = "Leading edge LFC, median",
  color_limits = c(-2, 2)
)
ggsave(file.path(OUTPUT_DIR, "gsea_plot_with_leadingEdge_LFC_overlay_median_non_sig.pdf"), bubble_plot_output$plot,
  height = 13,
  width = 6.7 + length(bubble_plot_output$col_order) / 27 * (11 - 6.7) / 0.8,
  device = "pdf"
)

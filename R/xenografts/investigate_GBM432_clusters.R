# Prepare the xenograft data after pipeline running for Figure 2 and supplements
# of the manuscript. Investigate clustering for the mouse microenvironment cells
# and tumor cells. Generate feature plots with marker genes for each of the cell
# groups we're calling. GBM432_5 is used to generate this.
#
# Create for manuscript preparation on 4/22/25. Last updated 5/5/25.

library(Seurat)
library(dplyr)
library(SCpubr)
library(ggplot2)
library(reticulate)
library(EnhancedVolcano)

pipelines_dir <- "output/xenografts/pipelines"
tumor_dir <- "output/xenografts/pipelines/gbm432_5/tumor_analysis"
menv_dir <- "output/xenografts/pipelines/gbm432_5/menv_analysis"
output_dir <- "output/xenografts/playground/gbm432_tumor_menv_clustering"

utils_deseq <- new.env()
source("R/utils/deseq.R", local = utils_deseq)

# =========================================================
# Investigate tumor clustering
# =========================================================
# The tumor clusters look clearest for louvain resolution 0.3.
# We will export markers to the D-SPIN package so we can run GPT annotation on
# them, and compare them to Enrichr.
tumor_data <- readRDS(file.path(tumor_dir, "postClustering_tumor_merged.rds"))

# Attempt to use the existing markers to calculate clusters
louvain_03_markers <- read.table(file.path(tumor_dir, "louvain_03_markers.txt"))

# What are the top markers? Right now they're grouped by p-value
top_markers <- louvain_03_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 30) %>%
  ungroup()

# Write these into a CSV format that can be taken up by the D-SPIN program
top_30_genes_only <- top_markers[c("gene", "cluster")]
clusters_list <- split(top_30_genes_only$gene, top_30_genes_only$cluster)
dspin_df <- as.data.frame(do.call(cbind, clusters_list))
write.table(dspin_df,
  file = file.path(
    output_dir, "gbm432_5_tumor",
    "tumor_louvain_03_markers_dspin_format.csv"
  ),
  sep = ",", quote = FALSE, row.names = FALSE
)

# Results from Enrichr and from GPT annotations:
# 0: Non-proliferative, metabolic stress
# 1: NPC-like, RNA-processing state
# 2: DLK1, CTHRC-drive, stem-like?
# 3: Proliferative and transcriptionally active
# 4: Non-proliferative, synthesis active state
# 5: MES-state with high extracellular matrix remodeling and adhesion
# # See how much is AC and how much is mitotic

# These aren't quite good enough yet.
for (cluster in c("0", "1", "2", "3", "4", "5")) {
  feature_plot <- FeaturePlot(tumor_data, dspin_df[[cluster]])
  ggsave(
    file.path(
      output_dir, "gbm432_5_tumor",
      sprintf("tumor_feature_plot_%s.png", cluster)
    ),
    plot = feature_plot, width = 20, height = 40
  )
}

# Note that the markers we calculated are not sufficient to fully separate
# out cell populations. We know this because we tried:
# (see louvain_03_markers_failed_separation directory). Thus, we will
# recalculate the markers with certain restrictions.
Idents(tumor_data) <- "leiden_0.3" # Note: misnomer. It's actually Louvain
louvain_03_markers_025pct025lfc <- FindAllMarkers(
  tumor_data,
  logfc.threshold = 0.25,
  min.diff.pct = 0.25,
)
top_markers <- louvain_03_markers_025pct025lfc %>%
  group_by(cluster) %>%
  slice_head(n = 30) %>%
  ungroup()

# What happens if we plot these ones only on a DimPlot?
dimplot <- DoHeatmap(
  tumor_data,
  features = top_markers$gene
) + NoLegend()
ggsave(
  file.path(
    output_dir, "gbm432_5_tumor",
    "louvain_03_markers_025pct025lfc_heatmap.png"
  ),
  height = 20, width = 20
)

# Save the markers
write.table(top_markers,
  file = file.path(
    output_dir,
    "gbm432_5_tumor",
    "louvain_03_025pct_025lfc_markers.txt"
  ),
  row.names = FALSE
)
# Write these into a CSV format that can be taken up by the D-SPIN program
top_30_genes_only <- top_markers[c("gene", "cluster")]
clusters_list <- split(top_30_genes_only$gene, top_30_genes_only$cluster)
dspin_df <- as.data.frame(do.call(cbind, clusters_list))
write.table(dspin_df,
  file = file.path(
    output_dir, "gbm432_5_tumor",
    "tumor_louvain_03_025pct_025lfc_markers_dspin_format.csv"
  ),
  sep = ",", quote = FALSE, row.names = FALSE
)

# We can see that interestingly enough, we can identify cluster 4 for example
# via the ABSENCE of some of the other ubiquitous markers. Let's try plotting
# the revised markers
for (cluster in c("0", "1", "2", "3", "4", "5")) {
  feature_plot <- FeaturePlot(tumor_data, dspin_df[[cluster]])
  ggsave(
    file.path(
      output_dir, "gbm432_5_tumor",
      sprintf("tumor_feature_plot_%s.png", cluster)
    ),
    plot = feature_plot, width = 20, height = 40
  )
}

# One interesting thing we can try to do is plot the expression via the GBM6
# markers. After all, we want to call consistent cell states.
gbm6_markers <- read.table(
  "output/xenografts/playground/gbm6_tumor_menv_clustering/gbm6_tumor/tumor_louvain_0.2_markers_dspin_format.csv",
  header = TRUE
)

# Combine these markers together and plot a heatmap
all_markers <- unlist(as.list(gbm6_markers))
gbm6_marker_heatmap_of_gbm432_expression <- DoHeatmap(
  tumor_data,
  features = all_markers
)
ggsave(
  sprintf(
    "%s/gbm432_5_tumor/gbm6_marker_heatmap_of_gbm432_expression.png",
    output_dir
  ),
  height = 15, width = 10
)

# Based on these markers, there are some that are shared
# 5: MES-like, based on COL1A2, VCAM1, EPHA3
# 1: NPC-like, based on DDX17 and TCF4
# 2: AC-like/DLK1/CTHRC, based on PTN
# 0, 2, and 4 share the marker S100A4, which is AC-like

# So where do we go from here? Final cluster assignments:
# 0: Proliferative, cell stress
# 1: NPC-like
# 2: AC-like, DLK1/CTHRC1 program with PTN and S100A4
# 3: Proliferative
# 4: Stalled
# 5: MES-like

marker_genes <- c(
  "CCL2", # 0: stress
  "DNAJC15", # 0: stress
  "KLF6", # 0
  "MEG8", # 1
  "PLAG1", # 1
  "DDX17", # 1
  "TCF4", # 1
  "DLK1", # 2
  "CTHRC1", # 2
  "GCHFR", # 2
  "S100A4", # 2
  "PTN", # 2,
  "STMN1", # 0, 3, 4: proliferative marker
  "CKS2", # 0, 3, 4: proliferative marker
  "UBE2T", # 0, 3, 4: proliferative marker
  "CCNB2", # 0, 3, 4: proliferative marker
  "DDX5", # 4: Ribosome and splicing deficient
  "TUBB", # 4: Cytoskeleton deficient
  "EEF1G", # 4: Protein synthesis stalling
  "COL1A2", # 5: MES-like
  "PLOD2", # 5: MES-like
  "VCAM1", # 5: MES-like
  "FN1" # 5 : MES-like
)

for (gene in marker_genes) {
  print(gene)
  feature_plot <- SCpubr::do_FeaturePlot(tumor_data, gene,
    order = TRUE, plot.title = gene,
    plot.subtitle = "GBM432",
    legend.length = 5, legend.width = 0.5
  )
  ggsave(
    file.path(
      output_dir, "gbm432_5_tumor",
      "gbm432_marker_plots",
      sprintf("tumor_feature_plot_%s.png", gene)
    ),
    plot = feature_plot, width = 5, height = 5.5, units = "in"
  )
}

# Plot a new UMAP plot that labels the cells by their cluster identity
cluster_identity_mapping <- list(
  "0" = "GBM-stress",
  "1" = "GBM-NPC-like",
  "2" = "GBM-AC-like",
  "3" = "GBM-Mitotic",
  "4" = "GBM-Inactive",
  "5" = "GBM-MES-like"
)
map_cluster_to_cell_type <- function(cluster_number) {
  return(cluster_identity_mapping[[cluster_number]])
}
tumor_data$annotated_cell_state <- sapply(
  tumor_data$leiden_0.3,
  map_cluster_to_cell_type
)
dimplot <- DimPlot(
  tumor_data,
  reduction = "umap", group.by = "annotated_cell_state",
  label = TRUE, label.box = TRUE, label.size = 1.5
)
ggsave(file.path(output_dir, "gbm432_5_tumor", "tumor_cell_state_dimplot.pdf"),
  dimplot,
  height = 5, width = 7, device = "pdf"
)
saveRDS(
  tumor_data,
  file.path(output_dir, "gbm432_5_tumor", "finalized_cluster_object.rds")
)

# Cool, this is done. Let's generate the UMAP coordinates and table with cell
# counts by treatment condition.
tumor_data_umap <- tumor_data[["umap"]]@cell.embeddings
tumor_data_with_embeddings <- merge(tumor_data_umap, tumor_data@meta.data,
  by = 0
)
colnames(tumor_data_with_embeddings)[1] <- "Cell barcodes"
tumor_data_with_embeddings <- tumor_data_with_embeddings[c(
  "Cell barcodes",
  "umap_1",
  "umap_2",
  "annotated_cell_state"
)]
write.table(tumor_data_with_embeddings,
  file.path(
    output_dir, "gbm432_5_tumor",
    "tumor_embeddings_and_cell_types_only.txt"
  ),
  sep = "\t", quote = FALSE, row.names = FALSE
)

frequency_table <- as.data.frame(table(
  tumor_data$annotated_cell_state,
  tumor_data$treatment_cond
))
colnames(frequency_table) <- c("Cell type", "Treatment", "Count")
write.table(frequency_table,
  file = file.path(
    output_dir,
    "gbm432_5_tumor",
    "frequency_table_across_tumor_states.txt"
  ),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# =========================================================
# Investigate microenvironment clustering
# =========================================================
# We want to achieve sufficient resolution on the total microenvironment to
# characterize first: immune cells and all of the other cells; then characterize
# the immune cells.
menv_data <- readRDS(file.path(
  menv_dir,
  "postSCMRMA_microenvironment_merged.rds"
))

# Choose louvain 0.2 because it offers qualitative clean marker separation on
# the heatmap
louvain_02_markers <- read.table(file.path(menv_dir, "leiden_0.2_markers"))

# What are the top markers? Right now they're grouped by p-value
top_markers <- louvain_02_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 30) %>%
  ungroup()

# Write these into a CSV format that can be taken up by the D-SPIN program
top_30_genes_only <- top_markers[c("gene", "cluster")]
clusters_list <- split(top_30_genes_only$gene, top_30_genes_only$cluster)
dspin_df <- as.data.frame(do.call(cbind, clusters_list))
write.table(dspin_df,
  file = file.path(
    output_dir,
    "gbm432_5_menv",
    "menv_louvain_03_markers_dspin_format.csv"
  ),
  sep = ",", quote = FALSE, row.names = FALSE
)

# Let's assign identities to each of the clusters. There are 11 total:
cluster_to_cell_type_map <- list(
  "0" = "Neurons", # Nrgn, Nrg4, Snap25, Meg3
  "1" = "Neurons", # Nrgn, Nrg4, Snap25, Meg3
  "2" = "Immune cells", # Cd52, Tyrobp, Fcer1g
  "3" = "Immune cells", # P2ry12, C1qc, Csf1r
  "4" = "Endothelial Cells", # Cldn5, Flt1 here
  "5" = "Oligodendrocytes", # Plp1, Mobp, Mog - lots of myelination
  "6" = "Astrocytes", # Slc1a3, Gja1, Atp1a2 all here
  "7" = "Cerebellar Neurons", # Scg5, Chgb here
  "8" = "Oligodendrocyte Precursor Cell", # Pdgfra, Olig1
  "9" = "Erythrocytes", # Hbb-a1, Hbb-bs, Alas2
  "10" = "Ependymal Cells" # Ttr, Ecrg4
)
map_cluster_to_cell_type <- function(cluster_number) {
  return(cluster_to_cell_type_map[[cluster_number]])
}
menv_data$annotated_cell_type <- sapply(
  menv_data$leiden_0.2,
  map_cluster_to_cell_type
)
dimplot <- DimPlot(
  menv_data,
  reduction = "umap", group.by = "annotated_cell_type", label = TRUE,
  label.box = TRUE, label.size = 1
)
ggsave(
  file.path(
    output_dir, "gbm432_5_menv",
    "menv_cell_type_dimplot_immune_inset.pdf"
  ), dimplot,
  height = 5, width = 7, device = "pdf"
)

# Alright, let's generate marker plots for all of these

marker_genes <- c(
  # Neurons (clusters 0, 1)
  "Nrgn", "Ndrg4", "Snap25", "Meg3",
  # Immune cells (clusters 2, 3)
  "Csf1r", "Tyrobp", "Cd52", "Fcer1g", "C1qc", "P2ry12",
  # Endothelial Cells (cluster 4)
  "Cldn5", "Flt1",
  # Mature Oligodendrocytes (cluster 5)
  "Plp1", "Mobp", "Mog",
  # Astrocytes (cluster 6)
  "Slc1a3", "Gja1", "Atp1a2",
  # Secretory Neurons (cluster 7)
  "Scg5", "Chgb",
  # Oligodendrocyte Precursor Cell (cluster 8)
  "Pdgfra", "Olig1",
  # Erythrocytes (cluster 9)
  "Hba-a1", "Hbb-bs", "Alas2",
  # Ependymal Cells (cluster 10)
  "Ttr", "Ecrg4"
)
for (gene in marker_genes) {
  print(gene)
  feature_plot <- SCpubr::do_FeaturePlot(
    menv_data, gene,
    order = TRUE, plot.title = gene, plot.subtitle = "GBM432",
    legend.length = 5, legend.width = 0.5
  )
  ggsave(
    file.path(
      output_dir, "gbm432_5_menv",
      "louvain_0_2_marker_feature_plots",
      sprintf("marker_feature_%s.png", gene)
    ),
    plot = feature_plot, width = 5, height = 5.5, units = "in"
  )
}

# Subset to the immune data so we can do immune-specific downstream analysis
immune_data <- subset(menv_data, annotated_cell_type == "Immune cells")
immune_data <- SCTransform(immune_data)
immune_data <- RunPCA(immune_data)
immune_data <- FindNeighbors(immune_data)
resolutions <- c(0.1, 0.2, 0.3, 0.5, 1)
for (res in resolutions) {
  immune_data <- FindClusters(
    immune_data,
    algorithm = 1, # Louvain
    resolution = res,
    cluster.name = sprintf("louvain_%s", res)
  )
}
immune_data <- RunUMAP(immune_data, dims = 1:50)
saveRDS(immune_data, file = file.path(
  output_dir,
  "gbm432_5_immune",
  "immune_data.rds"
))

# Visualize the subset results
umap_dimplot <- DimPlot(
  immune_data,
  reduction = "umap",
  group.by = c(sprintf("louvain_%s", resolutions), "sample", "treatment_cond")
)
ggsave(
  file.path(
    output_dir, "gbm432_5_immune",
    "immune_umap_against_discrete_metadata.pdf"
  ),
  plot = umap_dimplot, device = "pdf", height = 15, width = 20
)

# Find markers for each of the clusters and identify them
for (res in resolutions) {
  Idents(immune_data) <- sprintf("louvain_%s", res)
  df <- FindAllMarkers(
    immune_data
  )
  write.table(df, file = file.path(
    output_dir, "gbm432_5_immune",
    sprintf("louvain_%s_markers.txt", res)
  ))
  df %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5) %>%
    slice_head(n = 20) %>%
    ungroup() -> top20
  marker_heatmap <- DoHeatmap(immune_data, features = top20$gene) + NoLegend()
  ggsave(
    file.path(
      output_dir, "gbm432_5_immune",
      sprintf("louvain_%s_marker_heatmap.pdf", res)
    ),
    plot = marker_heatmap, device = "pdf", height = 20, width = 20
  )
}

# Choose Louvain 0.1
louvain_0.1_markers <- read.table(file.path(
  output_dir,
  "gbm432_5_immune",
  "louvain_0.1_markers.txt"
))

# What are the top markers? Right now they're grouped by p-value
top_markers <- louvain_0.1_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 30) %>%
  ungroup()

# Write these into a CSV format that can be taken up by the D-SPIN program
top_30_genes_only <- top_markers[c("gene", "cluster")]
clusters_list <- split(top_30_genes_only$gene, top_30_genes_only$cluster)
dspin_df <- as.data.frame(do.call(cbind, clusters_list))
write.table(dspin_df,
  file = file.path(
    output_dir,
    "gbm432_5_immune",
    "immune_louvain_0.1_markers_dspin_format.csv"
  ),
  sep = ",", quote = FALSE, row.names = FALSE
)

# Assign labels
cluster_to_cell_type_map <- list(
  "0" = "Microglia", # C1qa, P2ry12
  "1" = "Macrophages", # Cd74, H2-Aa, Ccr2, Ifitm2, Lyz2
  "2" = "Neurons", # Nrgn, Snap25, Gria2
  "3" = "Macrophages" # Mrc1, Lyve1, F13a1
)
map_cluster_to_cell_type <- function(cluster_number) {
  return(cluster_to_cell_type_map[[cluster_number]])
}
immune_data$annotated_cell_type <- sapply(
  immune_data$louvain_0.1,
  map_cluster_to_cell_type
)
umap_dimplot <- DimPlot(
  immune_data,
  reduction = "umap",
  group.by = "annotated_cell_type"
)
ggsave(
  file.path(
    output_dir, "gbm432_5_immune",
    "immune_umap_against_annotated_cell_type.pdf"
  ),
  plot = umap_dimplot, device = "pdf", height = 5, width = 7
)

# Decide on a set of feature plots for the insets:
marker_genes <- c(
  "C1qa", # Microglia
  "P2ry12", # Microglia
  "Cd74", # Macrophages
  "H2-Aa", # Macrophages
  "Ccr2", # Macrophages
  "Ifitm2", # Macrophages
  "Lyz2", # Macrophages
  "Nrgn", # Microglia-like neurons
  "Snap25", # Microglia-like neurons
  "Gria2", # Microglia-like neurons
  "Mrc1", # Macrophages
  "Lyve1", # Macrophages
  "F13a1" # Macrophages
)

for (cluster in c("0", "1", "2", "3")) {
  cluster_genes <- dspin_df[[cluster]]
  feature_plot <- SCpubr::do_FeaturePlot(immune_data, cluster_genes,
    order = TRUE, individual.captions = cluster_genes
  )
  ggsave(file.path(
    output_dir, "gbm432_5_immune", "louvain_0.1_marker_feature_plots",
    sprintf("louvain_%s_markers.png", cluster)
  ), width = 20, height = 40)
}
for (gene in marker_genes) {
  feature_plot <- SCpubr::do_FeaturePlot(immune_data, gene,
    order = TRUE, plot.title = gene, plot.subtitle = "GBM432",
    legend.length = 5, legend.width = 0.5
  )
  ggsave(
    file.path(
      output_dir, "gbm432_5_immune",
      "louvain_0.1_marker_feature_plots",
      sprintf("marker_feature_%s.png", gene)
    ),
    plot = feature_plot, width = 5, height = 5.5, units = "in"
  )
}

# Run DESeq and generate volcano plots for macrophages and microglia
macrophages_microglias <- subset(
  immune_data,
  annotated_cell_type %in% c(
    "Microglia",
    "Macrophage",
    "Perivascular macrophages"
  )
)
macrophages_microglias_dir <- file.path(
  output_dir, "gbm432_5_immune",
  "macrophages_microglias"
)
saveRDS(macrophages_microglias, file.path(
  macrophages_microglias_dir,
  "macrophages_microglias.rds"
))
treatment_conds <- unique(macrophages_microglias$treatment_cond)
for (treatment in treatment_conds[treatment_conds != "untreated_noRT"]) {
  deseq_df <- utils_deseq$find_deseq_differential_genes(
    macrophages_microglias,
    treatment,
    "untreated_noRT",
    "treatment_cond",
    seed = 5220
  )
  write.table(deseq_df, file = file.path(
    macrophages_microglias_dir, sprintf("%s_untreated_noRT.csv", treatment)
  ), sep = "\t", row.names = FALSE, quote = FALSE)
  volcano_plot <- EnhancedVolcano(
    deseq_df,
    deseq_df$feature,
    "log_fc",
    "padj",
    pCutoff = 1e-5,
    labSize = 3,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    subtitle = "Padj < 1e-5, abs(LFC) > 1"
  )
  ggsave(file.path(macrophages_microglias_dir, sprintf("%s_volcano.pdf", treatment)),
    volcano_plot,
    height = 10, width = 10, device = "pdf"
  )
}

# Ok with the immune cells in mind, let's go back to the global UMAP and
# relabel the immune cells with their appropriate labels
menv_data$annotated_cell_type_with_immune_subcategories <-
  menv_data$annotated_cell_type
menv_data@meta.data[
  rownames(immune_data@meta.data),
  "annotated_cell_type_with_immune_subcategories"
] <-
  immune_data@meta.data[, "annotated_cell_type"]
umap_dimplot <- DimPlot(
  menv_data,
  reduction = "umap",
  group.by = "annotated_cell_type_with_immune_subcategories",
  label = TRUE,
  label.size = 2
)
ggsave(
  file.path(
    output_dir, "gbm432_5_menv",
    "menv_cell_type_dimplot_immune_subset.pdf"
  ),
  plot = umap_dimplot, device = "pdf", height = 5, width = 7
)

# Now that we've decided we're going to use an inset,
# let's get the cell frequencies.
frequency_table <- as.data.frame(table(
  menv_data$annotated_cell_type,
  menv_data$treatment_cond
))
colnames(frequency_table) <- c("Cell type", "Treatment", "Count")
write.table(frequency_table,
  file = file.path(
    output_dir, "gbm432_5_menv",
    "frequency_table_across_annotations_immune_inset.txt"
  ),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Get the cell frequencies with the immune subset split
frequency_table <- as.data.frame(table(menv_data$annotated_cell_type_with_immune_subcategories, menv_data$treatment_cond))
colnames(frequency_table) <- c("Cell type", "Treatment", "Count")
write.table(frequency_table,
  file = file.path(
    output_dir, "gbm432_5_menv",
    "frequency_table_across_annotations_including_immune_subcategories.txt"
  ),
  sep = "\t", row.names = FALSE, quote = FALSE
)
dimplot <- DimPlot(
  menv_data,
  reduction = "umap", group.by = "leiden_0.3", label = TRUE
)
ggsave(
  file.path(output_dir, "dimplot.png"),
  height = 10, width = 12
)

# Get the immune cell frequencies
frequency_table <- as.data.frame(table(immune_data$annotated_cell_type, immune_data$treatment_cond))
colnames(frequency_table) <- c("Cell type", "Treatment", "Count")
write.table(frequency_table,
  file = file.path(
    output_dir, "gbm432_5_immune",
    "frequency_table_across_cell_types.txt"
  ),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Save the metadata from the immune data and from the menv data, with UMAP embeddings
menv_data_umap <- menv_data[["umap"]]@cell.embeddings
menv_data_with_embeddings <- merge(menv_data_umap, menv_data@meta.data, by = 0)
colnames(menv_data_with_embeddings)[1] <- "Cell barcodes"
menv_data_with_embeddings <- menv_data_with_embeddings[c("Cell barcodes", "umap_1", "umap_2", "annotated_cell_type")]
write.table(menv_data_with_embeddings, file.path(output_dir, "gbm432_5_menv", "menv_embeddings_and_cell_types_only.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

immune_data_umap <- immune_data[["umap"]]@cell.embeddings
immune_data_with_embeddings <- merge(immune_data_umap, immune_data@meta.data, by = 0)
colnames(immune_data_with_embeddings)[1] <- "Cell barcodes"
immune_data_with_embeddings <- immune_data_with_embeddings[c("Cell barcodes", "umap_1", "umap_2", "annotated_cell_type")]
write.table(immune_data_with_embeddings, file.path(output_dir, "gbm432_5_menv", "immune_data_metadata.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)


# # =========================================================
# # Combined UMAP plotting
# # =========================================================
# # For both GBM6 and GBM43, do the following:
# # 1. Load the tumor and menv post-QC, pre-clustering objects
# # 2. Merge them together and run the clustering pipeline
# # 3. Cluster at various resolutions, then try to determine a set of cells for
# #    immunogenicity
#
# resolutions <- c(0.1, 0.2, 0.3, 0.5, 1)
#
# tumor <- readRDS(file.path(pipelines_dir, "gbm6_5", "tumor_analysis", "postQC_tumor_merged.rds"))
# menv <- readRDS(file.path(pipelines_dir, "gbm6_5", "menv_analysis", "postQC_menv_merged.rds"))
# combined <- merge(tumor, menv)
# combined <- JoinLayers(combined)
# combined$species <- ifelse(combined$percent_human > 90, "human", "mouse")
#
# # Note that because of the size of globals needed, we need to set the maximum
# # allowed global size. Specify a size of around 2 GB
# options(future.globals.maxSize = 3 * 1024^3)
#
# # Run the SCTransform based pipeline
# combined <- SCTransform(combined)
# combined <- RunPCA(combined)
# combined <- FindNeighbors(combined)
# for (res in resolutions) {
#   combined <- FindClusters(
#     combined,
#     algorithm = 1, # Use Louvain for now
#     resolution = res,
#     cluster.name = sprintf("louvain_%s", res)
#   )
# }
# combined <- RunUMAP(combined, dims = 1:50)
# saveRDS(combined, file.path(output_dir, "combined.rds"))
#
# # Now that we've found the clusters...find markers for each?
# dim_loadings <- VizDimLoadings(data, dims = 1:4, reduction = "pca")
# ggsave(file.path(output_dir, "pca_dim_loadings.pdf"),
#   plot = dim_loadings,
#   device = "pdf", height = 10, width = 10
# )
#
# # Verify that the PCs do not overly cluster by sample
# pc_dimplot <- DimPlot(combined, reduction = "pca", group.by = c(
#   "sample", "treatment_cond", "species"
# ))
# ggsave(file.path(output_dir, "pca_dim_plot_treatment_samples_species.pdf"),
#   plot = pc_dimplot,
#   device = "pdf", height = 10, width = 20
# )
#
# # Plot various metadata against UMAP: sample, cluster resolutions, phase,
# # ribosomal gene percentage, mitochondrial gene percentage
# umap_dimplot <- DimPlot(
#   combined,
#   reduction = "umap",
#   group.by = c(sprintf("louvain_%s", resolutions), "sample", "treatment_cond", "species")
# )
# ggsave(file.path(output_dir, "umap_against_discrete_metadata.pdf"),
#   plot = umap_dimplot, device = "pdf", height = 15, width = 20
# )
#
# # Find proportion of treatment_cond in clusters and clusters in treatment_cond
# for (res in RESOLUTIONS) {
#   cluster_res <- sprintf("leiden_%s", res)
#   proportions <- data@meta.data %>%
#     group_by(treatment_cond, !!sym(cluster_res)) %>%
#     dplyr::summarize(intersect_counts = n()) %>%
#     group_by(treatment_cond) %>%
#     mutate(proportions = intersect_counts / sum(intersect_counts))
#   proportion_plot <- ggplot(proportions, aes(x = treatment_cond, y = proportions, fill = !!sym(cluster_res))) +
#     geom_bar(stat = "identity", position = "fill") +
#     labs(
#       title = "Proportion of Leiden Clusters by Treatment Condition",
#       x = "Treatment Condition",
#       y = "Proportion",
#       fill = "Leiden Cluster"
#     ) +
#     theme_minimal() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1))
#   ggsave(file.path(TUMOR_output_dir, sprintf("%s_proportion_plot.pdf", cluster_res)),
#     plot = proportion_plot, height = 10, width = 30
#   )
# }
#
#
# # Find markers for the clusters
# for (res in RESOLUTIONS) {
#   Idents(data) <- sprintf("leiden_%s", res)
#   df <- FindAllMarkers(
#     data
#   )
#   write.table(df, file = file.path(TUMOR_output_dir, sprintf("leiden_%s_markers", res)))
#
#   df %>%
#     group_by(cluster) %>%
#     dplyr::filter(avg_log2FC > 0.5) %>%
#     slice_head(n = 20) %>%
#     ungroup() -> top20
#   marker_heatmap <- DoHeatmap(data, features = top20$gene) + NoLegend()
#   ggsave(file.path(TUMOR_output_dir, sprintf("leiden_%s_marker_heatmap.pdf", res)),
#     plot = marker_heatmap, device = "pdf", height = 20, width = 20
#   )
# }
#
# saveRDS(data, file.path(TUMOR_output_dir, "postClustering_tumor_merged.rds"))

# Prepare the xenograft data after pipeline running for Figure 2 and supplements
# of the manuscript. Investigate clustering for the mouse microenvironment cells
# and tumor cells. Generate feature plots with marker genes for each of the cell
# groups we're calling.
#
# Create for manuscript preparation (Figure 2) on 3/19/25. Last updated 3/25/25.

library(Seurat)
library(dplyr)
library(SCpubr)
library(ggplot2)
library(reticulate)
library(EnhancedVolcano)

PIPELINES_DIR <- "output/xenografts/pipelines"
TUMOR_DIR <- "output/xenografts/pipelines/gbm6_5/tumor_analysis"
MENV_DIR <- "output/xenografts/pipelines/gbm6_5/menv_analysis"
OUTPUT_DIR <- "output/xenografts/playground/gbm6_tumor_menv_clustering"

utils.deseq <- new.env()
source("R/utils/deseq.R", local = utils.deseq)

# =========================================================
# Investigate tumor clustering
# =========================================================
# Given that the tumor clusters look extremely clear for louvain resolution 0.2,
# let's start there. What are the different cell groups that we can call? Export
# markers to the D-SPIN package so we can run a GPT annotation on them. Also
# use Enrichr.
tumor_data <- readRDS(file.path(TUMOR_DIR, "postClustering_tumor_merged.rds"))
louvain_0.2_markers <- read.table(file.path(TUMOR_DIR, "louvain_0.2_markers.txt"))

# What are the top markers? Right now they're grouped by p-value
top_markers <- louvain_0.2_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 60) %>%
  ungroup()

# Write these into a CSV format that can be taken up by the D-SPIN program
top_30_genes_only <- top_markers[c("gene", "cluster")]
clusters_list <- split(top_30_genes_only$gene, top_30_genes_only$cluster)
dspin_df <- as.data.frame(do.call(cbind, clusters_list))
write.table(dspin_df, file = file.path(OUTPUT_DIR, "gbm6_tumor", "tumor_leiden_0.2_markers_dspin_format.csv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Interesting. Results from Enrichr and from GPT annotations:
# 0: Neural development and synaptic function (0.78) AC-like
# 1: Mitotic spindle assembly and mitosis (0.95) Mitotic
# 2: Developmental gene containing composite program (Neurogenesis) MES-like
# 3: Hypoxia response and metabolic adaptation (0.85) Hypoxic
# 4: Neural development and synaptic function (0.78) NPC-like
# 5: Immune response and interferon signaling (0.92) (Interferon signaling)

# I'm ok with these annotations. Let's generate a ton of feature plots so
# we can evaluate which marker genes will be best for showing the clusters.
for (cluster in c("0", "1", "2", "3", "4", "5")) {
  feature_plot <- FeaturePlot(tumor_data, dspin_df[[cluster]])
  ggsave(file.path(OUTPUT_DIR, "gbm6_tumor", sprintf("tumor_feature_plot_%s.pdf", cluster)), plot = feature_plot, device = "pdf", width = 20, height = 40)
}

# Now that we've identified a good marker for each, let's plot individual scpubr
# feature plots with optimized dimensions and features.
marker_genes <- c(
  "HEPN1",
  "CTNND2",
  "TOP2A",
  "CENPF",
  "MKI67",
  "IGFBP7",
  "ERBB4",
  "HSPB6",
  "DDIT4",
  "ENO2",
  "NDUFA4L2",
  "TCF4",
  "ADGRL2",
  "TENM3",
  "CXCL10",
  "CXCL11",
  "ISG20",
  "IFITM1"
)
for (gene in marker_genes) {
  feature_plot <- SCpubr::do_FeaturePlot(tumor_data, gene,
    order = TRUE, plot.title = gene, plot.subtitle = "GBM6",
    legend.length = 5, legend.width = 0.5
  )
  ggsave(file.path(OUTPUT_DIR, "gbm6_tumor", "gbm6_tumor_marker_plots", sprintf("tumor_feature_plot_%s.png", gene)),
    plot = feature_plot, width = 5, height = 5.5, units = "in"
  )
}

# =========================================================
# Investigate microenvironment clustering
# =========================================================
# We want to achieve sufficient resolution on the total microenvironment to
# characterize first: immune cells and all of the other cells; then characterize
# the immune cells.
menv_data <- readRDS(file.path(MENV_DIR, "postSCMRMA_microenvironment_merged.rds"))

# Choose louvain 0.3 because it has sufficient resolution to identify pericytes
louvain_0.3_markers <- read.table(file.path(MENV_DIR, "leiden_0.3_markers"))

# What are the top markers? Right now they're grouped by p-value
top_markers <- louvain_0.3_markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%
  slice_head(n = 30) %>%
  ungroup()

# Write these into a CSV format that can be taken up by the D-SPIN program
top_30_genes_only <- top_markers[c("gene", "cluster")]
clusters_list <- split(top_30_genes_only$gene, top_30_genes_only$cluster)
dspin_df <- as.data.frame(do.call(cbind, clusters_list))
write.table(dspin_df, file = file.path(OUTPUT_DIR, "gbm6_menv", "menv_louvain_0.3_markers_dspin_format.csv"), sep = ",", quote = FALSE, row.names = FALSE)

# Let's assign identities to each of the clusters. There are 14 total:
cluster_to_cell_type_map <- list(
  "0" = "Excitatory neurons",
  "1" = "Neurons",
  "2" = "Cerebellar neurons",
  "3" = "Dopaminergic neurons",
  "4" = "Interneurons",
  "5" = "Immune cells",
  "6" = "Astrocytes",
  "7" = "Immune cells",
  "8" = "Endothelial cells",
  "9" = "Oligodendrocytes",
  "10" = "Erythrocytes",
  "11" = "Oligodendrocytes",
  "12" = "Oligodendrocytes",
  "13" = "Pericytes"
)
map_cluster_to_cell_type <- function(cluster_number) {
  return(cluster_to_cell_type_map[[cluster_number]])
}
menv_data$annotated_cell_type <- sapply(menv_data$leiden_0.3, map_cluster_to_cell_type)
dimplot <- DimPlot(
  menv_data,
  reduction = "umap", group.by = "annotated_cell_type", label = TRUE, label.box = TRUE, label.size = 2
)
ggsave(file.path(OUTPUT_DIR, "gbm6_menv", "menv_cell_type_dimplot_immune_inset.pdf"), dimplot, height = 5, width = 7, device = "pdf")

# I'm ok with these annotations. Let's generate a ton of feature plots so
# we can evaluate which marker genes will be best for showing the clusters.
clusters <- c(
  "0", "1", "2", "3", "4", "5", "6", "7", "8", "9",
  "10", "11", "12", "13"
)
for (cluster in clusters) {
  feature_plot <- FeaturePlot(menv_data, dspin_df[[cluster]])
  ggsave(file.path(OUTPUT_DIR, "gbm6_menv", "louvain_0.3_marker_feature_plots", sprintf("menv_feature_plot_%s.png", cluster)), plot = feature_plot, width = 20, height = 40)
}

# Now that we've identified a good marker for each, let's plot individual scpubr
# feature plots with optimized dimensions and features.
marker_genes <- c(
  "Pdgfrb", # Pericytes, cluster 13
  "Notch3", # Pericytes, cluster 13
  "Hba-a1", # Erythrocytes, cluster 10
  "Hbb-bs", # Erythrocytes, cluster 10
  "Alas2", # Erythrocytes, cluster 10
  "Olig1", # Oligodendrocytes, clusters 9 11 12
  "Olig2", # Oligodendrocytes, clusters 9 11 12
  "Sox10", # Oligodendrocytes, clusters 9 11 12
  "Cldn5", # Epithelial cells, cluster 8
  "Flt1", # Epithelial cells, cluster 8
  "Slc1a3", # Astrocytes, cluster 6
  "Atp1a2", # Astrocytes, cluster 6
  "Gja1", # Astrocytes, cluster 6
  "Csf1r", # Immune cells, clusters 5 and 7
  "Tyrobp", # Immune cells, clusters 5 and 7
  "Cd52", # Immune cells, clusters 5 and 7
  "Fcer1g", # Immune cells, clusters 5 and 7
  "C1qc", # Immune cells, clusters 5 and 7
  "Rarb", # Neurons, cluster 3 - dopaminergic
  "Rgs9", # Neurons, cluster 3 - dopaminergic
  "Erbb4", # (Inhibitory) Interneurons, cluster 4
  "Kcnc2", # Interneurons, cluster 4
  "Gad1", # Interneurons, cluster 4
  "Gabra6", # Cerebellar neurons, cluster 2 (one half)
  "Cbln1", # Cerebellar neurons, cluster 2 (one half)
  "Scg5", # Secretory neurons, cluster 2 (other half)
  "Chgb", # Secretory neurons, cluster 2 (other half)
  "Pcsk1n", # Secretory neurons, cluster 2 (other half)
  "Dpp10", # Excitatory projection neurons, cluster 0
  "Nrg3", # Excitatory projection neurons, cluster 0
  "Nav3", # Excitatory projection neurons, cluster 0
  "Ptprd", # Excitatory projection neurons, cluster 0
  "Cacna2d1", # Neurons, cluster 1
  "Kcnt2", # Neurons, cluster 1
  "Cdh12" # Neurons, cluster 1
)

for (gene in marker_genes) {
  feature_plot <- SCpubr::do_FeaturePlot(
    menv_data, gene,
    order = TRUE, plot.title = gene, plot.subtitle = "GBM6",
    legend.length = 5, legend.width = 0.5
  )
  ggsave(file.path(OUTPUT_DIR, "gbm6_menv", "louvain_0.3_marker_feature_plots", sprintf("marker_feature_%s.png", gene)),
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
saveRDS(immune_data, file = file.path(OUTPUT_DIR, "immune_data.rds"))

# Visualize the subset results
umap_dimplot <- DimPlot(
  immune_data,
  reduction = "umap",
  group.by = c(sprintf("louvain_%s", resolutions), "sample", "treatment_cond")
)
ggsave(file.path(OUTPUT_DIR, "immune_umap_against_discrete_metadata.pdf"),
  plot = umap_dimplot, device = "pdf", height = 15, width = 20
)

# Find markers for each of the clusters and identify them
for (res in resolutions) {
  Idents(immune_data) <- sprintf("louvain_%s", res)
  df <- FindAllMarkers(
    immune_data
  )
  write.table(df, file = file.path(OUTPUT_DIR, "gbm6_immune", sprintf("louvain_%s_markers.txt", res)))
  df %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5) %>%
    slice_head(n = 20) %>%
    ungroup() -> top20
  marker_heatmap <- DoHeatmap(immune_data, features = top20$gene) + NoLegend()
  ggsave(file.path(OUTPUT_DIR, "gbm6_immune", sprintf("louvain_%s_marker_heatmap.pdf", res)),
    plot = marker_heatmap, device = "pdf", height = 20, width = 20
  )
}

# Choose Louvain 0.1
louvain_0.1_markers <- read.table(file.path(OUTPUT_DIR, "gbm6_immune", "louvain_0.1_markers.txt"))

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
write.table(dspin_df, file = file.path(OUTPUT_DIR, "gbm6_immune", "immune_louvain_0.1_markers_dspin_format.csv"), sep = ",", quote = FALSE, row.names = FALSE)

# Assign labels
map_cluster_to_cell_type <- function(cluster_number) {
  if (cluster_number == 0) {
    return("Microglia")
  } else if (cluster_number == 1) {
    return("Microglia-like neurons")
  } else if (cluster_number == 2) {
    return("Macrophage")
  } else if (cluster_number == 3) {
    return("Macrophage")
  } else {
    return("NK cell")
  }
}
immune_data$annotated_cell_type <- sapply(immune_data$louvain_0.1, map_cluster_to_cell_type)
umap_dimplot <- DimPlot(
  immune_data,
  reduction = "umap",
  group.by = "annotated_cell_type"
)
ggsave(file.path(OUTPUT_DIR, "gbm6_immune", "immune_umap_against_annotated_cell_type.pdf"),
  plot = umap_dimplot, device = "pdf", height = 5, width = 7
)

# Decide on a set of feature plots for the insets:
marker_genes <- c(
  "P2ry12", # Microglia
  "C1qa", # Activation
  "Nrxn1", # That neuronal cluster
  "Anks1b", # Neuronal cluster
  "Cxcl14", # Macrophage cluster 2
  "Spp1", # Macrophage cluster 2
  "Msr1", # Macrophage cluster 2
  "Cd74", # Macrophage cluster 3
  "H2-Aa", # Macrophage cluster 3
  "H2-Eb1", # Macrophage cluster 3
  "Ccr2", # Macrophage cluster 3
  "Nkg7", # NK cluster 4
  "Gzma" # NK cluster 4
)
for (cluster in c("0", "1", "2", "3", "4")) {
  cluster_genes <- dspin_df[[cluster]]
  feature_plot <- SCpubr::do_FeaturePlot(immune_data, cluster_genes,
    order = TRUE, individual.captions = cluster_genes
  )
  ggsave(file.path(
    OUTPUT_DIR, "gbm6_immune", "louvain_0.1_marker_feature_plots",
    sprintf("louvain_%s_markers.png", cluster)
  ), width = 20, height = 40)
}
for (gene in marker_genes) {
  feature_plot <- SCpubr::do_FeaturePlot(immune_data, gene,
    order = TRUE, plot.title = gene, plot.subtitle = "GBM6",
    legend.length = 5, legend.width = 0.5
  )
  ggsave(file.path(OUTPUT_DIR, "gbm6_immune", "louvain_0.1_marker_feature_plots", sprintf("marker_feature_%s.png", gene)),
    plot = feature_plot, width = 5, height = 5.5, units = "in"
  )
}

# Run DESeq and generate volcano plots for macrophages, microglia, and NKs
macrophages_microglia_nks <- subset(immune_data, annotated_cell_type %in% c("Microglia", "Macrophage", "NK cell"))
macrophage_microglia_nk_dir <- file.path(OUTPUT_DIR, "gbm6_immune", "macrophage_microglia_nk")
saveRDS(macrophages_microglia_nks, file.path(macrophage_microglia_nk_dir, "macrophages_microglia_nks.rds"))
treatment_conds <- unique(macrophages_microglia_nks$treatment_cond)
for (treatment in treatment_conds[treatment_conds != "untreated_noRT"]) {
  deseq_df <- utils.deseq$find_deseq_differential_genes(
    macrophages_microglia_nks,
    treatment,
    "untreated_noRT",
    "treatment_cond",
    seed = 5220
  )
  write.table(deseq_df, file = file.path(
    macrophage_microglia_nk_dir, sprintf("%s_untreated_noRT.csv", treatment)
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
  ggsave(file.path(macrophage_microglia_nk_dir, sprintf("%s_volcano.pdf", treatment)),
    volcano_plot,
    height = 10, width = 10, device = "pdf"
  )
}

# Run DESeq and generate volcano plots for macrophages and microglia only
macrophages_microglia <- subset(immune_data, annotated_cell_type %in% c("Microglia", "Macrophage"))
macrophage_microglia_dir <- file.path(OUTPUT_DIR, "gbm6_immune", "macrophage_microglia")
saveRDS(macrophages_microglia, file.path(macrophage_microglia_dir, "macrophages_microglia.rds"))
treatment_conds <- unique(macrophages_microglia$treatment_cond)
for (treatment in treatment_conds[treatment_conds != "untreated_noRT"]) {
  deseq_df <- utils.deseq$find_deseq_differential_genes(
    macrophages_microglia,
    treatment,
    "untreated_noRT",
    "treatment_cond",
    seed = 5220
  )
  write.table(deseq_df, file = file.path(
    macrophage_microglia_dir, sprintf("%s_untreated_noRT.csv", treatment)
  ))
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
  ggsave(file.path(macrophage_microglia_dir, sprintf("%s_volcano.pdf", treatment)),
    volcano_plot,
    height = 10, width = 10, device = "pdf"
  )
}

# Run DESeq and generate volcano plots for NKs only
nk <- subset(immune_data, annotated_cell_type %in% c("NK cell"))
nk_dir <- file.path(OUTPUT_DIR, "gbm6_immune", "nk")
saveRDS(nk, file.path(nk_dir, "nk.rds"))
treatment_conds <- unique(nk$treatment_cond)
for (treatment in treatment_conds[treatment_conds != "untreated_noRT"]) {
  deseq_df <- utils.deseq$find_deseq_differential_genes(
    nk,
    treatment,
    "untreated_noRT",
    "treatment_cond",
    seed = 5220
  )
  write.table(deseq_df, file = file.path(
    nk_dir, sprintf("%s_untreated_noRT.csv", treatment)
  ))
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
  ggsave(file.path(nk_dir, sprintf("%s_volcano.pdf", treatment)),
    volcano_plot,
    height = 10, width = 10, device = "pdf"
  )
}

# Ok with the immune cells in mind, let's go back to the global UMAP and
# relabel the immune cells with their appropriate labels
menv_data$annotated_cell_type_with_immune_subcategories <- menv_data$annotated_cell_type
menv_data@meta.data[rownames(immune_data@meta.data), "annotated_cell_type_with_immune_subcategories"] <- immune_data@meta.data[, "annotated_cell_type"]
umap_dimplot <- DimPlot(
  menv_data,
  reduction = "umap",
  group.by = "annotated_cell_type",
  label = TRUE,
  label.size = 2
)
ggsave(file.path(OUTPUT_DIR, "gbm6_menv", "menv_cell_type_dimplot_immune_subset.pdf"),
  plot = umap_dimplot, device = "pdf", height = 5, width = 7
)

# Now that we've decided we're going to use an inset, let's get the cell frequencies.
frequency_table <- as.data.frame(table(menv_data$annotated_cell_type, menv_data$treatment_cond))
colnames(frequency_table) <- c("Cell type", "Treatment", "Count")
write.table(frequency_table,
  file = file.path(OUTPUT_DIR, "gbm6_menv", "frequency_table_across_annotations_including_immune_inset.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Get the cell frequencies for the immune subset as well
frequency_table <- as.data.frame(table(menv_data$annotated_cell_type_with_immune_subcategories, menv_data$treatment_cond))
colnames(frequency_table) <- c("Cell type", "Treatment", "Count")
write.table(frequency_table,
  file = file.path(OUTPUT_DIR, "gbm6_menv", "frequency_table_across_annotations_including_immune_subcategories.txt"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
dimplot <- DimPlot(
  menv_data,
  reduction = "umap", group.by = "leiden_0.3", label = TRUE
)
ggsave(
  file.path(OUTPUT_DIR, "dimplot.png"),
  height = 10, width = 12
)

# Save the metadata from the immune data and from the menv data, with UMAP embeddings
menv_data_umap <- menv_data[["umap"]]@cell.embeddings
menv_data_with_embeddings <- merge(menv_data_umap, menv_data@meta.data, by = 0)
colnames(menv_data_with_embeddings)[1] <- "Cell barcodes"
menv_data_with_embeddings <- menv_data_with_embeddings[c("Cell barcodes", "umap_1", "umap_2", "annotated_cell_type")]
write.table(menv_data_with_embeddings, file.path(OUTPUT_DIR, "gbm6_menv", "menv_embeddings_and_cell_types_only.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

immune_data_umap <- immune_data[["umap"]]@cell.embeddings
immune_data_with_embeddings <- merge(immune_data_umap, immune_data@meta.data, by = 0)
colnames(immune_data_with_embeddings)[1] <- "Cell barcodes"
immune_data_with_embeddings <- immune_data_with_embeddings[c("Cell barcodes", "umap_1", "umap_2", "annotated_cell_type")]
write.table(immune_data_with_embeddings, file.path(OUTPUT_DIR, "gbm6_immune", "immune_data_metadata.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# =========================================================
# Combined UMAP plotting
# =========================================================
# For both GBM6 and GBM43, do the following:
# 1. Load the tumor and menv post-QC, pre-clustering objects
# 2. Merge them together and run the clustering pipeline
# 3. Cluster at various resolutions, then try to determine a set of cells for
#    immunogenicity

resolutions <- c(0.1, 0.2, 0.3, 0.5, 1)

tumor <- readRDS(file.path(PIPELINES_DIR, "gbm6_5", "tumor_analysis", "postQC_tumor_merged.rds"))
menv <- readRDS(file.path(PIPELINES_DIR, "gbm6_5", "menv_analysis", "postQC_menv_merged.rds"))
combined <- merge(tumor, menv)
combined <- JoinLayers(combined)
combined$species <- ifelse(combined$percent_human > 90, "human", "mouse")

# Note that because of the size of globals needed, we need to set the maximum
# allowed global size. Specify a size of around 2 GB
options(future.globals.maxSize = 3 * 1024^3)

# Run the SCTransform based pipeline
combined <- SCTransform(combined)
combined <- RunPCA(combined)
combined <- FindNeighbors(combined)
for (res in resolutions) {
  combined <- FindClusters(
    combined,
    algorithm = 1, # Use Louvain for now
    resolution = res,
    cluster.name = sprintf("louvain_%s", res)
  )
}
combined <- RunUMAP(combined, dims = 1:50)
saveRDS(combined, file.path(OUTPUT_DIR, "combined.rds"))

# Now that we've found the clusters...find markers for each?
dim_loadings <- VizDimLoadings(data, dims = 1:4, reduction = "pca")
ggsave(file.path(OUTPUT_DIR, "pca_dim_loadings.pdf"),
  plot = dim_loadings,
  device = "pdf", height = 10, width = 10
)

# Verify that the PCs do not overly cluster by sample
pc_dimplot <- DimPlot(combined, reduction = "pca", group.by = c(
  "sample", "treatment_cond", "species"
))
ggsave(file.path(OUTPUT_DIR, "pca_dim_plot_treatment_samples_species.pdf"),
  plot = pc_dimplot,
  device = "pdf", height = 10, width = 20
)

# Plot various metadata against UMAP: sample, cluster resolutions, phase,
# ribosomal gene percentage, mitochondrial gene percentage
umap_dimplot <- DimPlot(
  combined,
  reduction = "umap",
  group.by = c(sprintf("louvain_%s", resolutions), "sample", "treatment_cond", "species")
)
ggsave(file.path(OUTPUT_DIR, "umap_against_discrete_metadata.pdf"),
  plot = umap_dimplot, device = "pdf", height = 15, width = 20
)

# Find proportion of treatment_cond in clusters and clusters in treatment_cond
for (res in RESOLUTIONS) {
  cluster_res <- sprintf("leiden_%s", res)
  proportions <- data@meta.data %>%
    group_by(treatment_cond, !!sym(cluster_res)) %>%
    dplyr::summarize(intersect_counts = n()) %>%
    group_by(treatment_cond) %>%
    mutate(proportions = intersect_counts / sum(intersect_counts))
  proportion_plot <- ggplot(proportions, aes(x = treatment_cond, y = proportions, fill = !!sym(cluster_res))) +
    geom_bar(stat = "identity", position = "fill") +
    labs(
      title = "Proportion of Leiden Clusters by Treatment Condition",
      x = "Treatment Condition",
      y = "Proportion",
      fill = "Leiden Cluster"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(file.path(TUMOR_OUTPUT_DIR, sprintf("%s_proportion_plot.pdf", cluster_res)),
    plot = proportion_plot, height = 10, width = 30
  )
}


# Find markers for the clusters
for (res in RESOLUTIONS) {
  Idents(data) <- sprintf("leiden_%s", res)
  df <- FindAllMarkers(
    data
  )
  write.table(df, file = file.path(TUMOR_OUTPUT_DIR, sprintf("leiden_%s_markers", res)))

  df %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5) %>%
    slice_head(n = 20) %>%
    ungroup() -> top20
  marker_heatmap <- DoHeatmap(data, features = top20$gene) + NoLegend()
  ggsave(file.path(TUMOR_OUTPUT_DIR, sprintf("leiden_%s_marker_heatmap.pdf", res)),
    plot = marker_heatmap, device = "pdf", height = 20, width = 20
  )
}

saveRDS(data, file.path(TUMOR_OUTPUT_DIR, "postClustering_tumor_merged.rds"))

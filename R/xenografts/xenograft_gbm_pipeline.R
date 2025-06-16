# =========================================================
# Library imports and configuration
# =========================================================

lang <- Sys.getenv("LANG")
lc_type <- Sys.getenv("LC_CTYPE")

cat("LANG: ", lang, "\n")
cat("LC_CTYPE: ", lc_type, "\n")

library(Seurat)
library(SoupX)
library(ggplot2)
library(RColorBrewer)
library(SCpubr)
library(dplyr)
library(scMRMA)

set.seed(5220)
args <- commandArgs(trailingOnly = TRUE)

utils.pipelines <- new.env()
source("R/utils/pipelines.R", local = utils.pipelines)

utils.deseq <- new.env()
source("R/utils/deseq.R", local = utils.deseq)

RUN_ID <- "gbm6_5"
BASE_OUTPUT_DIR <- sprintf("output/xenografts/pipelines/%s", RUN_ID)
utils.pipelines$create_dir(BASE_OUTPUT_DIR)
MDATA_FILE <- file.path(BASE_OUTPUT_DIR, "metadata.txt")

utils.pipelines$save_to_metadata_file(
  MDATA_FILE, sprintf("Starting pipeline run at %s with ID %s", Sys.time(), RUN_ID)
)

EXTERNAL_DATA_IMPORT_DIR <- "data/xenograft_pipeline_data"
DATA_DIR <- "/raleighlab/data1/liuj/gbm_perturb/cellranger_out_combined"
EXPERIMENT_DIR <- "gbm_pdx_dnapk_gbm43"

# Parameters
RUN_SOUPX <- if (args[2] == "TRUE") TRUE else FALSE
MT_CUTOFF <- if (args[3] == "NULL") NULL else as.numeric(args[3])
REGRESS_RIBOSOMAL_GENES <- if (args[4] == "TRUE") TRUE else FALSE

# =========================================================
# Data imports
# =========================================================
import_from_dir <- file.path(DATA_DIR, EXPERIMENT_DIR)
utils.pipelines$save_to_metadata_file(
  MDATA_FILE,
  sprintf("Retrieving data from %s", import_from_dir)
)

raw_seurat_objects <- list()
for (sample_name in list.files(import_from_dir)) {
  print(sprintf("Processing sample %s", sample_name))
  count_matrix <- Read10X(
    file.path(import_from_dir, sample_name, "outs/filtered_feature_bc_matrix"),
  )
  raw_seurat_objects[[sample_name]] <- CreateSeuratObject(
    counts = count_matrix
  )
}

# =========================================================
# Run QC
# =========================================================
# We will run QC on every sample individually. QC steps:
# - Run SoupX to remove ambient RNA on the entire dataset
# - Generate human and mouse Seurat objects
# - Filter for mitochondrial counts and nCount_RNA

QC_OUTPUT_DIR <- file.path(BASE_OUTPUT_DIR, "qc")
utils.pipelines$create_dir(QC_OUTPUT_DIR)
MIN_FEATURES_PER_CELL <- 200
MIN_CELLS_PER_FEATURE <- 3
utils.pipelines$save_to_metadata_file(
  MDATA_FILE,
  sprintf(
    "Starting QC with:
                - output dir %s
                - run soupX %s
                - MT cutoff %s
                - Min features per cell %s
                - Min cells per feature %s
            ",
    QC_OUTPUT_DIR,
    RUN_SOUPX,
    ifelse(is.null(MT_CUTOFF), "NULL", MT_CUTOFF),
    MIN_FEATURES_PER_CELL,
    MIN_CELLS_PER_FEATURE
  )
)

run_soupX <- function(path_to_data, output_dir, use_auto_contamination = TRUE,
                      gene_label = NULL, gene_list = NULL) {
  utils.pipelines$create_dir(output_dir)

  # Read in the raw matrix data. Note that this has clusters already
  sc <- load10X(path_to_data)

  # Estimate the contamination from the data
  if (use_auto_contamination) {
    sc <- autoEstCont(sc)
    ggsave(file.path(output_dir, "automatic_contamination_estimation.png"), height = 8, width = 10)
  } else {
    if (is.null(gene_list) || is.null(gene_label)) {
      utils.pipelines$save_to_metadata_file(
        MDATA_FILE,
        "Terminated at SoupX due to no gene list or gene label."
      )
      return(NULL)
    }
    useToEst <- estimateNonExpressingCells(
      sc,
      nonExpressedGeneList = gene_list
    )
    marker_map_plot <- plotMarkerMap(
      sc,
      geneSet = list(gene_label = gene_list),
      useToEst = useToEst
    )
    ggsave(file.path(output_dir, sprintf("%s_marker_map_plot.pdf", gene_label)), device = "pdf", height = 8, width = 10)
    sc <- calculateContaminationFraction(sc, gene_list, useToEst = useToEst)
  }

  # Save the SoupX object
  saveRDS(sc, file = file.path(output_dir, "soupx_object.rds"))

  # Adjust the counts
  out <- adjustCounts(sc, roundToInt = TRUE)

  # Return a Seurat object
  return(CreateSeuratObject(out))
}

tumor_seurat_objects <- list()
microenvironment_seurat_objects <- list()
for (sample_name in names(raw_seurat_objects)) {
  sample_out_directory <- file.path(QC_OUTPUT_DIR, sample_name)
  utils.pipelines$create_dir(sample_out_directory)

  # Fetch the Seurat object for this sample
  seurat_object <- raw_seurat_objects[[sample_name]]

  # Run SoupX to remove ambient RNA
  print(sample_name)
  if (RUN_SOUPX) {
    utils.pipelines$save_to_metadata_file(
      MDATA_FILE,
      sprintf("Running SoupX for %s", sample_name)
    )
    seurat_object <- run_soupX(
      file.path(DATA_DIR, EXPERIMENT_DIR, sample_name, "outs"),
      file.path(sample_out_directory, "soupX")
    )
  }

  # Filter features to min cells per feature
  counts <- GetAssayData(seurat_object, assay = "RNA", layer = "counts")
  seurat_object <- CreateSeuratObject(
    counts = counts,
    min.cells = MIN_CELLS_PER_FEATURE
  )

  # Split into human and mouse objects
  seurat_object[["percent_human"]] <- PercentageFeatureSet(
    seurat_object,
    pattern = "GRCh38-"
  )
  seurat_object[["percent_mouse"]] <- PercentageFeatureSet(
    seurat_object,
    pattern = "mm10---"
  )
  human_genes <- rownames(seurat_object)[grep("GRCh38-", rownames(seurat_object))]
  mouse_genes <- rownames(seurat_object)[grep("mm10---", rownames(seurat_object))]
  tumor <- subset(seurat_object, percent_human > 90, features = human_genes)
  menv <- subset(seurat_object, percent_mouse > 90, features = mouse_genes)

  # Filter for min cells per gene and min genes per cell
  tumor <- subset(tumor, nFeature_RNA > 200)
  menv <- subset(menv, nFeature_RNA > 200)

  # Calculate mitochondrial counts
  tumor[["percent_mt"]] <- PercentageFeatureSet(tumor, pattern = "GRCh38-MT-")
  menv[["percent_mt"]] <- PercentageFeatureSet(menv, pattern = "mm10---mt-")

  # Plot the mitochondrial counts ahead of filtering
  tumor_mt_plot <- VlnPlot(tumor, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3)
  menv_mt_plot <- VlnPlot(menv, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), ncol = 3)
  ggsave(file.path(sample_out_directory, "tumor_mt_plot.pdf"),
    plot = tumor_mt_plot, device = "pdf", height = 10, width = 10
  )
  ggsave(file.path(sample_out_directory, "menv_mt_plot.pdf"),
    plot = menv_mt_plot, device = "pdf", height = 10, width = 10
  )

  # Mitochondrial filters (warning MT_CUTOFF is a percent!)
  if (!is.null(MT_CUTOFF)) {
    menv <- tryCatch(
      {
        subset(menv, percent_mt < MT_CUTOFF)
      },
      error = function(e) {
        utils.pipelines$save_to_metadata_file(
          MDATA_FILE, sprintf("No cells; skipping sample %s", sample_name)
        )
        NULL
      }
    )
    tumor <- tryCatch(
      {
        subset(tumor, percent_mt < MT_CUTOFF)
      },
      error = function(e) {
        utils.pipelines$save_to_metadata_file(
          MDATA_FILE, sprintf("No cells; skipping sample %s", sample_name)
        )
        NULL
      }
    )
  }

  treatment_cond_mappings <- read.csv(file.path(
    EXTERNAL_DATA_IMPORT_DIR,
    "sample_treatment_cond_mappings.csv"
  ))
  rownames(treatment_cond_mappings) <- treatment_cond_mappings$sample

  if (!is.null(tumor)) {
    rownames(tumor) <- gsub("GRCh38-", "", rownames(tumor))
    tumor$sample <- sample_name
    tumor$treatment_cond <- treatment_cond_mappings[sample_name, "treatment_cond"]
  }
  tumor_seurat_objects[[sample_name]] <- tumor

  if (!is.null(menv)) {
    rownames(menv) <- gsub("mm10---", "", rownames(menv))
    menv$sample <- sample_name
    menv$treatment_cond <- treatment_cond_mappings[sample_name, "treatment_cond"]
  }
  microenvironment_seurat_objects[[sample_name]] <- menv
}

# =========================================================
# Sample-integrated analysis of tumor and microenvironments
#
# Tumor and microenvironment outputs will be aggregated separately,
# so we start defining tumor and menv output directories at this point.
# =========================================================
# Merge all of the samples together and give them unique barcodes.
TUMOR_OUTPUT_DIR <- file.path(BASE_OUTPUT_DIR, "tumor_analysis")
MENV_OUTPUT_DIR <- file.path(BASE_OUTPUT_DIR, "menv_analysis")
utils.pipelines$create_dir(TUMOR_OUTPUT_DIR)
utils.pipelines$create_dir(MENV_OUTPUT_DIR)

utils.pipelines$save_to_metadata_file(
  MDATA_FILE,
  "Integrating and saving Seurat objects after running QC"
)

tumor <- merge(
  tumor_seurat_objects[[1]],
  y = tumor_seurat_objects[-1],
  add.cell.ids = names(tumor_seurat_objects)
)
tumor <- JoinLayers(tumor)
saveRDS(tumor, file.path(TUMOR_OUTPUT_DIR, "postQC_tumor_merged.rds"))

menv <- merge(
  microenvironment_seurat_objects[[1]],
  y = microenvironment_seurat_objects[-1],
  add.cell.ids = names(microenvironment_seurat_objects)
)
menv <- JoinLayers(menv)
saveRDS(menv, file.path(MENV_OUTPUT_DIR, "postQC_menv_merged.rds"))

# =========================================================
# Tumor analysis
#
# Normalize and calculate regressions for ribosomal and
# cell cycle scores
# =========================================================
data <- tumor

REGRESS_CELL_CYCLE <- FALSE
REGRESS_CELL_CYCLE_DIFFERENCE <- FALSE

utils.pipelines$save_to_metadata_file(
  MDATA_FILE,
  sprintf(
    "Running SCTransform with ribosomal regression: %s,
            cell cycle regression: %s, cell cycle diff regression: %s",
    REGRESS_RIBOSOMAL_GENES,
    REGRESS_CELL_CYCLE,
    REGRESS_CELL_CYCLE_DIFFERENCE
  )
)

# Calculate ribosomal gene percentages
data[["percent_ribosomal"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")
rb_plot <- VlnPlot(data, features = c("percent_ribosomal"), ncol = 1)
ggsave(file.path(TUMOR_OUTPUT_DIR, "ribosomal_percentage_plot.pdf"),
  device = "pdf",
  height = 10, width = 10
)

# Normalize the data, regressing out ribosomal variation if necessary
if (REGRESS_RIBOSOMAL_GENES) {
  vars <- c("percent_ribosomal")
} else {
  vars <- NULL
}
data <- SCTransform(data, vars.to.regress = vars, ncells = ncol(data))

# Calculate cell cycle scores
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
data <- CellCycleScoring(data, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
data$cc_diff <- data$S.Score - data$G2M.Score

if (REGRESS_CELL_CYCLE) {
  DefaultAssay(data) <- "RNA"
  data <- SCTransform(data, vars.to.regress = c("S.Score", "G2M.Score"), ncells = ncol(data))
} else if (REGRESS_CELL_CYCLE_DIFFERENCE) {
  DefaultAssay(data) <-
    data <- SCTransform(data, vars.to.regress = "cc_diff", ncells = ncol(data))
}

saveRDS(data, file.path(TUMOR_OUTPUT_DIR, "postNormalization_tumor_merged.rds"))

# =========================================================
# Tumor analysis
#
# Cluster and find markers at various resolutions
# =========================================================
RESOLUTIONS <- c(0.1, 0.2, 0.3, 0.5, 1)

utils.pipelines$save_to_metadata_file(
  MDATA_FILE,
  sprintf(
    "Clustering and finding markers with resolutions %s",
    paste(RESOLUTIONS, collapse = ", ")
  )
)

data <- RunPCA(data, seed.use = 5220)
data <- FindNeighbors(data, k.param = 20)
for (res in RESOLUTIONS) {
  data <- FindClusters(
    data,
    resolution = res,
    random.seed = 5220,
    method = "igraph",
    cluster.name = sprintf("leiden_%s", res)
  )
}
data <- RunUMAP(data, dims = 1:50)

# Output the PCA loading genes
dim_loadings <- VizDimLoadings(data, dims = 1:4, reduction = "pca")
ggsave(file.path(TUMOR_OUTPUT_DIR, "pca_dim_loadings.pdf"),
  plot = dim_loadings,
  device = "pdf", height = 10, width = 10
)

# Verify that the PCs do not overly cluster by sample
pc_dimplot <- DimPlot(data, reduction = "pca", group.by = c(
  "sample", "treatment_cond"
))
ggsave(file.path(TUMOR_OUTPUT_DIR, "pca_dim_plot_treatment_samples.pdf"),
  plot = pc_dimplot,
  device = "pdf", height = 10, width = 20
)
pc_featureplot <- FeaturePlot(data, reduction = "pca", features = c(
  "percent_mt", "percent_ribosomal"
))
ggsave(file.path(TUMOR_OUTPUT_DIR, "pca_dim_plot_qc.pdf"),
  plot = pc_featureplot,
  device = "pdf", height = 10, width = 20
)

# Get an idea for what the PCs are doing
pc_dimheatmap <- DimHeatmap(data,
  dims = 1:15, cells = 500, balanced = TRUE,
  fast = FALSE
)
pc_dimheatmap <- lapply(1:length(pc_dimheatmap), function(x) {
  plot.i <- pc_dimheatmap[[x]] + theme(legend.position = "none") + ggtitle(paste0("PC", x))
  return(plot.i)
})
pc_dimheatmap <- patchwork::wrap_plots(pc_dimheatmap,
  ncol = 3
)
ggsave(file.path(TUMOR_OUTPUT_DIR, "pca_dim_heatmap.pdf"),
  plot = pc_dimheatmap,
  device = "pdf", height = 20, width = 10
)

# Plot various metadata against UMAP: sample, cluster resolutions, phase,
# ribosomal gene percentage, mitochondrial gene percentage
umap_dimplot <- DimPlot(
  data,
  reduction = "umap",
  group.by = c(sprintf("leiden_%s", RESOLUTIONS), "sample", "treatment_cond")
)
ggsave(file.path(TUMOR_OUTPUT_DIR, "umap_against_discrete_metadata.pdf"),
  plot = umap_dimplot, device = "pdf", height = 15, width = 20
)

umap_featureplot <- FeaturePlot(
  data,
  reduction = "umap",
  features = c("percent_mt", "percent_ribosomal", "cc_difference", "S.Score", "G2M.Score"),
) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(file.path(TUMOR_OUTPUT_DIR, "umap_against_continuous_metadata.pdf"),
  plot = umap_featureplot, device = "pdf", height = 15, width = 20
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

# =========================================================
# Tumor analysis
#
# DESeq
# =========================================================
DESEQ_OUTPUT_DIR <- file.path(TUMOR_OUTPUT_DIR, "deseq")
utils.pipelines$create_dir(DESEQ_OUTPUT_DIR)

all_conditions <- unique(data$treatment_cond)
conditions <- all_conditions[all_conditions != "untreated_noRT"]

utils.pipelines$save_to_metadata_file(
  MDATA_FILE,
  sprintf(
    "Running DESeq2 on tumor data with conditions: %s",
    paste(conditions, collapse = ", ")
  )
)

DefaultAssay(data) <- "RNA"
for (condition in conditions) {
  print("Starting deseq analysis")
  deseq_df <- utils.deseq$find_deseq_differential_genes(
    data,
    condition,
    "untreated_noRT",
    "treatment_cond",
    seed = 5220
  )
  print("Finished deseq analysis")
  write.table(deseq_df, file.path(
    DESEQ_OUTPUT_DIR,
    sprintf("%s_untreated_noRT.csv", condition)
  ))
  print("Finished writing table")
}
deseq_df <- utils.deseq$find_deseq_differential_genes(
  data,
  "nedisertib_RT",
  "untreated_RT",
  "treatment_cond",
  seed = 5220
)
print("Finished deseq analysis")
write.table(deseq_df, file.path(
  DESEQ_OUTPUT_DIR,
  sprintf("%s_untreated_RT.csv", "nedisertib_RT")
))
print("Finished writing table")

DefaultAssay(data) <- "SCT"

# =========================================================
# Transition
# =========================================================
utils.pipelines$save_to_metadata_file(
  MDATA_FILE,
  sprintf("Finished running tumor analysis. Transitioning to microenv analysis")
)

# =========================================================
# Microenvironment analysis
#
# Normalize and calculate regressions for ribosomal and
# cell cycle scores
# =========================================================
data <- menv

REGRESS_CELL_CYCLE <- FALSE
REGRESS_CELL_CYCLE_DIFFERENCE <- FALSE

utils.pipelines$save_to_metadata_file(
  MDATA_FILE,
  sprintf(
    "Running SCTransform with ribosomal regression: %s,
            cell cycle regression: %s, cell cycle diff regression: %s",
    REGRESS_RIBOSOMAL_GENES,
    REGRESS_CELL_CYCLE,
    REGRESS_CELL_CYCLE_DIFFERENCE
  )
)

# Calculate ribosomal gene percentages
data[["percent_ribosomal"]] <- PercentageFeatureSet(data, pattern = "^Rp[sl]")
rb_plot <- VlnPlot(data, features = c("percent_ribosomal"), ncol = 1)
ggsave(file.path(MENV_OUTPUT_DIR, "ribosomal_percentage_plot.pdf"),
  device = "pdf",
  height = 10, width = 10
)

# Normalize the data, regressing out ribosomal variation if necessary
if (REGRESS_RIBOSOMAL_GENES) {
  vars <- c("percent_ribosomal")
} else {
  vars <- NULL
}
data <- SCTransform(data, vars.to.regress = vars, ncells = ncol(data))

# Calculate cell cycle scores
s_genes <- tools::toTitleCase(tolower(cc.genes$s.genes))
g2m_genes <- tools::toTitleCase(tolower(cc.genes$g2m.genes))
data <- CellCycleScoring(data, s.features = s_genes, g2m.features = g2m_genes, set.ident = TRUE)
data$cc_diff <- data$S.Score - data$G2M.Score

if (REGRESS_CELL_CYCLE) {
  DefaultAssay(data) <- "RNA"
  data <- SCTransform(data, vars.to.regress = c("S.Score", "G2M.Score"), ncells = ncol(data))
} else if (REGRESS_CELL_CYCLE_DIFFERENCE) {
  DefaultAssay(data) <-
    data <- SCTransform(data, vars.to.regress = "cc_diff", ncells = ncol(data))
}

saveRDS(data, file.path(MENV_OUTPUT_DIR, "postNormalization_menv_merged.rds"))

# =========================================================
# Microenvironment analysis
#
# Cluster and find markers at various resolutions
# =========================================================
RESOLUTIONS <- c(0.1, 0.2, 0.3, 0.5, 1)

utils.pipelines$save_to_metadata_file(
  MDATA_FILE,
  sprintf(
    "Clustering and finding markers with resolutions %s",
    paste(RESOLUTIONS, collapse = ", ")
  )
)

data <- RunPCA(data, seed.use = 5220)
data <- FindNeighbors(data, k.param = 20)
for (res in RESOLUTIONS) {
  data <- FindClusters(
    data,
    resolution = res,
    random.seed = 5220,
    method = "igraph",
    cluster.name = sprintf("leiden_%s", res)
  )
}
data <- RunUMAP(data, dims = 1:50)

# Output the PCA loading genes
dim_loadings <- VizDimLoadings(data, dims = 1:4, reduction = "pca")
ggsave(file.path(MENV_OUTPUT_DIR, "pca_dim_loadings.pdf"),
  plot = dim_loadings,
  device = "pdf", height = 10, width = 10
)

# Verify that the PCs do not overly cluster by sample
pc_dimplot <- DimPlot(data, reduction = "pca", group.by = c(
  "sample", "treatment_cond"
))
ggsave(file.path(MENV_OUTPUT_DIR, "pca_dim_plot_treatment_samples.pdf"),
  plot = pc_dimplot,
  device = "pdf", height = 10, width = 20
)
pc_featureplot <- FeaturePlot(data, reduction = "pca", features = c(
  "percent_mt", "percent_ribosomal"
))
ggsave(file.path(MENV_OUTPUT_DIR, "pca_dim_plot_qc.pdf"),
  plot = pc_featureplot,
  device = "pdf", height = 10, width = 20
)

# Get an idea for what the PCs are doing
pc_dimheatmap <- DimHeatmap(data,
  dims = 1:15, cells = 500, balanced = TRUE,
  fast = FALSE
)
pc_dimheatmap <- lapply(1:length(pc_dimheatmap), function(x) {
  plot.i <- pc_dimheatmap[[x]] + theme(legend.position = "none") + ggtitle(paste0("PC", x))
  return(plot.i)
})
pc_dimheatmap <- patchwork::wrap_plots(pc_dimheatmap,
  ncol = 3
)
ggsave(file.path(MENV_OUTPUT_DIR, "pca_dim_heatmap.pdf"),
  plot = pc_dimheatmap,
  device = "pdf", height = 20, width = 10
)

# Plot various metadata against UMAP: sample, cluster resolutions, phase,
# ribosomal gene percentage, mitochondrial gene percentage
umap_dimplot <- DimPlot(
  data,
  reduction = "umap",
  group.by = c(sprintf("leiden_%s", RESOLUTIONS), "sample", "treatment_cond")
)
ggsave(file.path(MENV_OUTPUT_DIR, "umap_against_discrete_metadata.pdf"),
  plot = umap_dimplot, device = "pdf", height = 15, width = 20
)

umap_featureplot <- FeaturePlot(
  data,
  reduction = "umap",
  features = c("percent_mt", "percent_ribosomal", "cc_difference", "S.Score", "G2M.Score"),
) & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
ggsave(file.path(MENV_OUTPUT_DIR, "umap_against_continuous_metadata.pdf"),
  plot = umap_featureplot, device = "pdf", height = 15, width = 20
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
  ggsave(file.path(MENV_OUTPUT_DIR, sprintf("%s_proportion_plot.pdf", cluster_res)),
    plot = proportion_plot, height = 10, width = 30
  )
}

# Find markers for the clusters
for (res in RESOLUTIONS) {
  Idents(data) <- sprintf("leiden_%s", res)
  df <- FindAllMarkers(
    data
  )
  write.table(df, file = file.path(MENV_OUTPUT_DIR, sprintf("leiden_%s_markers", res)))

  df %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 0.5) %>%
    slice_head(n = 20) %>%
    ungroup() -> top20
  marker_heatmap <- DoHeatmap(data, features = top20$gene) + NoLegend()
  ggsave(file.path(MENV_OUTPUT_DIR, sprintf("leiden_%s_marker_heatmap.pdf", res)),
    plot = marker_heatmap, device = "pdf", height = 20, width = 20
  )
}

saveRDS(data, file.path(MENV_OUTPUT_DIR, "postClustering_menv_merged.rds"))

# =========================================================
# Microenvironment analysis
#
# Run scMRMA to line up the clusters with cell types using
# the known clusters and with no bias
# =========================================================
utils.pipelines$save_to_metadata_file(
  MDATA_FILE,
  "Running scMRMA on microenvironment data"
)

# Try with no bias
scMRMA_result <- scMRMA(data, species = "Mm")
data[["scMRMA_unbiased"]] <- scMRMA_result$multiR$annotationResult[
  colnames(data),
  ncol(scMRMA_result$multiR$annotationResult)
]

# Try with known clusters
for (res in RESOLUTIONS) {
  Idents(data) <- sprintf("leiden_%s", res)
  scMRMA_result <- scMRMA(data, species = "Mm", selfClusters = Idents(data))
  data[[sprintf("scMRMA_leiden_%s", res)]] <- scMRMA_result$multiR$annotationResult[
    colnames(data),
    ncol(scMRMA_result$multiR$annotationResult)
  ]
}

# Save an scMRMA dim plot
umap_dimplot <- DimPlot(
  data,
  reduction = "umap",
  group.by = c(sprintf("scMRMA_leiden_%s", RESOLUTIONS), "scMRMA_unbiased"),
  label = TRUE,
  label.size = 2,
  label.box = TRUE
) & NoLegend()
ggsave(file.path(MENV_OUTPUT_DIR, "umap_against_scMRMA_annotations.pdf"),
  plot = umap_dimplot, device = "pdf", height = 12, width = 20
)

saveRDS(data, file.path(MENV_OUTPUT_DIR, "postSCMRMA_microenvironment_merged.rds"))

# =========================================================
# Microenvironment analysis
#
# DESeq
# =========================================================
DESEQ_OUTPUT_DIR <- file.path(MENV_OUTPUT_DIR, "deseq")
utils.pipelines$create_dir(DESEQ_OUTPUT_DIR)

all_conditions <- unique(data$treatment_cond)
conditions <- all_conditions[all_conditions != "untreated_noRT"]

utils.pipelines$save_to_metadata_file(
  MDATA_FILE,
  sprintf(
    "Running DESeq2 on microenvironment data with conditions: %s",
    paste(conditions, collapse = ", ")
  )
)

DefaultAssay(data) <- "RNA"
for (condition in conditions) {
  deseq_df <- utils.deseq$find_deseq_differential_genes(
    data,
    condition,
    "untreated_noRT",
    "treatment_cond",
    seed = 5220
  )
  write.table(deseq_df, file.path(
    DESEQ_OUTPUT_DIR,
    sprintf("%s_untreated_noRT.csv", condition)
  ))
}
DefaultAssay(data) <- "SCT"

# Leaving space here to decide which cluster resolution works best for DESeq
# comparisons. But we shouldn't need to run this at extremely high throughput.

# =========================================================
# Finish
# =========================================================
utils.pipelines$save_to_metadata_file(
  MDATA_FILE,
  sprintf("Run finished at %s.", Sys.time())
)

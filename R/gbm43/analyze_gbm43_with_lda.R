# =========================================================
# GBM43 in vivo perturb-seq analysis with Mixscape and LDA
# Refactored 2025-04-22
# =========================================================

# ==================
# Library imports
# ==================
library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)
library(knitr)
library(msigdbr)
library(fgsea)
library(dplyr)
library(data.table)
library(ggplot2)
library(grid)
library(patchwork)
library(scales)
library(reshape2)
library(SCpubr)
library(RColorBrewer)

OUTPUT_DIR <- "output/gbm43/lda"

# ==================
# Data import
# ==================
data <- readRDS("/raleighlab/data1/liuj/gbm_perturb/analysis/GBM43_1_malignant_only_annotated_20230817.Rds")
off_targets <- c("GRAMD4", "FBXL16", "OVCH1")
include_genes <- setdiff(unique(data$sgRNA), off_targets)
data <- subset(data, sgRNA %in% include_genes)

# ==================
# Normalization
# ==================
DefaultAssay(data) <- "RNA"
data <- NormalizeData(data) %>%
  FindVariableFeatures() %>%
  ScaleData()
data <- RunPCA(data)
ElbowPlot(data, ndims = 50)
data <- RunUMAP(data, dims = 1:40)

# ==================
# Cell cycle scoring
# ==================
cc.genes <- Seurat::cc.genes
s.genes <- paste0("GRCh38-", cc.genes$s.genes)
g2m.genes <- paste0("GRCh38-", cc.genes$g2m.genes)
data <- CellCycleScoring(
  data,
  s.features = s.genes,
  g2m.features = g2m.genes,
  set.ident = TRUE
)

# ==================
# Global plots
# ==================
plot <- DimPlot(
  data,
  group.by = "Phase",
  split.by = "cond",
  reduction = "umap",
  ncol = 1,
  pt.size = 3,
  label = TRUE,
  label.box = TRUE,
  label.color = "white"
)

plot1 <- DimPlot(
  data,
  group.by = "sgRNA",
  split.by = "cond",
  reduction = "umap",
  ncol = 1,
  pt.size = 3,
  label = TRUE,
  label.box = TRUE,
  label.color = "white"
)

ggsave(
  sprintf("%s/umap_cellcycle.png", OUTPUT_DIR),
  plot | plot1,
  width = 40,
  height = 20
)

for (i in c("non-targeting", "PRKDC", "CENPT")) {
  perturb_data <- subset(data, sgRNA == i)
  perturb_cells <- colnames(perturb_data)
  plot <- DimPlot(
    data,
    group.by = "Phase",
    split.by = "cond",
    reduction = "umap",
    ncol = 1,
    pt.size = 3
  )
  plot1 <- DimPlot(
    data,
    group.by = "sgRNA",
    split.by = "cond",
    reduction = "umap",
    cells.highlight = perturb_cells,
    ncol = 1,
    pt.size = 3
  )
  ggsave(
    sprintf("%s/umap_cellcycle_highlight_%s.png", OUTPUT_DIR, i),
    plot | plot1,
    width = 40,
    height = 20
  )
}

# ==================
# Perturbation signal and Mixscape
# ==================
data <- CalcPerturbSig(
  data,
  assay = "RNA",
  slot = "data",
  gd.class = "sgRNA",
  nt.cell.class = "non-targeting",
  reduction = "pca",
  ndims = 50,
  num.neighbors = 30,
  new.assay.name = "PRTB",
  split.by = "cond"
)

DefaultAssay(data) <- "PRTB"
VariableFeatures(data) <- VariableFeatures(data[["RNA"]])
data <- ScaleData(data, do.scale = FALSE, do.center = TRUE)
data <- RunPCA(data, reduction.key = "prtbpca", reduction.name = "prtbpca")
data <- RunUMAP(data, dims = 1:40, reduction = "prtbpca", reduction.key = "prtbumap", reduction.name = "prtbumap")

plot <- DimPlot(
  data,
  group.by = "Phase",
  pt.size = 2,
  split.by = "cond",
  reduction = "prtbumap",
  ncol = 1
)

plot1 <- DimPlot(
  data,
  group.by = "sgRNA",
  pt.size = 2,
  split.by = "cond",
  reduction = "prtbumap",
  ncol = 1,
  label = TRUE,
  label.box = TRUE,
  label.color = "white"
)

ggsave(
  sprintf("%s/umap_prtb_cellcycle.png", OUTPUT_DIR),
  plot | plot1,
  width = 40,
  height = 20
)

# ==================
# Mixscape call distributions
# ==================
data <- RunMixscape(
  data,
  assay = "PRTB",
  slot = "scale.data",
  labels = "sgRNA",
  nt.class.name = "non-targeting",
  min.de.genes = 5,
  logfc.threshold = 0.1,
  iter.num = 10,
  de.assay = "RNA",
  verbose = TRUE,
  prtb.type = "KO",
  split.by = "cond"
)

df <- melt(table(data$mixscape_class.global, data$mixscape_class))
df$Var2 <- as.character(df$Var2)
df$Var1 <- factor(df$Var1, levels = c("non-targeting", "NP", "KO"))
df$gene <- sapply(as.character(df$Var2), function(x) strsplit(x, split = " ")[[1]][1])
df <- df[!(df$gene == "non-targeting"), ]
top_genes <- df[df$Var1 == "KO", ] %>%
  arrange(desc(value)) %>%
  pull(Var2)
df$Var2 <- factor(df$Var2, levels = top_genes)

plot <- ggplot(df, aes(x = gene, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("grey49", "grey79", "coral1")) +
  theme_classic() +
  ylab("n cells") +
  xlab("sgRNA") +
  facet_wrap(vars(gene), ncol = 5, scales = "free") +
  labs(fill = "mixscape class") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 8),
    strip.text = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )

ggsave(
  sprintf("%s/mixscape_prtb_ncell_calls.png", OUTPUT_DIR),
  plot,
  width = 12,
  height = 8
)

saveRDS(data, file = sprintf("%s/mixscaped_data.rds", OUTPUT_DIR))

perturb_plot <- PlotPerturbScore(
  object = data,
  target.gene.ident = "PRKDC",
  target.gene.class = "sgRNA",
  mixscape.class = "mixscape_class",
  col = "coral2",
  split.by = "cond"
) + labs(fill = "mixscape class")

ggsave(
  sprintf("%s/prkdc_perturb_plot.png", OUTPUT_DIR),
  perturb_plot,
  width = 12,
  height = 8
)

# For remaining plots, split by RT and noRT
for (condition in c("noRT", "RT")) {
  cond_object <- subset(data, cond == condition)
  vln_plot <- VlnPlot(
    cond_object,
    "mixscape_class_p_ko",
    idents = c("non-targeting", "PRKDC KO", "PRKDC NP")
  ) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text = element_text(size = 16),
      plot.title = element_text(size = 20)
    ) +
    NoLegend() +
    ggtitle("mixscape posterior probabilities")

  ggsave(
    sprintf("%s/prkdc_vln_plot_%s.png", OUTPUT_DIR, condition),
    vln_plot,
    width = 12,
    height = 8
  )

  Idents(object = cond_object) <- "sgRNA"
  genepass <- c("PRKDC", "CENPT", "RBX1")

  for (i in seq_along(genepass)) {
    plot <- MixscapeHeatmap(
      object = cond_object,
      ident.1 = "non-targeting",
      ident.2 = genepass[i],
      balanced = FALSE,
      assay = "RNA",
      max.genes = 30,
      angle = 0,
      group.by = "mixscape_class",
      max.cells.group = 300,
      size = 6.5
    ) + NoLegend() +
      theme(axis.text.y = element_text(size = 16))

    ggsave(
      sprintf("%s/mixscape_prtb_%s_%s_DEgenes.png", OUTPUT_DIR, genepass[i], condition),
      plot,
      width = 15,
      height = 10
    )
  }
}

for (condition in c("noRT", "RT")) {
  cond_object <- subset(data, cond == condition)
  cond_object <- MixscapeLDA(
    object = cond_object,
    assay = "RNA",
    pc.assay = "PRTB",
    labels = "sgRNA",
    nt.label = "non-targeting",
    npcs = 10,
    logfc.threshold = 0.1,
    verbose = TRUE
  )

  cond_object <- RunUMAP(
    object = cond_object,
    dims = 1:10,
    reduction = "lda",
    reduction.key = "ldaumap",
    reduction.name = "ldaumap"
  )

  plot <- DimPlot(
    object = cond_object,
    group.by = "Phase",
    pt.size = 1,
    reduction = "ldaumap",
    ncol = 1
  )

  plot1 <- DimPlot(
    object = cond_object,
    group.by = "sgRNA",
    pt.size = 1,
    reduction = "ldaumap",
    ncol = 1,
    label = TRUE,
    label.box = TRUE,
    label.color = "white"
  )

  ggsave(
    sprintf("%s/mixscape_prtb_LDAumap_%s.png", OUTPUT_DIR, condition),
    plot | plot1,
    width = 30,
    height = 20,
    limitsize = FALSE
  )

  iter <- names(table(cond_object$sgRNA))
  for (i in seq_along(table(cond_object$sgRNA))) {
    plot <- DimPlot(
      object = cond_object,
      group.by = "sgRNA",
      pt.size = 1,
      reduction = "ldaumap",
      ncol = 1,
      cells.highlight = colnames(cond_object)[cond_object$sgRNA == iter[i]]
    ) +
      NoLegend() +
      ggtitle(iter[i])

    ggsave(
      sprintf("%s/mixscape_prtb_LDAumap_%s_%s.png", OUTPUT_DIR, iter[i], condition),
      plot,
      width = 10,
      height = 8
    )
  }

  # Make custom plot highlighting only select genes
  genepass <- c("PRKDC", "CENPT", "RBX1", "non-targeting")
  cond_object$sgRNAPass <- cond_object$sgRNA
  cond_object$sgRNAPass[!(cond_object$sgRNAPass %in% genepass)] <- "Other"

  # Set colors for each perturbation
  col <- brewer.pal(n = 4, name = "Dark2")
  col <- c(col[1:2], "gray80", "gray90", col[3:4])

  Idents(object = cond_object) <- "sgRNAPass"
  plot <- DimPlot(
    object = cond_object,
    group.by = "Phase",
    pt.size = 3,
    reduction = "ldaumap",
    ncol = 1
  )

  plot1 <- DimPlot(
    object = cond_object,
    group.by = "sgRNAPass",
    pt.size = 3,
    reduction = "ldaumap",
    ncol = 1
  ) + scale_color_manual(
    values = col,
    drop = FALSE
  )

  ggsave(
    sprintf("%s/mixscape_prtb_LDAumap_passgenes_%s.png", OUTPUT_DIR, condition),
    plot | plot1,
    width = 30,
    height = 20,
    limitsize = FALSE
  )

  saveRDS(cond_object, file = sprintf("%s/%s_data.rds", OUTPUT_DIR, condition))
}

# ==================
# Plot Nebulosa plots for the LDA UMAPs
# ==================
dat.LDAnoRT <- readRDS(file.path(OUTPUT_DIR, "noRT_data.rds"))
dat.LDART <- readRDS(file.path(OUTPUT_DIR, "RT_data.rds"))

SCpubr::do_NebulosaPlot(
  sample = dat.LDART,
  features = "sgRNA_logUMI",
  reduction = "ldaumap"
)

SCpubr::do_FeaturePlot(
  sample = dat.LDAnoRT,
  features = "sgRNA_logUMI",
  reduction = "ldaumap"
)

sgRNAs <- unique(dat.LDART@meta.data$sgRNA)

# Loop through each sgRNA for RT
for (sgRNA in sgRNAs) {
  # Create a binary column in meta.data for the current sgRNA
  dat.LDART@meta.data[[paste0("binary_", sgRNA)]] <- ifelse(
    dat.LDART@meta.data$sgRNA == sgRNA,
    1,
    0
  )

  # Generate the Nebulosa plot for the current sgRNA
  p <- SCpubr::do_NebulosaPlot(
    sample = dat.LDART,
    features = paste0("binary_", sgRNA),
    reduction = "ldaumap",
    pt.size = 0.2
  )

  ggsave(
    sprintf("%s/lda_density/GBM43_LDA_RT_nebulosa_%s.png", OUTPUT_DIR, sgRNA),
    p,
    width = 3,
    height = 3.8
  )
}

# Loop over no RT
sgRNAs <- unique(dat.LDAnoRT@meta.data$sgRNA)
for (sgRNA in sgRNAs) {
  # Create a binary column in meta.data for the current sgRNA
  dat.LDAnoRT@meta.data[[paste0("binary_", sgRNA)]] <- ifelse(
    dat.LDAnoRT@meta.data$sgRNA == sgRNA,
    1,
    0
  )

  # Generate the Nebulosa plot for the current sgRNA
  p <- SCpubr::do_NebulosaPlot(
    sample = dat.LDAnoRT,
    features = paste0("binary_", sgRNA),
    reduction = "ldaumap",
    pt.size = 0.2
  )

  ggsave(
    sprintf("%s/lda_density/GBM43_LDA_noRT_nebulosa_%s.png", OUTPUT_DIR, sgRNA),
    p,
    width = 3,
    height = 3.8
  )
}

# Now plot them all together on a page
# Initialize an empty list to store paired plots
paired_plots <- list()

# Loop through each sgRNA and generate paired plots
for (sgRNA in sgRNAs) {
  # Generate Nebulosa plots for RT
  plot_RT <- SCpubr::do_NebulosaPlot(
    sample = dat.LDART,
    features = paste0("binary_", sgRNA),
    reduction = "ldaumap"
  ) +
    ggtitle(paste(sgRNA, "(RT)")) +
    theme(plot.title = element_text(hjust = 0.5, size = 12))

  # Generate Nebulosa plots for noRT
  plot_noRT <- SCpubr::do_NebulosaPlot(
    sample = dat.LDAnoRT,
    features = paste0("binary_", sgRNA),
    reduction = "ldaumap"
  ) +
    ggtitle(paste(sgRNA, "(noRT)")) +
    theme(plot.title = element_text(hjust = 0.5, size = 12))

  # Combine the RT and noRT plots side by side
  paired_plot <- plot_noRT | plot_RT

  # Add the paired plot to the list
  paired_plots[[sgRNA]] <- paired_plot

  ggsave(
    filename = sprintf("%s/lda_density/GBM43_LDA_RT_noRT_nebulosa_%s.png", OUTPUT_DIR, sgRNA),
    plot = paired_plot,
    width = 6,
    height = 4,
    dpi = 300
  )
}

# ==================
# Generate nebulosa plots for the LDA UMAPs
# ==================

dat.LDAnoRT <- readRDS("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/noRT_data.rds")
dat.LDART <- readRDS("/raleighlab/data1/czou/gbm_perturb/gbm_perturb_gbm43_explore_outputs/heterogeneity/gbm_pdx_perturb_mixscape/RT_data.rds")

SCpubr::do_NebulosaPlot(sample = dat.LDART, features = "sgRNA_logUMI", reduction = "ldaumap")
SCpubr::do_FeaturePlot(sample = dat.LDAnoRT, features = "sgRNA_logUMI", reduction = "ldaumap")

sgRNAs <- unique(dat.LDART@meta.data$sgRNA)

# Loop through each sgRNA for RT
for (sgRNA in sgRNAs) {
  # Create a binary column in meta.data for the current sgRNA
  dat.LDART@meta.data[[paste0("binary_", sgRNA)]] <- ifelse(dat.LDART@meta.data$sgRNA == sgRNA, 1, 0)

  # Generate the Nebulosa plot for the current sgRNA
  p <- SCpubr::do_NebulosaPlot(
    sample = dat.LDART,
    features = paste0("binary_", sgRNA),
    reduction = "ldaumap", pt.size = 0.2
  )
  ggsave(paste("lda_density/GBM43_LDA_RT_nebulosa_", sgRNA, ".png", sep = ""), p, width = 3, height = 3.8)
}

# Loop over no RT
sgRNAs <- unique(dat.LDAnoRT@meta.data$sgRNA)
for (sgRNA in sgRNAs) {
  # Create a binary column in meta.data for the current sgRNA
  dat.LDAnoRT@meta.data[[paste0("binary_", sgRNA)]] <- ifelse(dat.LDAnoRT@meta.data$sgRNA == sgRNA, 1, 0)

  # Generate the Nebulosa plot for the current sgRNA
  p <- SCpubr::do_NebulosaPlot(
    sample = dat.LDAnoRT,
    features = paste0("binary_", sgRNA),
    reduction = "ldaumap", pt.size = 0.2
  )
  ggsave(paste("lda_density/GBM43_LDA_noRT_nebulosa_", sgRNA, ".png", sep = ""), p, width = 3, height = 3.8)
}

# Now plot them all together on a page
# Initialize an empty list to store paired plots
paired_plots <- list()

# Loop through each sgRNA and generate paired plots
for (sgRNA in sgRNAs) {
  # Generate Nebulosa plots for RT
  plot_RT <- SCpubr::do_NebulosaPlot(
    sample = dat.LDART,
    features = paste0("binary_", sgRNA),
    reduction = "ldaumap"
  ) +
    ggtitle(paste(sgRNA, "(RT)")) +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) # Center and adjust title size

  # Generate Nebulosa plots for noRT
  plot_noRT <- SCpubr::do_NebulosaPlot(
    sample = dat.LDAnoRT,
    features = paste0("binary_", sgRNA),
    reduction = "ldaumap"
  ) +
    ggtitle(paste(sgRNA, "(noRT)")) +
    theme(plot.title = element_text(hjust = 0.5, size = 12)) # Center and adjust title size

  # Combine the RT and noRT plots side by side
  paired_plot <- plot_noRT | plot_RT

  # Add the paired plot to the list
  paired_plots[[sgRNA]] <- paired_plot

  ggsave(
    filename = paste0("lda_density/sgRNA_", sgRNA, "_RT_vs_noRT.png"),
    plot = paired_plot,
    width = 6,
    height = 4,
    dpi = 300
  )
}


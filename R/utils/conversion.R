# Helper functions to convert Seurat to h5ad objects and vice versa

library(Seurat)
library(SeuratDisk)

gbm43 <- readRDS("/raleighlab/data1/liuj/gbm_perturb/analysis/GBM43_1_malignant_only_annotated_20230817.Rds")

# Process the Seurat object to contain only human genes and use the raw RNA
# counts.
DefaultAssay(gbm43) <- "RNA"
all_genes <- Features(gbm43)
human_genes <- all_genes[grepl("GRCh38-", all_genes)]
gbm43 <- subset(gbm43, features = human_genes)

# Save the Seurat object as an h5ad object.
SaveH5Seurat(gbm43, filename = "data/gbm43.h5Seurat")
Convert("data/gbm43.h5Seurat", dest = "h5ad")

# Export GL261 tumor data

gl261_data <- readRDS("/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_integrated_20230619.Rds")
gl261_data <- subset(gl261_data, scMRMA_manual %in% c("GL261 in vitro", "GL261 in vivo"))
gl261_data <- subset(gl261_data, sorted %in% c("MACSFACS", "invitro"))
gl261_data$sgRNACondSource <- paste(gl261_data$sgRNACond, gl261_data$source, sep = "_")
sgRNACond_counts <- table(gl261_data$sgRNACondSource)
gl261_data <- subset(gl261_data, sgRNACondSource %in%
  names(sgRNACond_counts[sgRNACond_counts >= 5]))
DefaultAssay(gl261_data) <- "RNA"

# Remove the SCT assay for now
gl261_data[["SCT"]] <- NULL

SaveH5Seurat(gl261_data, filename = "data/raw_gl261_sorted_only_tumor_coverage_filtered.h5Seurat")
Convert("data/raw_gl261_sorted_only_tumor_coverage_filtered.h5Seurat", dest = "h5ad")

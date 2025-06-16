# Run DESeq on the GL261 data. Updated 6/13/25.

library(Seurat)
library(org.Hs.eg.db)
library(msigdbr)
library(fgsea)
library(dplyr)

utils.deseq <- new.env()
source("R/utils/deseq.R", local = utils.deseq)

# Use sorted and coverage filters.
gl261_data <- readRDS("/raleighlab/data1/liuj/gbm_perturb/analysis/GL261_integrated_20230619.Rds")
gl261_data <- subset(gl261_data, sgRNA %in% c("non-targeting", "Prkdc"))
gl261_data <- subset(gl261_data, sorted %in% c("MACSFACS", "invitro"))
gl261_data$sgRNACondSource <- paste(gl261_data$sgRNACond, gl261_data$source, sep = "_")
sgRNACond_counts <- table(gl261_data$sgRNACondSource)
gl261_data <- subset(gl261_data, sgRNACondSource %in%
  names(sgRNACond_counts[sgRNACond_counts >= 5]))

# Run DESeq for each context x run type (condNormalized or noRTNormalized)
for (context in c("invitro", "preinf", "CED")) {
  gl261_context <- subset(gl261_data, source == context)
  for (guide in c("Prkdc_RT", "Prkdc_noRT", "non-targeting_RT")) {
    if (guide %in% unique(gl261_context$sgRNACond)) {
      deseq_df <- utils.deseq$find_deseq_differential_genes(
        gl261_context, guide, "non-targeting_noRT", "sgRNACond"
      )
      write.table(deseq_df, file.path("output/gl261_sb28_screens/deseq", sprintf("%s_%s_non-targeting_noRT.csv", context, guide)))
    }
  }
  if ("Prkdc_RT" %in% unique(gl261_context$sgRNACond)) {
    guide <- "Prkdc_RT"
    deseq_df <- utils.deseq$find_deseq_differential_genes(
      gl261_context, guide, "non-targeting_RT", "sgRNACond"
    )
    write.table(deseq_df, file.path("output/gl261_sb28_screens/deseq", sprintf("%s_%s_non-targeting_RT.csv", context, guide)))
  }
}

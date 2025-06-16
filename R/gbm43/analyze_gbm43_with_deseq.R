# DESeq analysis of GBM43 CED in vivo perturb-seq data.
# Updated: 6/13/25.

# =========================================================
# Library imports and constants
# =========================================================

library(Seurat)

utils.deseq <- new.env()
source("R/utils/deseq.R", local = utils.deseq)

output_dir <- "output/gbm43/deseq/noRTNormalized"

# =========================================================
# Import data
# =========================================================

# Import the data
gbm43 <- readRDS("/raleighlab/data1/liuj/gbm_perturb/analysis/GBM43_1_malignant_only_annotated_20230817.Rds")
DefaultAssay(gbm43) <- "RNA"
all_genes <- rownames(gbm43)
human_genes <- all_genes[grep("GRCh38-", all_genes)]
gbm43 <- subset(gbm43, features = human_genes)

# =========================================================
# Run DESeq
# ========================================================= 
all_conditions <- unique(gbm43$sgRNACond)
conditions <- all_conditions[all_conditions != "non-targeting_noRT"]

DefaultAssay(data) <- "RNA"
for (condition in conditions) {
  print(sprintf("Starting deseq analysis for %s", condition))
  deseq_df <- utils.deseq$find_deseq_differential_genes(
    gbm43,
    condition,
    "non-targeting_noRT",
    "sgRNACond",
    seed = 5220
  )
  print(sprintf("Finished deseq analysis for %s", condition))
  write.table(deseq_df, file.path(
    output_dir,
    sprintf("%s_non-targeting_noRT.csv", condition)
  ))
  print("Finished writing table")
}# A series of things to help us prepare for publication

# =========================================================
# Prepare for publication
# =========================================================

# First, take all the DESeq outputs and put them into one sheet
utils.deseq <- new.env()
source("R/utils/deseq.R", local = utils.deseq)

gbm43_noRTNormalized_deseq <- utils.deseq$extract_deseq_list(
  "output/gbm43/deseq/noRTNormalized",
  delim = " ",
  header = TRUE
)

wb <- openxlsx::createWorkbook()
for (name in names(gbm43_noRTNormalized_deseq)) {
  print(name)
  table <- gbm43_noRTNormalized_deseq[[name]]

  # Modify the table so that it's smaller
  table$feature <- gsub("GRCh38-", "", table$feature)
  table <- table %>%
    select(-c(group1, group2, stat, rate1, rate2))
  table <- table %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))

  # Remove any rows where the LFC is NA
  table <- table %>%
    filter(!is.na(log_fc)) %>%
    filter(!is.na(padj))
  sheet_name <- gsub("_non-targeting_noRT", "nt_noRT", name)
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, table)
}
openxlsx::saveWorkbook(wb, "output/manuscript/SuppTable_DESeq.xlsx",
  overwrite = TRUE
)
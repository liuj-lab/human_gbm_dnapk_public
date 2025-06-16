# A series of things to help us prepare for publication

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

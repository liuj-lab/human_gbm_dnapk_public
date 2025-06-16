# Helper functions to run DESeq on a Seurat object

library(DElegate)

find_deseq_differential_genes <- function(data.obj,
                                          group1,
                                          group2,
                                          group_column,
                                          seed = NULL,
                                          replicate_column = NULL) {
  # Takes in a normalized Seurat object with metadata
  # indicating group1 and group2 and a group_column argument
  # indicating the metadata column to be used. Sets a seed, then
  # finds differential expressed genes and outputs a data frame.
  # Note that log fold comparisons will be output as log2(group1/group2)
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Set up futures:
  future::plan(strategy = "future::multicore", workers = 12)

  # Set futures to greater than max capacity; you may need to tweak this
  options(future.globals.maxSize = 8 * 10^9)

  df <- findDE(
    object = data.obj,
    group_column = group_column,
    compare = c(group1, group2),
    method = "deseq",
    replicate_column = replicate_column
  )
  return(df)
}

extract_deseq_list <- function(file_location, delim = ",", header = FALSE) {
  # Given a file location, extract DESeq2 output from all .csv
  # files in that file location and put them in a list
  deseq_list <- list()
  for (file_name in list.files(file_location)) {
    list_name <- gsub(".csv", "", file_name)
    deseq_list[[list_name]] <- read.table(
      sprintf(
        "%s/%s",
        file_location, file_name
      ),
      sep = delim,
      header = header
    )
  }
  return(deseq_list)
}

condense_into_lfc_matrix <- function(list_of_dfs, column = "log_fc") {
  # Given a list of DESeq dataframes, combine them into a matrix
  column_dfs <- lapply(
    names(list_of_dfs), function(name) {
      df <- list_of_dfs[[name]]
      df <- df[, c("feature", column)]
      colnames(df)[colnames(df) == column] <- paste0(column, name)
      return(df)
    }
  )
  column_matrix <- Reduce(function(x, y) {
    merge(x, y, by = "feature", all = TRUE)
  }, column_dfs)
  return(column_matrix)
}

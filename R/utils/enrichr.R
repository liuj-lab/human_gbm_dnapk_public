# Utility functions to interact with the enrichR API. If you would
# like to run this in script form, we also provide a script in the
# scripts directory.

library(enrichR)
library(dplyr)

enrich_directory_of_genes <- function(input_path, output_path, dbs, enrichr_site = "Enrichr") {
  # This function processes a directory of gene lists, queries the Enrichr API for enrichment analysis,
  # and saves the filtered results into specified output directories for each database.
  # Params:
  # input_path: Path to the directory containing gene list files (.txt format)
  # output_path: Path to the directory where the results will be saved
  # dbs: A vector of database names to query in Enrichr
  # enrichr_site: The Enrichr website to use for queries. Defaults to "Enrichr"
  # Returns: None

  gene_files <- list.files(input_path, pattern = "\\.txt$", full.names = TRUE)
  gene_lists <- lapply(gene_files, readLines)

  names(gene_lists) <- sapply(gene_files, function(x) tools::file_path_sans_ext(basename(x)))

  enriched_results <- list()
  for (gene_list_name in names(gene_lists)) {
    genes <- gene_lists[[gene_list_name]]
    enriched_results[[gene_list_name]] <- generate_enrichr_output(genes, dbs)
    Sys.sleep(2)
  }
  filtered_results <- lapply(enriched_results, filter_enrichr_output)

  for (db in dbs) {
    db_output_path <- file.path(output_path, db)
    if (!dir.exists(db_output_path)) {
      dir.create(db_output_path, recursive = TRUE)
    }
    for (gene_list_name in names(filtered_results)) {
      if (!is.null(filtered_results[[gene_list_name]][[db]])) {
        write.csv(filtered_results[[gene_list_name]][[db]],
          file.path(db_output_path, paste0(gene_list_name, "_", db, ".csv")),
          row.names = FALSE
        )
      }
    }
  }
}

generate_enrichr_output <- function(genes, dbs, enrichr_site = "Enrichr") {
  # Take in a list of gene symbols query Enrichr databases for
  # enriched terms. Defaults to Human genes.
  # Params:
  # genes: list of gene symbols
  # dbs: list of enrichr databases
  # enrichr_site: Enrichr website to use. Defaults to "Enrichr"
  # Returns: list of enriched terms

  websiteLive <- getOption("enrichR.live")
  if (websiteLive) {
    listEnrichrSites()
    setEnrichrSite(enrichr_site) # Human genes
  }
  if (websiteLive) available_dbs <- listEnrichrDbs()
  dbs <- available_dbs$libraryName[available_dbs$libraryName %in% dbs]
  enriched_terms <- enrichr(c(genes), dbs)
  return(enriched_terms)
}

filter_enrichr_output <- function(enrichr_dfs) {
  # Takes in a list of dataframes containing Enrichr output and filters each dataframe
  # for terms with an Adjusted P-value less than 0.05.
  # Params:
  # enrichr_dfs: list of dataframes, each representing Enrichr output for a gene list
  # Returns: list of filtered dataframes, with each dataframe containing only the terms
  #          that have an Adjusted P-value < 0.05
  list_dfs <- list()
  for (name in names(enrichr_dfs)) {
    df <- enrichr_dfs[[name]]
    df <- df %>% filter(Adjusted.P.value < 0.05)
    list_dfs[[name]] <- df
  }
  return(list_dfs)
}

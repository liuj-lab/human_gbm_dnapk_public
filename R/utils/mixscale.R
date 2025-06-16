# Functions to run Mixscale analysis on a Seurat object

library(Mixscale)

get_mixscale_scores <- function(data,
                                nt_cell_class,
                                perturb_class,
                                split_by_class = NULL) {
  # Run Mixscale analysis on a Seurat object called data
  # Returns:
  #   A Seurat object with Mixscale scores

  data <- CalcPerturbSig(
    object = data,
    assay = "RNA",
    slot = "data",
    gd.class = perturb_class,
    nt.cell.class = nt_cell_class,
    reduction = "pca",
    ndims = 40,
    num.neighbors = 20,
    new.assay.name = "PRTB",
    split.by = split_by_class
  )
  data <- RunMixscale(
    object = data,
    assay = "PRTB",
    slot = "scale.data",
    labels = perturb_class,
    nt.class.name = nt_cell_class,
    min.de.genes = 5,
    logfc.threshold = 0.2,
    de.assay = "RNA",
    max.de.genes = 100,
    new.class.name = "mixscale_score",
    fine.mode = F,
    verbose = T,
    split.by = split_by_class
  )
  return(data)
}

get_mixscale_degs <- function(data,
                              nt_cell_class,
                              perturb_class,
                              split_by_class = NULL,
                              prtb_list = NULL) {
  # Run Mixscale analysis on a Seurat object called data
  # Returns:
  #   A list of DE genes dataframes
  de_dfs <- Run_wmvRegDE(
    object = data,
    PRTB_list = prtb_list,
    assay = "RNA",
    slot = "counts",
    labels = perturb_class,
    nt.class.name = nt_cell_class,
    logfc.threshold = 0.0,
    split.by = split_by_class
  )
  return(de_dfs)
}

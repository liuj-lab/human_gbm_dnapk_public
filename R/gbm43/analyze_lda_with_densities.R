# Analyze densities out of LDA UMAPs for GBM43 data.
# Original author: David Wu
# This is the latest round of analysis for the density plot
# segmentation. Last updated on 5/16/25

# Note that what we ended up going with was features:
# density_area_difference, z_mean_RT, z_mean_noRT, area_RT, area_noRT

library(here)
library(tidyverse)
library(Seurat)
library(rcartocolor)
library(wutilities)
library(readr)
library(ggplot2)
library(ggrepel)
library(EBImage)
library(ComplexHeatmap)

output_dir <- "output/gbm43/lda/density_analysis_dwu"

# =========================================================
# Import data
# =========================================================
# Load in the RT and noRT data objects
RT_data <- readRDS("output/gbm43/lda/RT_data.rds")
noRT_data <- readRDS("output/gbm43/lda/noRT_data.rds")

# Export the UMAP coordinates
RT_embeddings <- RT_data@reductions$ldaumap@cell.embeddings
noRT_embeddings <- noRT_data@reductions$ldaumap@cell.embeddings

# Join this data with the metadata
RT_metadata <- merge(RT_embeddings, RT_data@meta.data, by = 0)
noRT_metadata <- merge(noRT_embeddings, noRT_data@meta.data, by = 0)

# Do some initial processing
RT_metadata <- RT_metadata %>% rename(UMAP_1 = ldaumap_1, UMAP_2 = ldaumap_2)
noRT_metadata <- noRT_metadata %>% rename(
  UMAP_1 = ldaumap_1,
  UMAP_2 = ldaumap_2
)
metadata <- bind_rows(RT_metadata, noRT_metadata)
rownames(metadata) <- metadata$Row.names

# =========================================================
# Calculate densities
# =========================================================
metadata %>% ggplot(
  aes(x = UMAP_1, y = UMAP_2, color = cond)
) +
  geom_point(alpha = 0.2, size = 0.1) +
  scale_color_carto_d(palette = "Prism")

# Create density plots according to the David Wu style
map(metadata$cond %>% unique(), function(RT_cond) {
  data_subset <- metadata %>% filter(cond == RT_cond)
  map(data_subset$sgRNACond %>% unique(), function(i) {
    print(i)
    if (data_subset %>% filter(sgRNACond == i) %>% nrow() > 5) {
      p <- data_subset %>%
        density_plot(
          group = "sgRNACond",
          target = i,
          background = "all",
          bw = "nrd0",
          point_color = "grey90",
          point_size = 1,
          point_alpha = 0.8
        ) +
        theme(plot.background = element_rect(fill = "white")) +
        ggtitle(i)
      ggsave(
        plot = p,
        filename = paste0(i, ".png"),
        path = here(output_dir, RT_cond),
        dpi = 300,
        w = 4,
        h = 4
      )
    }
  })
})

# Create density tables
density_table <- map(metadata$cond %>% unique(), function(j) {
  seurat_subset <- metadata %>% filter(cond == j)

  map(seurat_subset$sgRNACond %>% unique(), function(i) {
    print2(i)
    if (seurat_subset %>% filter(sgRNACond == i) %>% nrow() > 3) {
      p <- seurat_subset %>%
        density_plot(
          group = "sgRNACond",
          target = i,
          background = "all",
          bw = "nrd0",
          point_color = "grey90",
          point_size = 1,
          point_alpha = 0.8, plot = FALSE, table_only = TRUE
        ) %>%
        mutate(
          sgRNACond = i,
          source = j
        )
    }
  }) %>% bind_rows()
}) %>% bind_rows()

density_table %>%
  select(source, sgRNACond, x, y, z = z_diff) %>%
  write_tsv(here(output_dir, "density_table.tsv.gz"))
max_density <- density_table %>%
  group_by(sgRNACond, source) %>%
  slice_max(z_diff, n = 1)
max_density %>%
  select(source, sgRNACond, max_z = z_diff) %>%
  unique() %>%
  write_tsv(here(output_dir, "max_density_table_2024.tsv.gz"))

# =========================================================
# Preprocessing of density grids
# =========================================================
max_density_grid <- read_tsv(file.path(output_dir, "max_density_table_2024.tsv.gz"))
max_density_grid <- max_density_grid %>%
  mutate(sgRNA = str_split_fixed(sgRNACond, "_", 2)[, 1])
density_grid <- read_tsv(file.path(output_dir, "density_table.tsv.gz"))
density_grid <- density_grid %>%
  mutate(sgRNA = str_split_fixed(sgRNACond, "_", 2)[, 1])

# Remove perturbations that we aren't going to work with
perturbs_to_remove <- c(
  "OVCH1", "GRAMD4", "FBXL16", "NUP205", "CDCA5"
)
max_density_grid <- max_density_grid %>%
  filter(!(sgRNA %in% perturbs_to_remove))
density_grid <- density_grid %>%
  filter(!(sgRNA %in% perturbs_to_remove))

# Normalize the results across samples so that we can at least compare
# in a scale invariant way?
normalized_density_grid <- density_grid %>%
  group_by(sgRNACond) %>%
  mutate(
    z_minmax = (z - min(z)) / (max(z) - min(z)),
    z_zscore = (z - mean(z)) / sd(z)
  ) %>%
  ungroup()

# =========================================================
# Full set of density related features, calculated with normalized density
# =========================================================
feature_df <- normalized_density_grid %>%
  group_by(sgRNA, source) %>%
  mutate(local_threshold = quantile(z_minmax, 0.5, na.rm = TRUE)) %>%
  summarise(
    x_max = x[which.max(z_zscore)],
    y_max = y[which.max(z_zscore)],
    z_max = max(z_zscore),
    z_mean = mean(z_minmax, na.rm = TRUE),
    area = sum(z_minmax > threshold),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = source,
    values_from = c(x_max, y_max, z_max, area, z_mean),
    names_sep = "_"
  ) %>%
  mutate(
    max_density_shift = sqrt(
      (x_max_RT - x_max_noRT)^2 +
        (y_max_RT - y_max_noRT)^2
    ),
    density_ratio = z_max_RT / (z_max_noRT + 1e-6),
    density_area_difference = area_noRT - area_RT,
    density_area_RT = area_RT,
    density_area_noRT = area_noRT
  )

# Now run K-means on the features
features <- feature_df %>%
  dplyr::select(c(
    density_area_difference,
    z_mean_RT,
    z_mean_noRT,
    area_RT,
    area_noRT
  ))
scaled_features <- scale(features)

# Run PCA
pca_result <- prcomp(scaled_features, scale. = FALSE)

# Extract first two PCs
pca_df <- as.data.frame(pca_result$x[, 1:2])
pca_df$sgRNA <- feature_df$sgRNA
pca_df$cluster <- as.factor(feature_df$cluster)

# Plot
plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 3) +
  geom_text_repel(aes(label = sgRNA), size = 3) +
  labs(title = "PCA of sgRNA Density Features") +
  theme_minimal()
ggsave(file.path(output_dir, "3Means_PCA_sgRNA_normalizedDensity_features.png"),
  plot = plot,
  height = 7, width = 7
)

# Out of curiosity, what are the PC loadings?
loadings <- pca_result$rotation
loading_df <- as.data.frame(loadings)
loading_df$feature <- rownames(loading_df)
loading_long <- loading_df %>%
  pivot_longer(
    cols = starts_with("PC"),
    names_to = "PC",
    values_to = "loading"
  )

top_pc1 <- loading_long %>%
  filter(PC == "PC1") %>%
  arrange(desc(abs(loading))) %>%
  slice(1:10)
print(top_pc1)
top_pc2 <- loading_long %>%
  filter(PC == "PC2") %>%
  arrange(desc(abs(loading))) %>%
  slice(1:10)
print(top_pc2)

# Try simply running hierarchical clustering on the features
rownames(scaled_features) <- feature_df$sgRNA
clustering <- hclust(dist(scaled_features), method = "average")
plot(clustering)

# Hey this seems to work!
# Cut the tree at 3 clusters and get a sense of what the clusters are
clusterCut <- cutree(clustering, 3)
clusters <- split(names(clusterCut), clusterCut)

RT_INDEPENDENT_DISTRIBUTED_PERTURBATIONS <- c(
  "BORA", "C20orf196", "CAMKK2", "ERCC6L2",
  "RNF8", "ATRIP", "CCT6A", "CCT7",
  "DNAJC17", "MMS22L", "TFIP11", "LRR1"
)
RT_DEPENDENT_PERTURBATIONS <- c(
  "CENPT", "PRKDC", "RIF1"
)
RT_INDEPENDENT_CLUSTERED_PERTURBATIONS <- c(
  "MED17", "SFPQ", "CENPJ", "DHX16", "DDOST",
  "EIF2B2", "GKAP1", "POLD3", "POLR2C", "PSMC5",
  "RBX1", "RPLP0", "RBM42", "RRN3", "SLC39A7"
)

# Plot a heatmap of the scaled feature matrix
ht <- Heatmap(
  as.matrix(scaled_features),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  row_dend_width = unit(6, "cm"),
  width = unit(10, "cm"),
  height = unit(12, "cm"),
  clustering_method_rows = "average"
)
pdf(file.path(output_dir, "feature_heatmap_with_average_linkage.pdf"),
  height = 10, width = 12
)
draw(ht)
dev.off()
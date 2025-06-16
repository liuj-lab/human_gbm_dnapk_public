# Utility script for plotting DSPIN data in R friendly ways

library(ComplexHeatmap)
library(colorRamp2)

OUTPUT_DIR <- "output/gbm43/dspin/combined_oNMF_outputs/dual_model_analysis_RTnoRTNorm_noRT"

# Construct a table that categorizes the perturbations
RT_INDEPENDENT_DISTRIBUTED_PERTURBATIONS <- c(
  "LRR1", "RNF8", "MMS22L", "C20orf196", "TFIP11", "ATRIP", "DNAJC17",
  "BORA", "ERCC6L2", "GKAP1", "CAMKK2", "CCT6A", "CCT7", "DDOST"
)
RT_DEPENDENT_PERTURBATIONS <- c("RIF1", "PRKDC", "CENPT")
RT_INDEPENDENT_CLUSTERED_PERTURBATIONS <- c(
  "RBM42", "DHX16", "POLD3", "SLC39A7", "CENPJ", "EIF2B2", "RRN3", "PSMC5",
  "RBX1", "MED17", "POLR2C", "SFPQ", "RPLP0"
)
NT_PERTURBATIONS <- c("non.targeting")
perturbation_class_table <- data.frame(
  Perturbations = c(
    RT_DEPENDENT_PERTURBATIONS,
    RT_INDEPENDENT_DISTRIBUTED_PERTURBATIONS,
    RT_INDEPENDENT_CLUSTERED_PERTURBATIONS,
    NT_PERTURBATIONS
  ),
  Class = c(
    rep("RT_DEP", length(RT_DEPENDENT_PERTURBATIONS)),
    rep("RT_IND_DIST", length(RT_INDEPENDENT_DISTRIBUTED_PERTURBATIONS)),
    rep("RT_IND_CLUSTER", length(RT_INDEPENDENT_CLUSTERED_PERTURBATIONS)),
    rep("NT", length(NT_PERTURBATIONS))
  )
)
perturbation_lookup <- setNames(
  perturbation_class_table$Class,
  perturbation_class_table$Perturbation
)

# ======================================
# Plot the changes in edges
# ======================================
# Read in the data
diff_matrix <- read.csv(
  file.path(OUTPUT_DIR, "diff_matrix.csv"),
  header = TRUE,
  row.names = 1
)
sign_change_matrix <- read.csv(
  file.path(OUTPUT_DIR, "sign_change_matrix.csv"),
  header = TRUE,
  row.names = 1
)

# Plot as a complex heatmap with splits for perturbations
mat <- as.matrix(diff_matrix)
heatmap <- Heatmap(mat,
  name = "RTnoRTNorm - noRT",
  column_split = perturbation_lookup[colnames(mat)]
)
pdf(file.path(OUTPUT_DIR, "diff_matrix_heatmap_split.pdf"),
  width = 20, height = 10
)
draw(heatmap, heatmap_legend_side = "top", padding = unit(c(0, 2, 0, 8), "cm"))
dev.off()

# Now try ways to incorporate sign change information
# First, let's try using numbers inside the heatmap to show magnitude and
# colors to show the type of change.
diff_mat <- as.matrix(diff_matrix)
type_mat <- as.matrix(sign_change_matrix)
col_fun <- function(sign_value) {
  ifelse(sign_value == 1, "pink",
    ifelse(sign_value == -1, "lightblue", "lavender")
  )
}
heatmap <- Heatmap(
  diff_mat,
  name = "RTnoRTNorm - noRT",
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(fill = col_fun(type_mat[i, j]), col = NA))
    grid.text(sprintf("%.1f", diff_mat[i, j]), x, y, gp = gpar(fontsize = 10))
  },
  column_split = perturbation_lookup[colnames(type_mat)],
)
pdf(file.path(OUTPUT_DIR, "diff_matrix_heatmap_split_colors_type_numbers_diff.pdf"), width = 20, height = 10)
draw(heatmap, heatmap_legend_side = "top", padding = unit(c(0, 2, 0, 8), "cm"))
dev.off()

# Next, let's try using symbols in the heatmap to delineate the type of change
diff_mat <- as.matrix(diff_matrix)
type_mat <- as.matrix(sign_change_matrix)
heatmap <- Heatmap(
  diff_mat,
  name = "RTnoRTNorm - noRT",
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (type_mat[i, j] == 1) {
      grid.points(x, y, pch = 3, gp = gpar(col = "black", cex = 0.8))
    } else if (type_mat[i, j] == -1) {
      grid.points(x, y, pch = 45, gp = gpar(col = "black", cex = 1))
    } else {
      grid.points(x, y, pch = 5, gp = gpar(col = "black", cex = 1))
    }
  },
  column_split = perturbation_lookup[colnames(type_mat)],
  width = unit(30, "cm"),
  height = unit(15, "cm"),
)
legends <- list(
  Legend(labels = c("Both positive", "Both negative", "Sign flip"), type = "points", pch = c(3, 45, 5))
)
pdf(file.path(OUTPUT_DIR, "diff_matrix_heatmap_split_symbols_type_numbers_diff.pdf"), width = 20, height = 10)
draw(heatmap,
  heatmap_legend_side = "top", padding = unit(c(0, 2, 0, 8), "cm"),
  heatmap_legend_list = legends
)
dev.off()

# What if we try doing this using a 3D visualization, where the height represents
# magnitude of change and color represents type of change? We'll have to separate
# out the negative and positive signal.
# Positive signal
diff_mat <- as.matrix(diff_matrix)
diff_mat[diff_mat < 0] <- 0
type_mat <- as.matrix(sign_change_matrix)
col_fun <- function(sign_value) {
  ifelse(sign_value == 1, "pink",
    ifelse(sign_value == -1, "lightblue", "lavender")
  )
}
bar_angle <- 60
bar_rel_width <- 0.6
bar_rel_height <- 0.6
bar_max_length <- unit(1, "cm")
heatmap <- Heatmap(
  diff_mat,
  name = "RTnoRTNorm - noRT",
  col = colorRamp2(c(0, 2), c("white", "red")),
  layer_fun = function(j, i, x, y, w, h, f) {
    v <- pindex(diff_mat, i, j)
    od <- rank(order(-as.numeric(y), -as.numeric(x)))
    grid.rect(x[od], y[od], w[od], h[od], gp = gpar(col = "white", fill = col_fun(type_mat[i, j])))
    bar3D(x[od], y[od], w[od] * bar_rel_width, h[od] * bar_rel_height, v[od] / max(diff_mat) * bar_max_length, theta = bar_angle, fill = f[od])
  },
  column_split = perturbation_lookup[colnames(type_mat)],
)
pdf(file.path(OUTPUT_DIR, "diff_matrix_heatmap_split_color_type_height_diff_pos.pdf"), width = 20, height = 10)
draw(heatmap, heatmap_legend_side = "top", padding = unit(c(0, 2, 0, 8), "cm"))
dev.off()

# ======================================
# Plot the perturbation vectors in relation to one another
# ======================================
# Import the data
combined_rel_h <- read.csv(
  file.path(OUTPUT_DIR, "combined_rel_h_df.csv"),
  header = TRUE,
  row.names = 1
)

# Now let's plot this data, split by the perturbations
mat <- as.matrix(combined_rel_h)
perturbation_lookup_terms <- sapply(
  strsplit(colnames(mat), split = "_"),
  `[`, 1
)
perturbation_split <- factor(perturbation_lookup[perturbation_lookup_terms],
  levels = c("RT_DEP", "RT_IND_DIST", "RT_IND_CLUSTER", "NT")
)
perturb_order <- unlist(
  lapply(names(perturbation_lookup), function(x) paste0(x, c("_RT", "_noRT")))
)
perturb_order <- perturb_order[1:length(perturb_order) - 1]
col_fun <- colorRamp2(c(-1, 0, 1), c("#3285CC", "white", "#E84B23"))
heatmap <- Heatmap(mat,
  name = "Perturbation vectors", cluster_rows = TRUE, cluster_columns = FALSE,
  column_split = perturbation_split, col = col_fun,
  column_order = perturb_order,
  width = unit(30, "cm"),
  height = unit(12, "cm"),
)
pdf(file.path(OUTPUT_DIR, "rel_h_matrix_heatmap.pdf"), width = 20, height = 10)
ht <- draw(heatmap, heatmap_legend_side = "top", padding = unit(c(0, 2, 0, 8), "cm"))
dev.off()

# Extract the matrix
heatmap_matrix <- heatmap@matrix
r_order <- row_order(ht)
c_order <- column_order(ht)
c_order <- unlist(c_order)
heatmap_matrix <- heatmap_matrix[r_order, c_order]
write.table(heatmap_matrix, "output/manuscript/Fig1E_dspin_rel_h_matrix_dataframe.csv")

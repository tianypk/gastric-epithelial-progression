project_root <- "/Users/huahaijiao/Desktop/02_study/0206_artical/Hp_gastric_cancer"
result_root <- file.path(
  project_root,
  "gastric_cancer", "results", "GSE134520", "epithelial", "pseudotime_main"
)
scripts_dir <- file.path(result_root, "scripts")
tables_dir <- file.path(result_root, "tables")
figures_dir <- file.path(result_root, "figures")
manuscript_dir <- file.path(result_root, "manuscript")

dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(data.table)
  library(igraph)
  library(BiocNeighbors)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(svglite)
})

normalize_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

mean_expr_score <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(mat)))
  }
  Matrix::colMeans(mat[genes, , drop = FALSE])
}

scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

compute_cluster_centers <- function(df, x = "UMAP1", y = "UMAP2", cluster = "epi_cluster_id") {
  as.data.table(df)[, .(
    center_x = median(get(x), na.rm = TRUE),
    center_y = median(get(y), na.rm = TRUE),
    n_cells = .N
  ), by = cluster]
}

build_knn_mst <- function(coords, k = 12) {
  nn <- BiocNeighbors::findKNN(coords, k = k)
  edge_list <- data.table(
    from = rep(seq_len(nrow(coords)), k),
    to = as.vector(nn$index),
    weight = as.vector(nn$distance)
  )
  edge_list <- edge_list[from != to]
  edge_list[, key := ifelse(from < to, paste(from, to, sep = "_"), paste(to, from, sep = "_"))]
  edge_list <- edge_list[, .(from = first(from), to = first(to), weight = min(weight)), by = key]
  g <- graph_from_data_frame(edge_list[, .(from, to, weight)], directed = FALSE, vertices = data.frame(name = seq_len(nrow(coords))))
  mst(g, weights = E(g)$weight)
}

compute_rooted_pseudotime <- function(coords, clusters, root_cluster = "E18", k = 12) {
  mst_graph <- build_knn_mst(coords, k = k)
  root_candidates <- which(clusters == root_cluster)
  if (length(root_candidates) == 0) {
    stop("No cells found for root cluster: ", root_cluster)
  }
  root_center <- colMeans(coords[root_candidates, , drop = FALSE])
  root_mat <- matrix(
    root_center,
    nrow = length(root_candidates),
    ncol = ncol(coords),
    byrow = TRUE
  )
  root_id <- root_candidates[which.min(rowSums((coords[root_candidates, , drop = FALSE] - root_mat)^2))]
  dist_vec <- distances(mst_graph, v = as.character(root_id), to = V(mst_graph), weights = E(mst_graph)$weight)
  pseudotime <- as.numeric(dist_vec[1, ])
  pseudotime
}

sce_path <- file.path(
  project_root,
  "gastric_cancer", "results", "GSE134520", "epithelial", "rds",
  "gse134520_epithelial_sce_reclustered_annotated.rds"
)
gain_union_path <- file.path(
  project_root,
  "gastric_cancer", "results", "Figure1_integrated_candidate_selection_figure", "r_repro_main47", "tables",
  "01_gain_union_main.csv"
)
loss_union_path <- file.path(
  project_root,
  "gastric_cancer", "results", "Figure1_integrated_candidate_selection_figure", "r_repro_main47", "tables",
  "01_loss_union_main.csv"
)

required_files <- c(sce_path, gain_union_path, loss_union_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

sce <- readRDS(sce_path)
meta <- as.data.frame(SummarizedExperiment::colData(sce))
progression_clusters <- c("E10", "E11", "E13", "E17", "E18", "E19")
keep <- meta$epi_cluster_id %in% progression_clusters
sce_prog <- sce[, keep]
meta_prog <- as.data.table(as.data.frame(SummarizedExperiment::colData(sce_prog)))

umap <- reducedDim(sce_prog, "UMAP")
colnames(umap) <- c("UMAP1", "UMAP2")
pca <- reducedDim(sce_prog, "PCA")[, seq_len(min(10, ncol(reducedDim(sce_prog, "PCA")))), drop = FALSE]

gain_union <- fread(gain_union_path)
loss_union <- fread(loss_union_path)
gain_genes <- normalize_gene(gain_union$gene)
loss_genes <- normalize_gene(loss_union$gene)

log_mat <- logcounts(sce_prog)
rownames(log_mat) <- normalize_gene(rownames(log_mat))

gain_score <- mean_expr_score(log_mat, gain_genes)
loss_score <- mean_expr_score(log_mat, loss_genes)
identity_loss_score <- -loss_score

pseudotime_raw <- compute_rooted_pseudotime(pca, clusters = meta_prog$epi_cluster_id, root_cluster = "E18", k = 12)
pseudotime_scaled <- scale01(pseudotime_raw)

meta_prog[, UMAP1 := umap[, 1]]
meta_prog[, UMAP2 := umap[, 2]]
meta_prog[, pseudotime_raw := pseudotime_raw]
meta_prog[, pseudotime := pseudotime_scaled]
meta_prog[, gain_union_score := gain_score]
meta_prog[, loss_union_score := loss_score]
meta_prog[, identity_loss_score := identity_loss_score]
meta_prog[, lesion_group := factor(lesion_group, levels = c("NAG", "CAG", "IM", "EGC"))]
meta_prog[, epi_cluster_id := factor(epi_cluster_id, levels = progression_clusters)]

cluster_centers <- compute_cluster_centers(meta_prog, cluster = "epi_cluster_id")
setnames(cluster_centers, old = names(cluster_centers)[1], new = "epi_cluster_id")
cluster_pt <- meta_prog[, .(
  pseudotime_median = median(pseudotime, na.rm = TRUE),
  pseudotime_mean = mean(pseudotime, na.rm = TRUE),
  gain_union_score = median(gain_union_score, na.rm = TRUE),
  identity_loss_score = median(identity_loss_score, na.rm = TRUE),
  lesion_group_major = names(sort(table(lesion_group), decreasing = TRUE))[1]
), by = epi_cluster_id]
cluster_centers <- merge(cluster_centers, cluster_pt, by = "epi_cluster_id", all.x = TRUE)
setorder(cluster_centers, pseudotime_median)
cluster_centers[, next_cluster := shift(epi_cluster_id, type = "lead")]
cluster_edges <- cluster_centers[!is.na(next_cluster), .(
  from = epi_cluster_id,
  to = next_cluster
)]
cluster_edges <- merge(cluster_edges, cluster_centers[, .(epi_cluster_id, from_x = center_x, from_y = center_y)], by.x = "from", by.y = "epi_cluster_id", all.x = TRUE)
cluster_edges <- merge(cluster_edges, cluster_centers[, .(epi_cluster_id, to_x = center_x, to_y = center_y)], by.x = "to", by.y = "epi_cluster_id", all.x = TRUE)

stage_summary <- meta_prog[, .(
  n_cells = .N,
  pseudotime_median = median(pseudotime, na.rm = TRUE),
  gain_union_score = median(gain_union_score, na.rm = TRUE),
  identity_loss_score = median(identity_loss_score, na.rm = TRUE)
), by = lesion_group][order(lesion_group)]

cluster_stage_summary <- meta_prog[, .(
  n_cells = .N,
  proportion = .N / nrow(meta_prog)
), by = .(lesion_group, epi_cluster_id)]

fwrite(meta_prog, file.path(tables_dir, "gse134520_progression_cells_pseudotime_metadata.csv"))
fwrite(cluster_centers, file.path(tables_dir, "gse134520_progression_cluster_pseudotime_summary.csv"))
fwrite(cluster_edges, file.path(tables_dir, "gse134520_progression_cluster_pseudotime_edges.csv"))
fwrite(stage_summary, file.path(tables_dir, "gse134520_progression_stage_pseudotime_summary.csv"))
fwrite(cluster_stage_summary, file.path(tables_dir, "gse134520_progression_cluster_stage_summary.csv"))

cluster_colors <- c(
  E10 = "#D98C2B",
  E11 = "#C85D48",
  E13 = "#9E2F3F",
  E17 = "#5A9A73",
  E18 = "#4688B8",
  E19 = "#3A5F8A"
)
group_colors <- c(
  NAG = "#8EA8C7",
  CAG = "#C5D6EA",
  IM = "#E7C58A",
  EGC = "#D67A74"
)

p_umap_cluster <- ggplot(meta_prog, aes(UMAP1, UMAP2)) +
  geom_point(aes(color = epi_cluster_id), size = 0.35, alpha = 0.65) +
  geom_segment(
    data = cluster_edges,
    aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.8,
    arrow = grid::arrow(length = unit(0.12, "inches"), type = "closed")
  ) +
  geom_label(
    data = cluster_centers,
    aes(x = center_x, y = center_y, label = epi_cluster_id),
    inherit.aes = FALSE,
    size = 4,
    fill = "white",
    label.size = 0.2
  ) +
  scale_color_manual(values = cluster_colors) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(title = "Progression-related epithelial states and inferred pseudotime backbone", color = "Cluster")

p_umap_pseudotime <- ggplot(meta_prog, aes(UMAP1, UMAP2)) +
  geom_point(aes(color = pseudotime), size = 0.35, alpha = 0.7) +
  geom_segment(
    data = cluster_edges,
    aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
    inherit.aes = FALSE,
    color = "white",
    linewidth = 0.7,
    arrow = grid::arrow(length = unit(0.10, "inches"), type = "closed")
  ) +
  scale_color_viridis_c(option = "magma", end = 0.95) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(title = "Inferred epithelial pseudotime", color = "Pseudotime")

p_score_gain <- ggplot(meta_prog, aes(pseudotime, gain_union_score, color = lesion_group)) +
  geom_point(size = 0.25, alpha = 0.22) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.1) +
  scale_color_manual(values = group_colors) +
  theme_bw(base_size = 14) +
  labs(title = "Gain-union score increases along epithelial pseudotime", x = "Scaled pseudotime", y = "Gain-union score", color = "Group")

p_score_loss <- ggplot(meta_prog, aes(pseudotime, identity_loss_score, color = lesion_group)) +
  geom_point(size = 0.25, alpha = 0.22) +
  geom_smooth(method = "loess", se = TRUE, linewidth = 1.1) +
  scale_color_manual(values = group_colors) +
  theme_bw(base_size = 14) +
  labs(title = "Identity-loss score increases along epithelial pseudotime", x = "Scaled pseudotime", y = "Identity-loss score", color = "Group")

p_stage_box <- ggplot(meta_prog, aes(lesion_group, pseudotime, fill = lesion_group)) +
  geom_boxplot(outlier.size = 0.15, width = 0.72) +
  scale_fill_manual(values = group_colors) +
  theme_bw(base_size = 14) +
  labs(title = "Pseudotime shifts across histologic stages", x = NULL, y = "Scaled pseudotime") +
  theme(legend.position = "none")

p_cluster_order <- ggplot(cluster_centers, aes(reorder(epi_cluster_id, pseudotime_median), pseudotime_median, fill = epi_cluster_id)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = cluster_colors) +
  theme_bw(base_size = 14) +
  labs(title = "Cluster-level ordering along pseudotime", x = NULL, y = "Median pseudotime") +
  theme(legend.position = "none")

main_fig <- (p_umap_cluster + p_umap_pseudotime) / (p_score_gain + p_score_loss) / (p_stage_box + p_cluster_order)

png(file.path(figures_dir, "gse134520_progression_pseudotime_gain_loss_main.png"), width = 3200, height = 4200, res = 320)
print(main_fig)
dev.off()

pdf(file.path(figures_dir, "gse134520_progression_pseudotime_gain_loss_main.pdf"), width = 12, height = 16)
print(main_fig)
dev.off()

svglite::svglite(file.path(figures_dir, "gse134520_progression_pseudotime_gain_loss_main.svg"), width = 12, height = 16)
print(main_fig)
dev.off()

tiff(file.path(figures_dir, "gse134520_progression_pseudotime_gain_loss_main.tiff"), width = 3200, height = 4200, res = 320, compression = "lzw")
print(main_fig)
dev.off()

writeLines(
  c(
    "# Epithelial progression pseudotime and gain/loss dynamics",
    "",
    "## Design",
    "- Cells were restricted to the six progression-related epithelial states: E10, E11, E13, E17, E18, and E19.",
    "- A graph-based pseudotime was inferred from the progression-cell PCA space using a k-nearest-neighbor minimum spanning tree.",
    "- The root state was anchored at E18 because this cluster is IM-dominant and represents an intestinalized transition state that is more proximal to early metaplastic progression than E13/E11.",
    "- Gain-union score was calculated as the average log-expression of the 468 single-cell gain-union genes.",
    "- Identity-loss score was defined as the negative mean expression of the 194 loss-union genes, such that higher values indicate stronger erosion of gastric identity.",
    "",
    "## Main observations",
    paste0("- Progression cells included: ", ncol(sce_prog)),
    paste0("- Root cluster: E18"),
    paste0("- Terminal high-pseudotime states were dominated by E13 and E11, whereas E17/E18/E19 and part of E10 occupied earlier to intermediate regions."),
    "- Both gain-union score and identity-loss score increased with pseudotime, supporting a coupled process of acquired malignant-like programming and gastric identity erosion."
  ),
  file.path(manuscript_dir, "gse134520_progression_pseudotime_gain_loss_note.md")
)

cat("Pseudotime analysis completed.\n")
cat("Progression cells:", ncol(sce_prog), "\n")
cat("Saved main figure to:", file.path(figures_dir, "gse134520_progression_pseudotime_gain_loss_main.pdf"), "\n")

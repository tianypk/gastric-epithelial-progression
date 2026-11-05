project_root <- "/Users/huahaijiao/Desktop/02_study/0206_artical/Hp_gastric_cancer"
result_root <- file.path(
  project_root,
  "gastric_cancer", "results", "Figure5_regulon_pathway_bulk"
)
pyscenic_run_tag <- "pyscenic_out_cluster_aware_server_0121_20261001"
pyscenic_prefix <- "gastric_progression_cluster_aware_0121"
pyscenic_dir <- file.path(result_root, pyscenic_run_tag)
tables_dir <- file.path(result_root, "tables")
figures_dir <- file.path(result_root, "figures")
rds_dir <- file.path(result_root, "rds")
logs_dir <- file.path(result_root, "logs")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

auc_path <- file.path(pyscenic_dir, paste0(pyscenic_prefix, ".auc_mtx.csv"))
cell_meta_path <- file.path(tables_dir, "01_pyscenic_cell_metadata_cluster_aware.csv")
input_rds_path <- file.path(rds_dir, "01b_progression_epithelial_for_regulon_cluster_aware.rds")

required_files <- c(auc_path, cell_meta_path, input_rds_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

cluster_order <- c("E10", "E11", "E13", "E17", "E18", "E19")
top_n_per_cluster <- 6L
min_selected_regulons <- 24L

zscore_rows <- function(mat) {
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds <- apply(mat, 1, sd, na.rm = TRUE)
  row_sds[row_sds == 0 | is.na(row_sds)] <- 1
  z <- sweep(mat, 1, row_means, FUN = "-")
  z <- sweep(z, 1, row_sds, FUN = "/")
  z[z > 2] <- 2
  z[z < -2] <- -2
  z
}

add_regulon_metrics <- function(activity_dt, mean_mat, cluster_levels) {
  out <- copy(activity_dt)
  out[, mean_diff_vs_others := {
    cur_cluster <- regulon_group[1]
    if (length(cluster_levels) == 1) {
      mean_auc
    } else {
      mean_auc - rowMeans(mean_mat[regulon, setdiff(cluster_levels, cur_cluster), drop = FALSE], na.rm = TRUE)
    }
  }, by = regulon_group]
  out[, mean_range_across_clusters := {
    vals <- mean_mat[regulon, cluster_levels, drop = FALSE]
    apply(vals, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }, by = regulon_group]
  out[, cluster_rank_by_diff := frank(-mean_diff_vs_others, ties.method = "first"), by = regulon_group]
  out
}

compute_regulon_correlations <- function(auc_mat, score_vec, metric_label) {
  valid_score <- is.finite(score_vec)
  if (sum(valid_score) < 10) {
    stop("Fewer than 10 cells have non-missing values for ", metric_label)
  }

  corr_dt <- rbindlist(lapply(
    seq_len(nrow(auc_mat)),
    function(i) {
      auc_values <- auc_mat[i, ]
      ok <- is.finite(auc_values) & valid_score
      n_ok <- sum(ok)
      if (n_ok < 10) {
        return(data.table(
          regulon = rownames(auc_mat)[i],
          n_cells = n_ok,
          spearman_rho = NA_real_,
          p_value = NA_real_
        ))
      }

      auc_ok <- auc_values[ok]
      score_ok <- score_vec[ok]
      if (sd(auc_ok) == 0 || sd(score_ok) == 0) {
        return(data.table(
          regulon = rownames(auc_mat)[i],
          n_cells = n_ok,
          spearman_rho = NA_real_,
          p_value = NA_real_
        ))
      }

      test <- suppressWarnings(cor.test(
        auc_ok,
        score_ok,
        method = "spearman",
        exact = FALSE
      ))

      data.table(
        regulon = rownames(auc_mat)[i],
        n_cells = n_ok,
        spearman_rho = unname(test$estimate),
        p_value = test$p.value
      )
    }
  ))

  corr_dt[, fdr := p.adjust(p_value, method = "BH")]
  setnames(
    corr_dt,
    c("n_cells", "spearman_rho", "p_value", "fdr"),
    paste0(c("n_cells_", "rho_", "p_", "fdr_"), metric_label)
  )
  corr_dt
}

cat("Reading AUCell matrix...\n")
auc_dt <- fread(auc_path)
if (ncol(auc_dt) < 3 || names(auc_dt)[1] != "Regulon") {
  stop("Unexpected AUCell matrix format: first column should be 'Regulon'")
}

regulon_names <- auc_dt[["Regulon"]]
cell_ids_in_auc <- setdiff(names(auc_dt), "Regulon")

cat("Reading cell metadata...\n")
cell_meta <- fread(cell_meta_path)
required_meta_cols <- c("cell_id", "regulon_group")
if (!all(required_meta_cols %in% names(cell_meta))) {
  stop("Cell metadata is missing required columns: ", paste(setdiff(required_meta_cols, names(cell_meta)), collapse = ", "))
}

cell_meta <- cell_meta[regulon_group %in% cluster_order]
cell_meta[, regulon_group := factor(as.character(regulon_group), levels = cluster_order)]

common_cells <- intersect(cell_ids_in_auc, cell_meta$cell_id)
if (length(common_cells) == 0) {
  stop("No overlapping cells between AUCell matrix and cell metadata")
}

cell_meta <- cell_meta[match(common_cells, cell_id)]
if (anyNA(cell_meta$cell_id)) {
  stop("Failed to align cell metadata to AUCell matrix columns")
}

present_clusters <- cluster_order[cluster_order %in% unique(as.character(cell_meta$regulon_group))]
if (length(present_clusters) < 2) {
  stop("Fewer than two progression clusters are available in cell metadata")
}

auc_mat <- as.matrix(auc_dt[, ..common_cells])
mode(auc_mat) <- "numeric"
rownames(auc_mat) <- regulon_names

cluster_cell_counts <- cell_meta[, .N, by = regulon_group]
setorder(cluster_cell_counts, regulon_group)
fwrite(cluster_cell_counts, file.path(tables_dir, "02_regulon_cluster_cell_counts.csv"))

cat("Summarizing regulon activity by cluster...\n")
mean_list <- list()
median_list <- list()
pct_active_list <- list()

for (cl in present_clusters) {
  idx <- cell_meta$regulon_group == cl
  cl_mat <- auc_mat[, idx, drop = FALSE]
  mean_list[[cl]] <- rowMeans(cl_mat, na.rm = TRUE)
  median_list[[cl]] <- apply(cl_mat, 1, median, na.rm = TRUE)
  pct_active_list[[cl]] <- rowMeans(cl_mat > 0, na.rm = TRUE)
}

cluster_mean_mat <- do.call(cbind, mean_list)
cluster_median_mat <- do.call(cbind, median_list)
cluster_pct_active_mat <- do.call(cbind, pct_active_list)
rownames(cluster_mean_mat) <- regulon_names
rownames(cluster_median_mat) <- regulon_names
rownames(cluster_pct_active_mat) <- regulon_names

activity_dt <- rbindlist(lapply(
  present_clusters,
  function(cl) {
    data.table(
      regulon = regulon_names,
      regulon_group = cl,
      n_cells = sum(cell_meta$regulon_group == cl),
      mean_auc = cluster_mean_mat[, cl],
      median_auc = cluster_median_mat[, cl],
      pct_active = cluster_pct_active_mat[, cl]
    )
  }
))

activity_dt <- add_regulon_metrics(activity_dt, cluster_mean_mat, present_clusters)

dominant_cluster_idx <- max.col(cluster_mean_mat[, present_clusters, drop = FALSE], ties.method = "first")
dominant_cluster_map <- data.table(
  regulon = regulon_names,
  dominant_cluster = present_clusters[dominant_cluster_idx],
  dominant_mean_auc = apply(cluster_mean_mat[, present_clusters, drop = FALSE], 1, max, na.rm = TRUE)
)

activity_dt <- merge(activity_dt, dominant_cluster_map, by = "regulon", all.x = TRUE)
setorder(activity_dt, regulon_group, -mean_diff_vs_others, -mean_auc)

fwrite(activity_dt, file.path(tables_dir, "02_regulon_activity_by_cluster.csv"))
fwrite(
  as.data.table(cluster_mean_mat, keep.rownames = "regulon"),
  file.path(tables_dir, "02_regulon_activity_by_cluster_mean_wide.csv")
)
fwrite(
  as.data.table(cluster_median_mat, keep.rownames = "regulon"),
  file.path(tables_dir, "02_regulon_activity_by_cluster_median_wide.csv")
)
fwrite(
  as.data.table(cluster_pct_active_mat, keep.rownames = "regulon"),
  file.path(tables_dir, "02_regulon_activity_by_cluster_pct_active_wide.csv")
)

selected_candidates <- rbindlist(lapply(
  present_clusters,
  function(cl) {
    sub <- activity_dt[regulon_group == cl][order(-mean_diff_vs_others, -mean_auc, -pct_active)]
    head(sub, top_n_per_cluster)
  }
))

selected_dt <- selected_candidates[order(-mean_diff_vs_others, -mean_auc)][!duplicated(regulon)]

if (nrow(selected_dt) < min_selected_regulons) {
  variable_dt <- unique(activity_dt[, .(regulon, dominant_cluster, dominant_mean_auc, mean_range_across_clusters)])
  setorder(variable_dt, -mean_range_across_clusters, -dominant_mean_auc)
  variable_dt <- variable_dt[!regulon %in% selected_dt$regulon]
  need_n <- min_selected_regulons - nrow(selected_dt)
  if (need_n > 0) {
    extra_regulons <- head(variable_dt, need_n)$regulon
    extra_dt <- activity_dt[regulon %in% extra_regulons][order(-mean_range_across_clusters, -dominant_mean_auc)]
    extra_dt <- extra_dt[!duplicated(regulon)]
    selected_dt <- rbind(selected_dt, extra_dt, fill = TRUE)
  }
}

selected_dt <- unique(selected_dt[, .(
  regulon,
  dominant_cluster,
  dominant_mean_auc,
  mean_range_across_clusters,
  representative_cluster = regulon_group,
  representative_mean_auc = mean_auc,
  representative_mean_diff_vs_others = mean_diff_vs_others,
  representative_pct_active = pct_active
)])

selected_dt[, dominant_cluster := factor(dominant_cluster, levels = present_clusters)]
setorder(selected_dt, dominant_cluster, -representative_mean_diff_vs_others, -representative_mean_auc)
fwrite(selected_dt, file.path(tables_dir, "02_regulon_heatmap_selected_regulons.csv"))

regulon_annotation_dt <- unique(activity_dt[, .(
  regulon,
  dominant_cluster,
  dominant_mean_auc,
  mean_range_across_clusters
)])
regulon_annotation_dt[, selected_for_heatmap := regulon %in% selected_dt$regulon]

heatmap_mean_mat <- cluster_mean_mat[selected_dt$regulon, present_clusters, drop = FALSE]
heatmap_z_mat <- zscore_rows(heatmap_mean_mat)

heatmap_dt <- as.data.table(heatmap_z_mat, keep.rownames = "regulon")
heatmap_dt <- melt(
  heatmap_dt,
  id.vars = "regulon",
  variable.name = "cluster",
  value.name = "z_score"
)

selected_order <- rev(selected_dt$regulon)
heatmap_dt[, regulon := factor(regulon, levels = selected_order)]
heatmap_dt[, cluster := factor(cluster, levels = present_clusters)]

heatmap_plot <- ggplot(heatmap_dt, aes(x = cluster, y = regulon, fill = z_score)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-2, 2),
    name = "Row z-score"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "C  Cluster-level regulon activity across progression-related epithelial states"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

heatmap_height <- max(6, min(14, 2 + 0.22 * nrow(selected_dt)))

ggsave(
  file.path(figures_dir, "Figure5_regulon_heatmap_main.pdf"),
  heatmap_plot,
  width = 7.2,
  height = heatmap_height,
  dpi = 300
)
ggsave(
  file.path(figures_dir, "Figure5_regulon_heatmap_main.png"),
  heatmap_plot,
  width = 7.2,
  height = heatmap_height,
  dpi = 300
)

cat("Computing cell-level regulon correlations with CNV burden and malignant program score...\n")
cell_meta[, malignant_program_score := as.numeric(malignant_program_score)]
cell_meta[, raw_cnv_burden := as.numeric(raw_cnv_burden)]

cnv_corr_dt <- compute_regulon_correlations(
  auc_mat = auc_mat,
  score_vec = cell_meta$raw_cnv_burden,
  metric_label = "raw_cnv_burden"
)
program_corr_dt <- compute_regulon_correlations(
  auc_mat = auc_mat,
  score_vec = cell_meta$malignant_program_score,
  metric_label = "malignant_program_score"
)

correlation_dt <- merge(cnv_corr_dt, program_corr_dt, by = "regulon", all = TRUE)
correlation_dt <- merge(correlation_dt, regulon_annotation_dt, by = "regulon", all.x = TRUE)
correlation_dt[, shared_positive_score :=
  ifelse(is.na(rho_raw_cnv_burden), 0, pmax(rho_raw_cnv_burden, 0)) +
  ifelse(is.na(rho_malignant_program_score), 0, pmax(rho_malignant_program_score, 0))
]
correlation_dt[, shared_negative_score :=
  ifelse(is.na(rho_raw_cnv_burden), 0, pmax(-rho_raw_cnv_burden, 0)) +
  ifelse(is.na(rho_malignant_program_score), 0, pmax(-rho_malignant_program_score, 0))
]
correlation_dt[, cnv_positive_rank := frank(-rho_raw_cnv_burden, ties.method = "first", na.last = "keep")]
correlation_dt[, program_positive_rank := frank(-rho_malignant_program_score, ties.method = "first", na.last = "keep")]
correlation_dt[, shared_positive_rank := cnv_positive_rank + program_positive_rank]
correlation_dt[, cnv_significant := !is.na(fdr_raw_cnv_burden) & fdr_raw_cnv_burden < 0.05]
correlation_dt[, program_significant := !is.na(fdr_malignant_program_score) & fdr_malignant_program_score < 0.05]
correlation_dt[, shared_significant_positive :=
  cnv_significant & program_significant &
  rho_raw_cnv_burden > 0 & rho_malignant_program_score > 0
]

setorder(
  correlation_dt,
  shared_positive_rank,
  -shared_positive_score,
  -mean_range_across_clusters,
  -dominant_mean_auc
)
fwrite(
  correlation_dt,
  file.path(tables_dir, "02_regulon_correlation_with_cnv_program.csv")
)

top_shared_dt <- copy(correlation_dt)[
  order(
    -shared_significant_positive,
    shared_positive_rank,
    -shared_positive_score,
    -mean_range_across_clusters,
    -dominant_mean_auc
  )
]
top_shared_dt <- top_shared_dt[, .(
  regulon,
  dominant_cluster,
  selected_for_heatmap,
  n_cells_raw_cnv_burden,
  rho_raw_cnv_burden,
  fdr_raw_cnv_burden,
  n_cells_malignant_program_score,
  rho_malignant_program_score,
  fdr_malignant_program_score,
  shared_positive_score,
  shared_positive_rank,
  shared_significant_positive,
  mean_range_across_clusters,
  dominant_mean_auc
)]
fwrite(
  head(top_shared_dt, 20L),
  file.path(tables_dir, "02_regulon_correlation_top_hits.csv")
)

label_regulons <- head(top_shared_dt$regulon, 12L)
correlation_plot_dt <- copy(correlation_dt)
correlation_plot_dt[, dominant_cluster := factor(dominant_cluster, levels = present_clusters)]
correlation_plot_dt[, label := ifelse(regulon %in% label_regulons, regulon, NA_character_)]

cluster_colors <- c(
  E10 = "#1b9e77",
  E11 = "#d95f02",
  E13 = "#7570b3",
  E17 = "#e7298a",
  E18 = "#66a61e",
  E19 = "#e6ab02"
)

correlation_plot <- ggplot(
  correlation_plot_dt,
  aes(
    x = rho_malignant_program_score,
    y = rho_raw_cnv_burden,
    color = dominant_cluster,
    shape = selected_for_heatmap,
    size = mean_range_across_clusters
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_point(alpha = 0.9) +
  geom_text(
    data = correlation_plot_dt[!is.na(label)][order(shared_positive_rank)],
    aes(label = label),
    size = 3,
    vjust = -0.8,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(values = cluster_colors, drop = FALSE) +
  scale_shape_manual(
    values = c(`FALSE` = 16, `TRUE` = 17),
    name = "Selected for\nheatmap"
  ) +
  scale_size_continuous(
    range = c(2.2, 5.0),
    name = "Mean range\nacross clusters"
  ) +
  labs(
    x = "Spearman rho with malignant_program_score",
    y = "Spearman rho with raw_cnv_burden",
    title = "D  Regulon activity tracks malignant program and CNV burden"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(figures_dir, "Figure5_regulon_correlation_scatter.pdf"),
  correlation_plot,
  width = 7.2,
  height = 6.2,
  dpi = 300
)
ggsave(
  file.path(figures_dir, "Figure5_regulon_correlation_scatter.png"),
  correlation_plot,
  width = 7.2,
  height = 6.2,
  dpi = 300
)

saveRDS(
  list(
    auc_matrix = auc_mat,
    cell_metadata = cell_meta,
    cluster_mean_auc = cluster_mean_mat,
    cluster_median_auc = cluster_median_mat,
    cluster_pct_active = cluster_pct_active_mat,
    cluster_cell_counts = cluster_cell_counts,
    activity_table = activity_dt,
    selected_regulons = selected_dt,
    regulon_annotations = regulon_annotation_dt,
    regulon_correlation_table = correlation_dt,
    regulon_correlation_top_hits = top_shared_dt
  ),
  file.path(rds_dir, "02_regulon_auc_matrix.rds")
)

cat("Regulon cluster summary completed.\n")
cat("Cells used:", ncol(auc_mat), "\n")
cat("Regulons summarized:", nrow(auc_mat), "\n")
cat("Selected regulons for heatmap:", nrow(selected_dt), "\n")
cat("Main heatmap:", file.path(figures_dir, "Figure5_regulon_heatmap_main.pdf"), "\n")
cat("Correlation scatter:", file.path(figures_dir, "Figure5_regulon_correlation_scatter.pdf"), "\n")

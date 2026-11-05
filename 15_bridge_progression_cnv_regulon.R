args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  Sys.getenv("HP_GC_PROJECT_ROOT", unset = "/Users/huahaijiao/Desktop/02_study/0206_artical/Hp_gastric_cancer")
}

result_root <- file.path(
  project_root,
  "gastric_cancer", "results", "Figure5_regulon_pathway_bulk"
)
figure2_root <- file.path(
  project_root,
  "gastric_cancer", "results", "figure2_epithelial_progression_pseudotime"
)
tables_dir <- file.path(result_root, "tables")
figures_dir <- file.path(result_root, "figures")
manuscript_dir <- file.path(result_root, "manuscript")
rds_dir <- file.path(result_root, "rds")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
})

cluster_order <- c("E10", "E11", "E13", "E17", "E18", "E19")
requested_key_regulons <- c("MAX(+)", "NR1H4(+)", "HNF4G(+)")
key_regulons <- requested_key_regulons
cluster_colors <- c(
  E10 = "#1b9e77",
  E11 = "#d95f02",
  E13 = "#7570b3",
  E17 = "#e7298a",
  E18 = "#66a61e",
  E19 = "#e6ab02"
)

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 10) {
    return(data.table(n = length(x), rho = NA_real_, pvalue = NA_real_))
  }
  test <- tryCatch(
    suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(test)) {
    return(data.table(n = length(x), rho = NA_real_, pvalue = NA_real_))
  }
  data.table(
    n = length(x),
    rho = unname(test$estimate),
    pvalue = test$p.value
  )
}

base_theme <- theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

required_files <- c(
  file.path(rds_dir, "02_regulon_auc_matrix.rds"),
  file.path(tables_dir, "01_pyscenic_cell_metadata_cluster_aware.csv"),
  file.path(figure2_root, "tables", "gse134520_progression_cells_pseudotime_metadata.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

auc_obj <- readRDS(file.path(rds_dir, "02_regulon_auc_matrix.rds"))
cell_meta_dt <- fread(file.path(tables_dir, "01_pyscenic_cell_metadata_cluster_aware.csv"))
pseudotime_dt <- fread(file.path(
  figure2_root, "tables", "gse134520_progression_cells_pseudotime_metadata.csv"
))

setnames(pseudotime_dt, "cell_barcode", "cell_id")

merge_cols <- c(
  "cell_id", "pseudotime_raw", "pseudotime",
  "gain_union_score", "loss_union_score", "identity_loss_score",
  "UMAP1", "UMAP2"
)
pseudotime_dt <- pseudotime_dt[, ..merge_cols]

bridge_meta <- merge(
  cell_meta_dt,
  pseudotime_dt,
  by = "cell_id",
  all.x = TRUE
)
bridge_meta[, epi_cluster_id := factor(epi_cluster_id, levels = cluster_order)]

auc_mat <- auc_obj$auc_matrix
key_regulons <- requested_key_regulons[requested_key_regulons %in% rownames(auc_mat)]
missing_key_regulons <- setdiff(requested_key_regulons, key_regulons)
if (length(missing_key_regulons) > 0) {
  warning("Skipping key regulons absent from the cluster-aware pySCENIC result: ",
          paste(missing_key_regulons, collapse = ", "))
}
if (length(key_regulons) == 0) {
  stop("None of the requested key regulons are present in the regulon AUC matrix.")
}
if (!identical(colnames(auc_mat), bridge_meta$cell_id)) {
  bridge_meta <- bridge_meta[match(colnames(auc_mat), cell_id)]
  if (anyNA(bridge_meta$cell_id)) {
    stop("Failed to align bridge metadata to regulon AUC matrix")
  }
}

feature_list <- c(
  "malignant_program_score",
  "raw_cnv_burden",
  "pseudotime",
  "gain_union_score",
  "identity_loss_score"
)

cell_level_cor_dt <- rbindlist(lapply(rownames(auc_mat), function(reg) {
  auc_vec <- as.numeric(auc_mat[reg, ])
  rbindlist(lapply(feature_list, function(feat) {
    stat <- safe_spearman(auc_vec, bridge_meta[[feat]])
    data.table(
      regulon = reg,
      feature = feat,
      n = stat$n,
      rho = stat$rho,
      pvalue = stat$pvalue
    )
  }))
}), fill = TRUE)
cell_level_cor_dt[, fdr := p.adjust(pvalue, method = "BH"), by = feature]

key_cell_cor_dt <- cell_level_cor_dt[regulon %in% key_regulons]
fwrite(
  cell_level_cor_dt,
  file.path(tables_dir, "06_bridge_regulon_vs_progression_features_cell_level.csv")
)
fwrite(
  key_cell_cor_dt,
  file.path(tables_dir, "06_bridge_key_regulon_vs_progression_features_cell_level.csv")
)

regulon_long_dt <- rbindlist(lapply(key_regulons, function(reg) {
  data.table(
    cell_id = colnames(auc_mat),
    regulon = reg,
    regulon_auc = as.numeric(auc_mat[reg, ])
  )
}))

bridge_long_dt <- merge(regulon_long_dt, bridge_meta, by = "cell_id", all.x = TRUE)

cluster_level_dt <- bridge_long_dt[, .(
  n_cells = .N,
  mean_regulon_auc = mean(regulon_auc, na.rm = TRUE),
  mean_malignant_program_score = mean(malignant_program_score, na.rm = TRUE),
  mean_raw_cnv_burden = mean(raw_cnv_burden, na.rm = TRUE),
  mean_pseudotime = mean(pseudotime, na.rm = TRUE),
  mean_gain_union_score = mean(gain_union_score, na.rm = TRUE),
  mean_identity_loss_score = mean(identity_loss_score, na.rm = TRUE)
), by = .(regulon, epi_cluster_id)]

fwrite(
  cluster_level_dt,
  file.path(tables_dir, "06_bridge_key_regulon_cluster_level_summary.csv")
)

cluster_feature_cor_dt <- rbindlist(lapply(key_regulons, function(reg) {
  sub <- cluster_level_dt[regulon == reg]
  rbindlist(lapply(c(
    "mean_malignant_program_score",
    "mean_raw_cnv_burden",
    "mean_pseudotime",
    "mean_gain_union_score",
    "mean_identity_loss_score"
  ), function(feat) {
    stat <- safe_spearman(sub$mean_regulon_auc, sub[[feat]])
    data.table(
      regulon = reg,
      feature = feat,
      n = stat$n,
      rho = stat$rho,
      pvalue = stat$pvalue
    )
  }))
}))
cluster_feature_cor_dt[, fdr := p.adjust(pvalue, method = "BH")]

fwrite(
  cluster_feature_cor_dt,
  file.path(tables_dir, "06_bridge_key_regulon_vs_progression_features_cluster_level.csv")
)

heatmap_feature_labels <- c(
  malignant_program_score = "Malignant\nprogram",
  raw_cnv_burden = "Raw CNV\nburden",
  pseudotime = "Pseudotime",
  gain_union_score = "Gain-union\nscore",
  identity_loss_score = "Identity-loss\nscore"
)

heatmap_dt <- copy(key_cell_cor_dt)
heatmap_dt[, feature_label := heatmap_feature_labels[feature]]
heatmap_dt[, regulon := factor(regulon, levels = rev(key_regulons))]
heatmap_dt[, feature_label := factor(
  feature_label,
  levels = heatmap_feature_labels[feature_list]
)]
heatmap_dt[, rho_label := sprintf("%.2f", rho)]

p_heatmap <- ggplot(heatmap_dt, aes(x = feature_label, y = regulon, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(aes(label = rho_label), size = 3) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-0.5, 0.5),
    name = "Spearman rho"
  ) +
  labs(
    title = "A  Cell-level bridge between progression, CNV, and regulons",
    x = NULL,
    y = NULL
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

scatter_dt <- bridge_meta[, .(
  cell_id,
  epi_cluster_id,
  malignant_program_score,
  raw_cnv_burden,
  pseudotime,
  gain_union_score,
  identity_loss_score
)]
for (reg in key_regulons) {
  scatter_dt[[reg]] <- as.numeric(auc_mat[reg, scatter_dt$cell_id])
}

make_scatter <- function(dt, xcol, ycol, title_txt, xlab, ylab) {
  ggplot(dt[is.finite(get(xcol)) & is.finite(get(ycol))], aes(
    x = get(xcol),
    y = get(ycol),
    color = epi_cluster_id
  )) +
    geom_point(size = 0.5, alpha = 0.55) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
    scale_color_manual(values = cluster_colors, drop = FALSE) +
    labs(title = title_txt, x = xlab, y = ylab, color = "Cluster") +
    base_theme +
    theme(legend.position = "right")
}

scatter_specs <- data.table(
  regulon = c("NR1H4(+)", "HNF4G(+)", "MAX(+)"),
  xcol = c("malignant_program_score", "gain_union_score", "malignant_program_score"),
  title_txt = c(
    "B  NR1H4(+) rises with malignant program",
    "C  HNF4G(+) tracks acquired progression program",
    "D  MAX(+) links malignant program activity"
  ),
  xlab = c(
    "Malignant program score",
    "Gain-union score",
    "Malignant program score"
  )
)
scatter_specs <- scatter_specs[regulon %in% key_regulons]
scatter_specs <- head(scatter_specs, 3L)
scatter_plots <- lapply(seq_len(nrow(scatter_specs)), function(i) {
  make_scatter(
    scatter_dt,
    xcol = scatter_specs$xcol[i],
    ycol = scatter_specs$regulon[i],
    title_txt = scatter_specs$title_txt[i],
    xlab = scatter_specs$xlab[i],
    ylab = paste0(scatter_specs$regulon[i], " regulon AUC")
  )
})
while (length(scatter_plots) < 3L) {
  scatter_plots[[length(scatter_plots) + 1L]] <- ggplot() +
    theme_void(base_size = 10) +
    labs(title = "Key regulon unavailable")
}

cluster_plot_dt <- melt(
  cluster_level_dt,
  id.vars = c("regulon", "epi_cluster_id", "n_cells", "mean_regulon_auc"),
  measure.vars = c(
    "mean_malignant_program_score",
    "mean_raw_cnv_burden",
    "mean_pseudotime",
    "mean_gain_union_score",
    "mean_identity_loss_score"
  ),
  variable.name = "feature",
  value.name = "cluster_feature_mean"
)
cluster_plot_dt[, feature := factor(feature, levels = c(
  "mean_malignant_program_score",
  "mean_raw_cnv_burden",
  "mean_pseudotime",
  "mean_gain_union_score",
  "mean_identity_loss_score"
), labels = c(
  "Malignant program",
  "Raw CNV burden",
  "Pseudotime",
  "Gain-union score",
  "Identity-loss score"
))]

p_cluster <- ggplot(
  cluster_plot_dt,
  aes(x = cluster_feature_mean, y = mean_regulon_auc, color = epi_cluster_id)
) +
  geom_point(size = 2.2, alpha = 0.95) +
  geom_smooth(method = "lm", se = FALSE, color = "grey35", linewidth = 0.45) +
  facet_grid(regulon ~ feature, scales = "free") +
  scale_color_manual(values = cluster_colors, drop = FALSE) +
  labs(
    title = "E  Cluster-level bridge summary",
    x = "Cluster-level progression / CNV feature mean",
    y = "Cluster-level mean regulon AUC",
    color = "Cluster"
  ) +
  base_theme +
  theme(
    legend.position = "right",
    strip.text = element_text(size = 7)
  )

bridge_figure <- (p_heatmap | scatter_plots[[1]]) /
  (scatter_plots[[2]] | scatter_plots[[3]]) /
  p_cluster +
  plot_layout(heights = c(0.95, 0.95, 1.3))

ggsave(
  file.path(figures_dir, "Figure5_bridge_progression_cnv_regulon_main.pdf"),
  bridge_figure,
  width = 13.5,
  height = 15.5,
  dpi = 320
)
ggsave(
  file.path(figures_dir, "Figure5_bridge_progression_cnv_regulon_main.png"),
  bridge_figure,
  width = 13.5,
  height = 15.5,
  dpi = 320
)

writeLines(
  c(
    "# Figure 5 bridge analysis note",
    "",
    "This bridge analysis links the Figure 1-3 progression axis to the Figure 4-5 CNV/regulon axis at the single-cell level.",
    "",
    "Main inputs:",
    "- Figure 2 pseudotime metadata",
    "- Figure 5 pySCENIC regulon AUC matrix",
    "- Figure 4/5 raw CNV burden metadata",
    "",
    "Main outputs:",
    "- Cell-level regulon vs progression/CNV correlation table",
    "- Cluster-level key regulon bridge summary",
    "- Main bridge figure"
  ),
  file.path(manuscript_dir, "Figure5_bridge_progression_cnv_regulon_note.md")
)

cat("Bridge analysis completed.\n")
cat("Cell-level table:", file.path(tables_dir, "06_bridge_regulon_vs_progression_features_cell_level.csv"), "\n")
cat("Cluster-level table:", file.path(tables_dir, "06_bridge_key_regulon_cluster_level_summary.csv"), "\n")
cat("Main figure:", file.path(figures_dir, "Figure5_bridge_progression_cnv_regulon_main.pdf"), "\n")

project_root <- "/Users/huahaijiao/Desktop/02_study/0206_artical/Hp_gastric_cancer"
result_root <- file.path(
  project_root,
  "gastric_cancer", "results", "Figure4_malignant_program_mechanism"
)
tables_dir <- file.path(result_root, "tables")
figures_dir <- file.path(result_root, "figures")
logs_dir <- file.path(result_root, "logs")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

raw_cna_path <- "/Users/huahaijiao/Documents/New project 3/GSE134520_epithelial_copykat_CNA_raw_results_gene_by_cell.txt"
meta_path <- file.path(tables_dir, "01_epithelial_cells_for_cnv_metadata.csv")
input_summary_path <- file.path(tables_dir, "02_copykat_input_summary.csv")

required_files <- c(raw_cna_path, meta_path, input_summary_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

meta <- fread(meta_path)
meta[, cell_id := as.character(cell_id)]
input_summary <- fread(input_summary_path)

raw_cna <- fread(raw_cna_path, showProgress = TRUE)
meta_cols <- c(
  "abspos", "chromosome_name", "start_position", "end_position",
  "ensembl_gene_id", "hgnc_symbol", "band"
)
cell_cols <- setdiff(colnames(raw_cna), meta_cols)

if (length(cell_cols) == 0) {
  stop("No cell columns found in raw copykat CNA matrix: ", raw_cna_path)
}

cnv_mat <- as.matrix(raw_cna[, ..cell_cols])
storage.mode(cnv_mat) <- "double"

cnv_burden <- colMeans(abs(cnv_mat), na.rm = TRUE)
cnv_gain_fraction_005 <- colMeans(cnv_mat > 0.05, na.rm = TRUE)
cnv_loss_fraction_005 <- colMeans(cnv_mat < -0.05, na.rm = TRUE)
cnv_shift_mean <- colMeans(cnv_mat, na.rm = TRUE)

cell_dt <- data.table(
  cell_id = cell_cols,
  raw_cnv_burden = as.numeric(cnv_burden),
  raw_cnv_gain_fraction_005 = as.numeric(cnv_gain_fraction_005),
  raw_cnv_loss_fraction_005 = as.numeric(cnv_loss_fraction_005),
  raw_cnv_shift_mean = as.numeric(cnv_shift_mean)
)

cell_dt <- merge(meta, cell_dt, by = "cell_id", all.y = TRUE)

ref_burden <- cell_dt[cnv_reference_flag == TRUE, raw_cnv_burden]
if (length(ref_burden) < 20) {
  stop("Too few reference cells retained in raw copykat output")
}

reference_threshold_q95 <- as.numeric(quantile(ref_burden, probs = 0.95, na.rm = TRUE))
reference_threshold_median_mad3 <- median(ref_burden, na.rm = TRUE) +
  3 * mad(ref_burden, constant = 1, na.rm = TRUE)

cell_dt[, raw_cnv_state_ref95 := fifelse(
  raw_cnv_burden > reference_threshold_q95,
  "high_cnv_like",
  "reference_like"
)]
cell_dt[, raw_cnv_state_mad3 := fifelse(
  raw_cnv_burden > reference_threshold_median_mad3,
  "high_cnv_like",
  "reference_like"
)]

cell_dt[, cnv_result_source := "copykat_raw_matrix_salvage"]

cluster_coverage <- meta[, .(input_n = .N), by = .(epi_cluster_id)]
cluster_output <- cell_dt[, .(output_n = .N), by = .(epi_cluster_id)]
cluster_summary <- merge(cluster_coverage, cluster_output, by = "epi_cluster_id", all = TRUE)
cluster_summary[is.na(input_n), input_n := 0L]
cluster_summary[is.na(output_n), output_n := 0L]

cluster_stats <- cell_dt[, .(
  n_cells = .N,
  lesion_group_major = names(sort(table(lesion_group), decreasing = TRUE))[1],
  raw_cnv_burden_median = median(raw_cnv_burden, na.rm = TRUE),
  raw_cnv_burden_mean = mean(raw_cnv_burden, na.rm = TRUE),
  raw_cnv_gain_fraction_005_median = median(raw_cnv_gain_fraction_005, na.rm = TRUE),
  raw_cnv_loss_fraction_005_median = median(raw_cnv_loss_fraction_005, na.rm = TRUE),
  raw_cnv_shift_mean_median = median(raw_cnv_shift_mean, na.rm = TRUE),
  high_cnv_like_fraction_ref95 = mean(raw_cnv_state_ref95 == "high_cnv_like", na.rm = TRUE),
  high_cnv_like_fraction_mad3 = mean(raw_cnv_state_mad3 == "high_cnv_like", na.rm = TRUE),
  progression_cluster_flag = unique(progression_cluster_flag)[1],
  cnv_reference_flag = unique(cnv_reference_flag)[1]
), by = .(epi_cluster_id)]

cluster_summary <- merge(cluster_summary, cluster_stats, by = "epi_cluster_id", all = TRUE)
cluster_summary[, output_rate_vs_input := fifelse(input_n > 0, output_n / input_n, NA_real_)]

group_summary <- cell_dt[, .(
  n_cells = .N,
  raw_cnv_burden_median = median(raw_cnv_burden, na.rm = TRUE),
  raw_cnv_burden_mean = mean(raw_cnv_burden, na.rm = TRUE),
  high_cnv_like_fraction_ref95 = mean(raw_cnv_state_ref95 == "high_cnv_like", na.rm = TRUE),
  high_cnv_like_fraction_mad3 = mean(raw_cnv_state_mad3 == "high_cnv_like", na.rm = TRUE)
), by = .(lesion_group)][order(factor(lesion_group, levels = c("NAG", "CAG", "IM", "EGC")))]

progression_vs_reference_summary <- cell_dt[, .(
  n_cells = .N,
  raw_cnv_burden_median = median(raw_cnv_burden, na.rm = TRUE),
  raw_cnv_burden_mean = mean(raw_cnv_burden, na.rm = TRUE),
  high_cnv_like_fraction_ref95 = mean(raw_cnv_state_ref95 == "high_cnv_like", na.rm = TRUE),
  high_cnv_like_fraction_mad3 = mean(raw_cnv_state_mad3 == "high_cnv_like", na.rm = TRUE)
), by = .(cnv_group = fifelse(
  cnv_reference_flag == TRUE, "reference",
  fifelse(progression_cluster_flag == TRUE, "progression", "other")
))][order(match(cnv_group, c("reference", "progression", "other")))]

threshold_summary <- data.table(
  metric = c("reference_burden_q95", "reference_burden_median_plus_3mad"),
  value = c(reference_threshold_q95, reference_threshold_median_mad3),
  reference_cells_retained = length(ref_burden),
  total_cells_retained = nrow(cell_dt),
  input_cells_requested = input_summary$total_cells[1]
)

cluster_levels <- unique(meta$epi_cluster_id)
cluster_levels <- cluster_levels[!is.na(cluster_levels)]
cell_dt[, epi_cluster_id := factor(epi_cluster_id, levels = cluster_levels)]
cluster_summary[, epi_cluster_id := factor(epi_cluster_id, levels = cluster_levels)]

p1 <- ggplot(
  cell_dt[epi_cluster_id %in% unique(cluster_summary$epi_cluster_id)],
  aes(x = epi_cluster_id, y = raw_cnv_burden, fill = progression_cluster_flag)
) +
  geom_boxplot(outlier.size = 0.15, width = 0.65) +
  scale_fill_manual(values = c("FALSE" = "#B8C4D6", "TRUE" = "#C95F3C")) +
  labs(
    x = NULL,
    y = "Raw CNV burden",
    fill = "Progression cluster"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

p2 <- ggplot(
  cluster_summary[!is.na(high_cnv_like_fraction_ref95)],
  aes(x = epi_cluster_id, y = high_cnv_like_fraction_ref95, fill = progression_cluster_flag)
) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("FALSE" = "#8FA8C9", "TRUE" = "#D96B45")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = NULL,
    y = "High-CNV-like fraction",
    fill = "Progression cluster"
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

ggsave(
  filename = file.path(figures_dir, "03_copykat_raw_matrix_cnv_burden_by_cluster.pdf"),
  plot = p1, width = 8.8, height = 4.8, dpi = 300
)
ggsave(
  filename = file.path(figures_dir, "03_copykat_raw_matrix_cnv_burden_by_cluster.png"),
  plot = p1, width = 8.8, height = 4.8, dpi = 300
)
ggsave(
  filename = file.path(figures_dir, "03_copykat_raw_matrix_high_cnv_fraction_by_cluster.pdf"),
  plot = p2, width = 8.8, height = 4.8, dpi = 300
)
ggsave(
  filename = file.path(figures_dir, "03_copykat_raw_matrix_high_cnv_fraction_by_cluster.png"),
  plot = p2, width = 8.8, height = 4.8, dpi = 300
)

fwrite(cell_dt, file.path(tables_dir, "03_copykat_raw_matrix_cell_level_results.csv"))
fwrite(cluster_summary, file.path(tables_dir, "03_copykat_raw_matrix_cluster_summary.csv"))
fwrite(group_summary, file.path(tables_dir, "03_copykat_raw_matrix_group_summary.csv"))
fwrite(
  progression_vs_reference_summary,
  file.path(tables_dir, "03_copykat_raw_matrix_progression_vs_reference_summary.csv")
)
fwrite(threshold_summary, file.path(tables_dir, "03_copykat_raw_matrix_threshold_summary.csv"))

cat("Raw copykat matrix salvage summary completed.\n")
cat("Raw CNA matrix:", raw_cna_path, "\n")
cat("Retained cells:", nrow(cell_dt), "\n")
cat("Reference cells retained:", length(ref_burden), "\n")
cat("Reference 95th percentile threshold:", reference_threshold_q95, "\n")
cat("Saved cell-level results:", file.path(tables_dir, "03_copykat_raw_matrix_cell_level_results.csv"), "\n")

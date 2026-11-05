project_root <- "/Users/huahaijiao/Desktop/02_study/0206_artical/Hp_gastric_cancer"
result_root <- file.path(
  project_root,
  "gastric_cancer", "results", "Figure4_malignant_program_mechanism"
)
tables_dir <- file.path(result_root, "tables")
figures_dir <- file.path(result_root, "figures")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

cell_result_path <- file.path(tables_dir, "03_copykat_raw_matrix_cell_level_results.csv")
required_files <- c(cell_result_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

dt <- fread(cell_result_path)
dt[, lesion_group := factor(as.character(lesion_group), levels = c("NAG", "CAG", "IM", "EGC"))]
dt[, epi_cluster_id := factor(as.character(epi_cluster_id), levels = sort(unique(as.character(epi_cluster_id))))]

progression_dt <- dt[progression_cluster_flag == TRUE]

cor_dt <- data.table(
  metric_x = c(
    "pseudotime",
    "progression_score_47g",
    "core_progression_score_12g_auxiliary",
    "malignant_program_score",
    "gastric_homeostasis_score"
  ),
  rho = c(
    suppressWarnings(cor(progression_dt$raw_cnv_burden, progression_dt$pseudotime, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(progression_dt$raw_cnv_burden, progression_dt$progression_score_47g, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(progression_dt$raw_cnv_burden, progression_dt$core_progression_score_12g, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(progression_dt$raw_cnv_burden, progression_dt$malignant_program_score, method = "spearman", use = "complete.obs")),
    suppressWarnings(cor(progression_dt$raw_cnv_burden, progression_dt$gastric_homeostasis_score, method = "spearman", use = "complete.obs"))
  )
)

fwrite(cor_dt, file.path(tables_dir, "04_cnv_vs_progression_correlations.csv"))

p1 <- ggplot(
  progression_dt[!is.na(pseudotime)],
  aes(x = pseudotime, y = raw_cnv_burden, color = epi_cluster_id)
) +
  geom_point(size = 0.35, alpha = 0.35) +
  geom_smooth(method = "loess", se = FALSE, color = "black", linewidth = 0.9) +
  labs(x = "Pseudotime", y = "Raw CNV burden", color = "Cluster") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

p2 <- ggplot(
  progression_dt[!is.na(progression_score_47g)],
  aes(x = progression_score_47g, y = raw_cnv_burden, color = epi_cluster_id)
) +
  geom_point(size = 0.35, alpha = 0.35) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.9) +
  labs(x = "47-gene progression score", y = "Raw CNV burden", color = "Cluster") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

p3 <- ggplot(
  dt[!is.na(lesion_group)],
  aes(x = lesion_group, y = raw_cnv_burden, fill = lesion_group)
) +
  geom_boxplot(outlier.size = 0.15, width = 0.65) +
  scale_fill_manual(values = c("NAG" = "#8FA8C9", "CAG" = "#B8C4D6", "IM" = "#D59B5B", "EGC" = "#C95F3C")) +
  labs(x = NULL, y = "Raw CNV burden") +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")

ggsave(
  file.path(figures_dir, "04_cnv_burden_vs_pseudotime.pdf"),
  p1, width = 6.6, height = 5.0, dpi = 300
)
ggsave(
  file.path(figures_dir, "04_cnv_burden_vs_pseudotime.png"),
  p1, width = 6.6, height = 5.0, dpi = 300
)
ggsave(
  file.path(figures_dir, "04_cnv_burden_vs_progression_score_47g.pdf"),
  p2, width = 6.6, height = 5.0, dpi = 300
)
ggsave(
  file.path(figures_dir, "04_cnv_burden_vs_progression_score_47g.png"),
  p2, width = 6.6, height = 5.0, dpi = 300
)
ggsave(
  file.path(figures_dir, "04_cnv_burden_by_lesion_group.pdf"),
  p3, width = 5.8, height = 4.6, dpi = 300
)
ggsave(
  file.path(figures_dir, "04_cnv_burden_by_lesion_group.png"),
  p3, width = 5.8, height = 4.6, dpi = 300
)

cat("CNV vs progression analysis completed.\n")
print(cor_dt)

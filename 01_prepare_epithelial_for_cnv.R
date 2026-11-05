project_root <- "/Users/huahaijiao/Desktop/02_study/0206_artical/Hp_gastric_cancer"
result_root <- file.path(
  project_root,
  "gastric_cancer", "results", "Figure4_malignant_program_mechanism"
)
scripts_dir <- file.path(result_root, "scripts")
tables_dir <- file.path(result_root, "tables")
rds_dir <- file.path(result_root, "rds")
logs_dir <- file.path(result_root, "logs")

dir.create(scripts_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(data.table)
  library(Matrix)
})

normalize_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

mean_score <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(mat)))
  }
  sub_mat <- as.matrix(mat[genes, , drop = FALSE])
  row_means <- rowMeans(sub_mat, na.rm = TRUE)
  row_sds <- apply(sub_mat, 1, sd, na.rm = TRUE)
  row_sds[!is.finite(row_sds) | row_sds == 0] <- 1
  zmat <- sweep(sub_mat, 1, row_means, "-")
  zmat <- sweep(zmat, 1, row_sds, "/")
  colMeans(zmat, na.rm = TRUE)
}

annotated_sce_path <- file.path(
  project_root,
  "gastric_cancer", "results", "GSE134520", "epithelial", "rds",
  "gse134520_epithelial_sce_reclustered_annotated.rds"
)
pseudotime_meta_path <- file.path(
  project_root,
  "gastric_cancer", "results", "GSE134520", "epithelial", "pseudotime_main", "tables",
  "gse134520_progression_cells_pseudotime_metadata.csv"
)
signature_path <- file.path(
  project_root,
  "gastric_cancer", "results", "core_progression_score_dfs_priority",
  "biology_constrained_round2", "tables", "02_biology_constrained_core_signature.csv"
)
main47_signature_path <- file.path(
  project_root,
  "gastric_cancer", "results", "progression_score_main47",
  "tables", "progression_score_47g_definition.csv"
)

required_files <- c(
  annotated_sce_path,
  pseudotime_meta_path,
  signature_path,
  main47_signature_path
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

sce <- readRDS(annotated_sce_path)
meta <- as.data.table(as.data.frame(SummarizedExperiment::colData(sce)))
meta[, cell_barcode := as.character(cell_barcode)]
meta <- meta[major_lineage == "Epithelial"]

pseudotime_meta <- fread(pseudotime_meta_path)
if (!"cell_barcode" %in% colnames(pseudotime_meta)) {
  stop("pseudotime metadata must contain a cell_barcode column: ", pseudotime_meta_path)
}
if (anyDuplicated(pseudotime_meta$cell_barcode)) {
  pseudotime_meta <- pseudotime_meta[!duplicated(cell_barcode)]
}

signature_dt <- fread(signature_path)
core_up <- normalize_gene(signature_dt[direction == "up_gain", gene])
core_down <- normalize_gene(signature_dt[direction == "down_loss", gene])
main47_dt <- fread(main47_signature_path)
main47_up <- normalize_gene(main47_dt[direction == "up_gain", gene])
main47_down <- normalize_gene(main47_dt[direction == "down_loss", gene])

rownames(sce) <- normalize_gene(rownames(sce))
sce <- sce[, match(meta$cell_barcode, colnames(sce))]

log_mat <- logcounts(sce)
rownames(log_mat) <- normalize_gene(rownames(log_mat))
log_mat <- log_mat[!duplicated(rownames(log_mat)), , drop = FALSE]
core_progression_score_12g <- mean_score(log_mat, core_up) - mean_score(log_mat, core_down)
progression_score_47g <- mean_score(log_mat, main47_up) - mean_score(log_mat, main47_down)

meta[, progression_score_47g := progression_score_47g[match(cell_barcode, colnames(sce))]]
meta[, core_progression_score_12g := core_progression_score_12g[match(cell_barcode, colnames(sce))]]
meta[, cell_name := cell_barcode]

progression_clusters <- c("E10", "E11", "E13", "E17", "E18", "E19")
reference_clusters <- c("E01", "E02", "E03", "E04", "E05", "E06")

meta[, progression_cluster_flag := epi_cluster_id %in% progression_clusters]
meta[, cnv_reference_flag := lesion_group %in% c("NAG", "CAG") & epi_cluster_id %in% reference_clusters]
meta[, cnv_reference_group := fifelse(cnv_reference_flag, "reference_NAG_CAG_gastric", "non_reference")]

keep_cols <- intersect(
  c(
    "cell_barcode", "cell_name", "sample_id", "stage_code", "lesion_group",
    "cluster_id", "epi_cluster_id", "major_lineage", "detailed_celltype",
    "annotation_confidence", "malignant_program_score", "gastric_homeostasis_score",
    "epithelial_major_program", "epithelial_detailed_celltype", "progression_role",
    "progression_cluster_flag", "cnv_reference_flag", "cnv_reference_group",
    "progression_score_47g",
    "core_progression_score_12g"
  ),
  colnames(meta)
)

epi_meta <- copy(meta[, ..keep_cols])
setnames(epi_meta, old = "cell_barcode", new = "cell_id")

if ("cell_barcode" %in% colnames(pseudotime_meta)) {
  pt_keep <- intersect(
    c(
      "cell_barcode", "UMAP1", "UMAP2", "pseudotime_raw", "pseudotime",
      "gain_union_score", "loss_union_score", "identity_loss_score"
    ),
    colnames(pseudotime_meta)
  )
  pt_dt <- copy(pseudotime_meta[, ..pt_keep])
  setnames(pt_dt, "cell_barcode", "cell_id")
  epi_meta <- merge(epi_meta, pt_dt, by = "cell_id", all.x = TRUE)
}

epi_meta[, lesion_group := factor(as.character(lesion_group), levels = c("NAG", "CAG", "IM", "EGC"))]
epi_meta[, stage_code := factor(as.character(stage_code), levels = c("NAG", "CAG", "IMW", "IMS", "EGC"))]
epi_meta[, epi_cluster_id := factor(as.character(epi_cluster_id), levels = sprintf("E%02d", 1:19))]

epi_meta[, pseudotime_available := !is.na(pseudotime)]

epithelial_cells <- epi_meta$cell_id
sce_epi <- sce[, match(epithelial_cells, colnames(sce))]
SummarizedExperiment::colData(sce_epi)$progression_score_47g <- epi_meta$progression_score_47g
SummarizedExperiment::colData(sce_epi)$core_progression_score_12g <- epi_meta$core_progression_score_12g
SummarizedExperiment::colData(sce_epi)$progression_cluster_flag <- epi_meta$progression_cluster_flag
SummarizedExperiment::colData(sce_epi)$cnv_reference_flag <- epi_meta$cnv_reference_flag
SummarizedExperiment::colData(sce_epi)$cnv_reference_group <- epi_meta$cnv_reference_group
if ("pseudotime" %in% colnames(epi_meta)) {
  SummarizedExperiment::colData(sce_epi)$pseudotime <- epi_meta$pseudotime
}
if ("pseudotime_raw" %in% colnames(epi_meta)) {
  SummarizedExperiment::colData(sce_epi)$pseudotime_raw <- epi_meta$pseudotime_raw
}

summary_by_cluster <- epi_meta[, .(
  n_cells = .N,
  n_reference = sum(cnv_reference_flag, na.rm = TRUE),
  progression_cluster = unique(progression_cluster_flag)[1],
  pseudotime_available_n = sum(pseudotime_available, na.rm = TRUE),
  progression_score_47g_median = median(progression_score_47g, na.rm = TRUE),
  core_progression_score_median = median(core_progression_score_12g, na.rm = TRUE),
  malignant_program_score_median = median(malignant_program_score, na.rm = TRUE),
  gastric_homeostasis_score_median = median(gastric_homeostasis_score, na.rm = TRUE)
), by = .(epi_cluster_id)][order(epi_cluster_id)]

summary_by_group <- epi_meta[, .(
  n_cells = .N,
  n_reference = sum(cnv_reference_flag, na.rm = TRUE),
  n_progression_cluster = sum(progression_cluster_flag, na.rm = TRUE),
  progression_score_47g_median = median(progression_score_47g, na.rm = TRUE),
  core_progression_score_median = median(core_progression_score_12g, na.rm = TRUE)
), by = .(lesion_group)][order(lesion_group)]

reference_samples <- epi_meta[cnv_reference_flag == TRUE, .N, by = .(sample_id, lesion_group, epi_cluster_id)][order(sample_id, epi_cluster_id)]

fwrite(epi_meta, file.path(tables_dir, "01_epithelial_cells_for_cnv_metadata.csv"))
fwrite(summary_by_cluster, file.path(tables_dir, "01_epithelial_cluster_cnv_prep_summary.csv"))
fwrite(summary_by_group, file.path(tables_dir, "01_epithelial_group_cnv_prep_summary.csv"))
fwrite(reference_samples, file.path(tables_dir, "01_cnv_reference_cells_by_sample_cluster.csv"))

saveRDS(sce_epi, file.path(rds_dir, "01_epithelial_for_cnv.rds"))

cat("Prepared epithelial object for CNV analysis.\n")
cat("Cells:", ncol(sce_epi), "\n")
cat("Genes:", nrow(sce_epi), "\n")
cat("Reference cells:", sum(epi_meta$cnv_reference_flag, na.rm = TRUE), "\n")
cat("Progression-cluster cells:", sum(epi_meta$progression_cluster_flag, na.rm = TRUE), "\n")
cat("Output RDS:", file.path(rds_dir, "01_epithelial_for_cnv.rds"), "\n")
cat("Output metadata:", file.path(tables_dir, "01_epithelial_cells_for_cnv_metadata.csv"), "\n")

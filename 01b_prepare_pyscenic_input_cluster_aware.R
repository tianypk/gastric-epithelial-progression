project_root <- "/Users/huahaijiao/Desktop/02_study/0206_artical/Hp_gastric_cancer"
result_root <- file.path(
  project_root,
  "gastric_cancer", "results", "Figure5_regulon_pathway_bulk"
)
tables_dir <- file.path(result_root, "tables")
rds_dir <- file.path(result_root, "rds")
reference_dir <- file.path(result_root, "reference")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

normalize_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

input_rds_path <- file.path(rds_dir, "01_progression_epithelial_for_regulon.rds")
tf_path <- file.path(reference_dir, "allTFs_hg38.txt")

required_files <- c(input_rds_path, tf_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

obj <- readRDS(input_rds_path)
counts_mat <- obj$counts
meta_dt <- as.data.table(obj$metadata)
gene_meta_dt <- as.data.table(obj$gene_metadata)

if (!"cell_id" %in% colnames(meta_dt)) {
  stop("`cell_id` not found in metadata.")
}
if (!"epi_cluster_id" %in% colnames(meta_dt)) {
  stop("`epi_cluster_id` not found in metadata.")
}

meta_dt[, cell_id := as.character(cell_id)]
meta_dt[, epi_cluster_id := as.character(epi_cluster_id)]
meta_dt <- meta_dt[match(colnames(counts_mat), cell_id)]

progression_clusters <- c("E10", "E11", "E13", "E17", "E18", "E19")
meta_dt[, progression_cluster_flag := epi_cluster_id %in% progression_clusters]

if (!all(meta_dt$progression_cluster_flag)) {
  keep_idx <- which(meta_dt$progression_cluster_flag)
  meta_dt <- meta_dt[keep_idx]
  counts_mat <- counts_mat[, keep_idx, drop = FALSE]
}

rownames(counts_mat) <- normalize_gene(rownames(counts_mat))
counts_mat <- counts_mat[!duplicated(rownames(counts_mat)), , drop = FALSE]

tf_dt <- fread(tf_path, header = FALSE)
tf_genes <- unique(normalize_gene(tf_dt[[1]]))

global_detect_n <- Matrix::rowSums(counts_mat > 0)
global_detect_frac <- global_detect_n / ncol(counts_mat)

cluster_detect_list <- lapply(progression_clusters, function(cluster_id) {
  cell_idx <- which(meta_dt$epi_cluster_id == cluster_id)
  frac <- if (length(cell_idx) == 0) {
    rep(0, nrow(counts_mat))
  } else {
    Matrix::rowSums(counts_mat[, cell_idx, drop = FALSE] > 0) / length(cell_idx)
  }
  data.table(
    gene = rownames(counts_mat),
    epi_cluster_id = cluster_id,
    detected_fraction_cluster = frac
  )
})
cluster_detect_dt <- rbindlist(cluster_detect_list, fill = TRUE)
cluster_detect_wide <- dcast(
  cluster_detect_dt,
  gene ~ epi_cluster_id,
  value.var = "detected_fraction_cluster"
)

cluster_frac_mat <- as.matrix(cluster_detect_wide[, -1])
rownames(cluster_frac_mat) <- cluster_detect_wide$gene
max_cluster_detect_frac <- apply(cluster_frac_mat, 1, max, na.rm = TRUE)
cluster_ge5_flag <- max_cluster_detect_frac >= 0.05
global_ge1_flag <- global_detect_frac >= 0.01
tf_flag <- rownames(counts_mat) %in% tf_genes

keep_flag <- cluster_ge5_flag | global_ge1_flag | tf_flag
counts_use <- counts_mat[keep_flag, , drop = FALSE]

gene_filter_dt <- data.table(
  gene = rownames(counts_mat),
  detected_cells_global = as.integer(global_detect_n),
  detected_fraction_global = as.numeric(global_detect_frac),
  max_detected_fraction_cluster = as.numeric(max_cluster_detect_frac),
  keep_cluster_ge5 = as.logical(cluster_ge5_flag),
  keep_global_ge1 = as.logical(global_ge1_flag),
  keep_tf = as.logical(tf_flag),
  keep_final = as.logical(keep_flag)
)

gene_filter_dt <- merge(
  gene_filter_dt,
  cluster_detect_wide,
  by = "gene",
  all.x = TRUE
)

gene_meta_out <- gene_filter_dt[keep_final == TRUE]

expr_out_path <- file.path(tables_dir, "01_pyscenic_expression_matrix_cluster_aware.tsv")
meta_out_path <- file.path(tables_dir, "01_pyscenic_cell_metadata_cluster_aware.csv")
gene_out_path <- file.path(tables_dir, "01_pyscenic_gene_metadata_cluster_aware.csv")
filter_out_path <- file.path(tables_dir, "01_pyscenic_cluster_aware_filter_summary.csv")

fwrite(as.data.table(as.matrix(counts_use), keep.rownames = "gene"), expr_out_path, sep = "\t")
fwrite(meta_dt, meta_out_path)
fwrite(gene_meta_out, gene_out_path)

summary_dt <- data.table(
  n_cells = ncol(counts_mat),
  n_genes_input = nrow(counts_mat),
  n_genes_kept = nrow(counts_use),
  n_genes_cluster_ge5 = sum(cluster_ge5_flag),
  n_genes_global_ge1 = sum(global_ge1_flag),
  n_tf_total_in_reference = length(tf_genes),
  n_tf_present_in_input = sum(tf_flag),
  n_tf_kept = sum(tf_flag & keep_flag),
  kept_by_cluster_only = sum(cluster_ge5_flag & !global_ge1_flag & !tf_flag),
  kept_by_global_only = sum(!cluster_ge5_flag & global_ge1_flag & !tf_flag),
  kept_by_tf_only = sum(!cluster_ge5_flag & !global_ge1_flag & tf_flag)
)
fwrite(summary_dt, filter_out_path)

saveRDS(
  list(
    counts = counts_use,
    metadata = meta_dt,
    gene_metadata = gene_meta_out,
    filter_summary = summary_dt
  ),
  file.path(rds_dir, "01b_progression_epithelial_for_regulon_cluster_aware.rds")
)

cat("Prepared cluster-aware pySCENIC input.\n")
cat("Cells:", ncol(counts_use), "\n")
cat("Genes input:", nrow(counts_mat), "\n")
cat("Genes kept:", nrow(counts_use), "\n")
cat("Expression matrix:", expr_out_path, "\n")
cat("Cell metadata:", meta_out_path, "\n")
cat("Gene metadata:", gene_out_path, "\n")
cat("Filter summary:", filter_out_path, "\n")

project_root <- "/Users/huahaijiao/Desktop/02_study/0206_artical/Hp_gastric_cancer"
result_root <- file.path(
  project_root,
  "gastric_cancer", "results", "Figure4_malignant_program_mechanism"
)
tables_dir <- file.path(result_root, "tables")
rds_dir <- file.path(result_root, "rds")
logs_dir <- file.path(result_root, "logs")
copykat_raw_dir <- file.path(result_root, "copykat_raw")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(copykat_raw_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(data.table)
  library(Matrix)
})

if (!requireNamespace("copykat", quietly = TRUE)) {
  stop(
    "Package 'copykat' is not installed. ",
    "Please install copykat first, then rerun this script."
  )
}

suppressPackageStartupMessages({
  library(copykat)
})

reference_cap_per_stratum <- 20L
min_detected_cells_per_gene <- 20L

normalize_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

calc_cnv_burden <- function(cnv_mat) {
  cnv_mat <- as.matrix(cnv_mat)
  cnv_mat[!is.finite(cnv_mat)] <- NA_real_
  apply(cnv_mat, 2, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) {
      return(NA_real_)
    }
    mean(abs(x), na.rm = TRUE)
  })
}

extract_prediction_table <- function(copykat_obj) {
  pred <- copykat_obj$prediction
  if (is.null(pred)) {
    stop("copykat output does not contain $prediction")
  }
  pred <- as.data.table(pred, keep.rownames = "cell_id")

  cn_cols <- colnames(pred)
  if (!"copykat.pred" %in% cn_cols) {
    alt_pred <- cn_cols[grepl("pred", cn_cols, ignore.case = TRUE)]
    if (length(alt_pred) > 0) {
      setnames(pred, alt_pred[1], "copykat.pred")
    }
  }
  if (!"copykat.pred" %in% colnames(pred)) {
    stop("Could not identify prediction column in copykat output")
  }

  if (!"cell_id" %in% colnames(pred)) {
    stop("Could not identify cell_id column in copykat prediction table")
  }

  pred
}

extract_cnv_matrix <- function(copykat_obj) {
  for (nm in c("CNAmat", "cnv", "copykat.test", "test.cna")) {
    if (!is.null(copykat_obj[[nm]])) {
      return(copykat_obj[[nm]])
    }
  }
  stop("Could not identify CNV matrix in copykat output object")
}

sce_path <- file.path(rds_dir, "01_epithelial_for_cnv.rds")
meta_path <- file.path(tables_dir, "01_epithelial_cells_for_cnv_metadata.csv")

required_files <- c(sce_path, meta_path)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

sce <- readRDS(sce_path)
meta <- fread(meta_path)
meta[, cell_id := as.character(cell_id)]

counts_mat <- counts(sce)
rownames(counts_mat) <- normalize_gene(rownames(counts_mat))
counts_mat <- counts_mat[!duplicated(rownames(counts_mat)), , drop = FALSE]

# copykat prefers raw count matrix and often runs more stably after removing
# extremely low-detection genes.
gene_detect_n <- Matrix::rowSums(counts_mat > 0)
keep_gene <- gene_detect_n >= min_detected_cells_per_gene
counts_use <- counts_mat[keep_gene, , drop = FALSE]

reference_cells <- meta[cnv_reference_flag == TRUE, cell_id]
reference_cells <- intersect(reference_cells, colnames(counts_use))

progression_cells <- meta[progression_cluster_flag == TRUE, cell_id]
progression_cells <- intersect(progression_cells, colnames(counts_use))

# First-pass CNV run: keep all progression-related epithelial cells, sharply
# downsample the NAG/CAG reference pool across sample-cluster strata, and drop
# other epithelial cells so we can obtain a usable cell-level / cluster-level
# result as quickly as possible on the current machine.
set.seed(20260503)

reference_dt <- meta[cnv_reference_flag == TRUE]
reference_dt[, stratum := paste(sample_id, epi_cluster_id, sep = "__")]
reference_cells_use <- unlist(lapply(split(reference_dt$cell_id, reference_dt$stratum), function(x) {
  n_take <- min(length(x), reference_cap_per_stratum)
  sample(x, n_take)
}), use.names = FALSE)

progression_dt <- meta[progression_cluster_flag == TRUE]
progression_cells_use <- progression_dt$cell_id

other_cells <- setdiff(colnames(counts_use), union(reference_cells, progression_cells))
other_cells_use <- character(0)

cells_use <- unique(c(reference_cells_use, progression_cells_use, other_cells_use))
counts_use <- counts_use[, cells_use, drop = FALSE]
meta <- meta[cell_id %in% cells_use]
reference_cells <- intersect(reference_cells_use, colnames(counts_use))
progression_cells <- intersect(progression_cells_use, colnames(counts_use))

# copykat internally handles base matrices more reliably than sparse triplet
# matrices on this machine/runtime.
counts_use <- as.matrix(counts_use)

if (length(reference_cells) < 2000) {
  warning("Reference cell pool is smaller than expected: ", length(reference_cells))
}

copykat_input_summary <- data.table(
  min_detected_cells_per_gene = min_detected_cells_per_gene,
  total_genes = nrow(counts_use),
  total_cells = ncol(counts_use),
  reference_cells = length(reference_cells),
  progression_cells = length(progression_cells),
  downsampled_other_cells = length(other_cells_use)
)
fwrite(copykat_input_summary, file.path(tables_dir, "02_copykat_input_summary.csv"))

set.seed(20260503)
old_wd <- getwd()
setwd(copykat_raw_dir)
on.exit(setwd(old_wd), add = TRUE)
copykat_obj <- copykat::copykat(
  rawmat = counts_use,
  id.type = "S",
  ngene.chr = 5,
  win.size = 25,
  LOW.DR = 0.1,
  sam.name = "GSE134520_epithelial",
  distance = "euclidean",
  norm.cell.names = reference_cells,
  output.seg = FALSE,
  plot.genes = FALSE,
  genome = "hg20",
  n.cores = 1
)

saveRDS(copykat_obj, file.path(rds_dir, "02_copykat_result.rds"))

pred_dt <- extract_prediction_table(copykat_obj)
pred_dt[, cell_id := as.character(cell_id)]
pred_dt <- pred_dt[!duplicated(cell_id)]

cnv_mat <- extract_cnv_matrix(copykat_obj)
cnv_burden <- calc_cnv_burden(cnv_mat)
cnv_burden_dt <- data.table(
  cell_id = colnames(cnv_mat),
  cnv_burden = as.numeric(cnv_burden)
)

copykat_meta <- merge(meta, pred_dt, by = "cell_id", all.x = TRUE)
copykat_meta <- merge(copykat_meta, cnv_burden_dt, by = "cell_id", all.x = TRUE)

pred_col <- "copykat.pred"
if (!pred_col %in% colnames(copykat_meta)) {
  stop("copykat prediction column missing after merge")
}

copykat_meta[, copykat_pred_simple := fifelse(
  grepl("aneuploid", get(pred_col), ignore.case = TRUE), "aneuploid",
  fifelse(grepl("diploid", get(pred_col), ignore.case = TRUE), "diploid", as.character(get(pred_col)))
)]

cluster_summary <- copykat_meta[, .(
  n_cells = .N,
  n_predicted = sum(!is.na(copykat_pred_simple)),
  aneuploid_fraction = mean(copykat_pred_simple == "aneuploid", na.rm = TRUE),
  diploid_fraction = mean(copykat_pred_simple == "diploid", na.rm = TRUE),
  cnv_burden_median = median(cnv_burden, na.rm = TRUE),
  cnv_burden_mean = mean(cnv_burden, na.rm = TRUE),
  progression_cluster = unique(progression_cluster_flag)[1]
), by = .(epi_cluster_id)][order(epi_cluster_id)]

group_summary <- copykat_meta[, .(
  n_cells = .N,
  n_predicted = sum(!is.na(copykat_pred_simple)),
  aneuploid_fraction = mean(copykat_pred_simple == "aneuploid", na.rm = TRUE),
  cnv_burden_median = median(cnv_burden, na.rm = TRUE)
), by = .(lesion_group)][order(lesion_group)]

progression_vs_other <- copykat_meta[, .(
  n_cells = .N,
  n_predicted = sum(!is.na(copykat_pred_simple)),
  aneuploid_fraction = mean(copykat_pred_simple == "aneuploid", na.rm = TRUE),
  cnv_burden_median = median(cnv_burden, na.rm = TRUE)
), by = .(progression_cluster_flag)]

fwrite(copykat_meta, file.path(tables_dir, "02_copykat_cell_level_results.csv"))
fwrite(cluster_summary, file.path(tables_dir, "02_copykat_cluster_summary.csv"))
fwrite(group_summary, file.path(tables_dir, "02_copykat_group_summary.csv"))
fwrite(progression_vs_other, file.path(tables_dir, "02_copykat_progression_vs_other_summary.csv"))

cat("copykat analysis completed.\n")
cat("Genes used:", nrow(counts_use), "\n")
cat("Cells used:", ncol(counts_use), "\n")
cat("Reference cells:", length(reference_cells), "\n")
cat("Progression cells:", length(progression_cells), "\n")
cat("Saved RDS:", file.path(rds_dir, "02_copykat_result.rds"), "\n")
cat("Saved cell-level table:", file.path(tables_dir, "02_copykat_cell_level_results.csv"), "\n")

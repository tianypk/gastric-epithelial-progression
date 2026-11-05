args <- commandArgs(trailingOnly = TRUE)
project_root <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  Sys.getenv("HP_GC_PROJECT_ROOT", unset = "/Users/huahaijiao/Desktop/02_study/0206_artical/Hp_gastric_cancer")
}

result_root <- file.path(
  project_root,
  "gastric_cancer", "results", "Figure6_immunotherapy_validation"
)
figure5_root <- file.path(
  project_root,
  "gastric_cancer", "results", "Figure5_regulon_pathway_bulk"
)
pyscenic_run_tag <- "pyscenic_out_cluster_aware_server_0121_20261001"
pyscenic_prefix <- "gastric_progression_cluster_aware_0121"
progression_root <- file.path(
  project_root,
  "gastric_cancer", "results", "progression_score_main47"
)
input_root <- file.path(project_root, "immunotherapy_datasets", "PRJEB25780")
tables_dir <- file.path(result_root, "tables")
rds_dir <- file.path(result_root, "rds")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

read_tsv_with_optional_index <- function(path, key_col) {
  header_line <- readLines(path, n = 1, warn = FALSE)
  header_fields <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  dt <- fread(path, skip = 1, header = FALSE, fill = TRUE)
  if (ncol(dt) == length(header_fields) + 1) {
    setnames(dt, c(".row_index", header_fields))
  } else if (ncol(dt) == length(header_fields)) {
    setnames(dt, header_fields)
  } else {
    stop("Unexpected column count while reading: ", basename(path))
  }
  if (".row_index" %in% colnames(dt)) {
    dt[, .row_index := NULL]
  }
  if (!key_col %in% colnames(dt)) {
    stop("Key column not found after parsing: ", key_col)
  }
  dt
}

normalize_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

scale_by_gene <- function(expr_mat) {
  expr_mat <- as.matrix(expr_mat)
  z <- t(scale(t(expr_mat)))
  z[!is.finite(z)] <- NA_real_
  z
}

weighted_signature_score <- function(expr_z, gene_weights) {
  genes_use <- intersect(names(gene_weights), rownames(expr_z))
  if (length(genes_use) == 0) {
    return(list(score = rep(NA_real_, ncol(expr_z)), n_used = 0L))
  }
  w <- gene_weights[genes_use]
  score <- colSums(expr_z[genes_use, , drop = FALSE] * w, na.rm = TRUE) / sum(abs(w))
  list(score = score, n_used = length(genes_use))
}

calc_module_score <- function(expr_z, genes) {
  genes_use <- intersect(normalize_gene(genes), rownames(expr_z))
  if (length(genes_use) == 0) {
    return(list(score = rep(NA_real_, ncol(expr_z)), n_used = 0L))
  }
  list(
    score = colMeans(expr_z[genes_use, , drop = FALSE], na.rm = TRUE),
    n_used = length(genes_use)
  )
}

extract_target_genes <- function(x) {
  x <- as.character(x)
  if (!nzchar(x) || is.na(x)) {
    return(character(0))
  }
  m <- gregexpr("'([^']+)'", x, perl = TRUE)
  genes <- regmatches(x, m)[[1]]
  genes <- gsub("^'|'$", "", genes)
  unique(genes)
}

clean_regulon_table <- function(path) {
  raw <- fread(path, fill = TRUE, header = FALSE)
  header1 <- as.character(unlist(raw[1]))
  header2 <- as.character(unlist(raw[2]))
  header3 <- as.character(unlist(raw[3]))
  nms <- header3
  for (i in seq_along(nms)) {
    if (is.na(nms[i]) || nms[i] == "") {
      part2 <- header2[i]
      if (is.na(part2) || part2 == "") {
        part2 <- header1[i]
      }
      if (is.na(part2) || part2 == "") {
        part2 <- paste0("V", i)
      }
      nms[i] <- part2
    }
  }
  dt <- raw[-c(1:3)]
  setnames(dt, nms)
  dt[]
}

required_files <- c(
  file.path(tables_dir, "01_prjeb25780_cohort_manifest.csv"),
  file.path(input_root, "STAD-PRJEB25780_anti-PD-1.Response.expression.tsv"),
  file.path(progression_root, "tables", "progression_score_47g_definition.csv"),
  file.path(figure5_root, pyscenic_run_tag, paste0(pyscenic_prefix, ".reg.tsv"))
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

manifest_dt <- fread(file.path(tables_dir, "01_prjeb25780_cohort_manifest.csv"), fill = TRUE)
expr_dt <- read_tsv_with_optional_index(
  file.path(input_root, "STAD-PRJEB25780_anti-PD-1.Response.expression.tsv"),
  "GENE_SYMBOL"
)
expr_mat <- as.matrix(expr_dt[, -1, with = FALSE])
rownames(expr_mat) <- normalize_gene(expr_dt[["GENE_SYMBOL"]])
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
storage.mode(expr_mat) <- "numeric"
expr_z <- scale_by_gene(expr_mat)

progression_genes_dt <- fread(file.path(
  progression_root, "tables", "progression_score_47g_definition.csv"
))
progression_weights <- ifelse(progression_genes_dt$direction == "up_gain", 1, -1)
names(progression_weights) <- normalize_gene(progression_genes_dt$gene)
progression_weights <- progression_weights[!is.na(names(progression_weights))]

regulon_reg_dt <- clean_regulon_table(file.path(
  figure5_root, pyscenic_run_tag, paste0(pyscenic_prefix, ".reg.tsv")
))
key_regulons <- c("MAX(+)", "NR1H4(+)", "HNF4G(+)")
regulon_targets <- lapply(gsub("\\(\\+\\)$", "", key_regulons), function(tf) {
  sub <- regulon_reg_dt[TF == tf]
  unique(sort(unlist(lapply(sub$TargetGenes, extract_target_genes))))
})
names(regulon_targets) <- key_regulons

progression_score_obj <- weighted_signature_score(expr_z, progression_weights)

score_dt <- data.table(sample_id = colnames(expr_z))
score_dt[, progression_score_47g := progression_score_obj$score]
score_dt[, progression_score := progression_score_47g]
score_dt[, n_progression_genes_used := progression_score_obj$n_used]

gene_usage_dt <- data.table(
  score_name = "progression_score",
  gene_set_source = "Figure1_core_progression_signature",
  genes_in_signature = length(progression_weights),
  genes_used = progression_score_obj$n_used
)

for (reg_name in names(regulon_targets)) {
  mod_obj <- calc_module_score(expr_z, regulon_targets[[reg_name]])
  score_col <- sub("\\(\\+\\)", "", reg_name)
  score_dt[, (paste0(score_col, "_module_score")) := mod_obj$score]
  score_dt[, (paste0("n_", score_col, "_targets_used")) := mod_obj$n_used]
  gene_usage_dt <- rbind(
    gene_usage_dt,
    data.table(
      score_name = paste0(score_col, "_module_score"),
      gene_set_source = paste0("Figure5_", reg_name, "_targets"),
      genes_in_signature = length(regulon_targets[[reg_name]]),
      genes_used = mod_obj$n_used
    ),
    fill = TRUE
  )
}

score_dt <- merge(
  manifest_dt,
  score_dt,
  by = "sample_id",
  all.x = TRUE
)

scored_bundle <- list(
  scored_data = score_dt,
  gene_usage = gene_usage_dt,
  key_regulons = key_regulons,
  regulon_targets = regulon_targets
)

fwrite(score_dt, file.path(tables_dir, "02_prjeb25780_projected_scores.csv"))
fwrite(gene_usage_dt, file.path(tables_dir, "02_prjeb25780_score_gene_usage_summary.csv"))
saveRDS(scored_bundle, file.path(rds_dir, "02_prjeb25780_scored_dataset.rds"))

cat("Figure 6 step 02 completed.\n")
cat("Scored samples:", nrow(score_dt), "\n")

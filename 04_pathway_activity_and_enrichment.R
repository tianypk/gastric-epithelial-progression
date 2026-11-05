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
pyscenic_run_tag <- "pyscenic_out_cluster_aware_server_0121_20261001"
pyscenic_prefix <- "gastric_progression_cluster_aware_0121"
pyscenic_dir <- file.path(result_root, pyscenic_run_tag)
tables_dir <- file.path(result_root, "tables")
figures_dir <- file.path(result_root, "figures")
rds_dir <- file.path(result_root, "rds")
supp_dir <- file.path(result_root, "supplementary_tables")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(Matrix)
  library(fgsea)
})

has_gsva <- requireNamespace("GSVA", quietly = TRUE)

cluster_order <- c("E10", "E11", "E13", "E17", "E18", "E19")
key_regulons <- c("MAX(+)", "NR1H4(+)", "HNF4G(+)")
hallmark_gmt_path <- file.path(
  result_root, "reference", "h.all.v2025.1.Hs.symbols.gmt"
)
hallmark_keep <- c(
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_KRAS_SIGNALING_UP",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_BILE_ACID_METABOLISM",
  "HALLMARK_XENOBIOTIC_METABOLISM",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
)

read_hallmark_gmt <- function(path) {
  gmt <- readLines(path, warn = FALSE)
  gmt <- gmt[nzchar(gmt)]
  split_sets <- lapply(gmt, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    genes <- unique(fields[-c(1, 2)])
    genes[nzchar(genes)]
  })
  names(split_sets) <- vapply(
    strsplit(gmt, "\t", fixed = TRUE),
    `[`,
    character(1),
    1
  )
  split_sets
}

zscore_rows <- function(mat) {
  row_means <- rowMeans(mat, na.rm = TRUE)
  row_sds <- apply(mat, 1, sd, na.rm = TRUE)
  row_sds[row_sds == 0 | is.na(row_sds)] <- 1
  z <- sweep(mat, 1, row_means, "-")
  z <- sweep(z, 1, row_sds, "/")
  z[z > 2] <- 2
  z[z < -2] <- -2
  z
}

pretty_hallmark <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  gsub("_", " ", x)
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
  if (nrow(raw) < 4) {
    stop("Regulon file appears too short: ", path)
  }
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

auc_rds_path <- file.path(rds_dir, "02_regulon_auc_matrix.rds")
expr_rds_path <- file.path(rds_dir, "01b_progression_epithelial_for_regulon_cluster_aware.rds")
cell_meta_path <- file.path(tables_dir, "01_pyscenic_cell_metadata_cluster_aware.csv")
regulon_cor_path <- file.path(tables_dir, "02_regulon_correlation_with_cnv_program.csv")
regulon_reg_path <- file.path(pyscenic_dir, paste0(pyscenic_prefix, ".reg.tsv"))

required_files <- c(
  auc_rds_path, expr_rds_path, cell_meta_path,
  regulon_cor_path, regulon_reg_path
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

cat("Reading progression expression object...\n")
expr_obj <- readRDS(expr_rds_path)
counts_mat <- expr_obj$counts
cell_meta <- fread(cell_meta_path)
cell_meta <- cell_meta[regulon_group %in% cluster_order]
cell_meta[, regulon_group := factor(as.character(regulon_group), levels = cluster_order)]

common_cells <- intersect(colnames(counts_mat), cell_meta$cell_id)
if (length(common_cells) == 0) {
  stop("No overlapping cells between counts matrix and cell metadata")
}
counts_mat <- counts_mat[, common_cells, drop = FALSE]
cell_meta <- cell_meta[match(common_cells, cell_id)]

log_expr_mat <- log1p(counts_mat)

cat("Preparing Hallmark pathways...\n")
if (file.exists(hallmark_gmt_path)) {
  hallmark_list <- read_hallmark_gmt(hallmark_gmt_path)
  hallmark_list <- hallmark_list[names(hallmark_list) %in% hallmark_keep]
  hallmark_dt <- rbindlist(lapply(names(hallmark_list), function(gs_name) {
    data.table(gs_name = gs_name, gene_symbol = hallmark_list[[gs_name]])
  }))
} else {
  suppressPackageStartupMessages(library(msigdbr))
  hallmark_dt <- as.data.table(msigdbr(species = "Homo sapiens", category = "H"))
  hallmark_dt <- hallmark_dt[gs_name %in% hallmark_keep]
  hallmark_list <- split(hallmark_dt$gene_symbol, hallmark_dt$gs_name)
  hallmark_list <- lapply(hallmark_list, unique)
}

present_clusters <- cluster_order[cluster_order %in% unique(as.character(cell_meta$regulon_group))]
cluster_expr_list <- lapply(present_clusters, function(cl) {
  idx <- cell_meta$regulon_group == cl
  Matrix::rowMeans(log_expr_mat[, idx, drop = FALSE])
})
names(cluster_expr_list) <- present_clusters
cluster_expr_mat <- do.call(cbind, cluster_expr_list)
storage.mode(cluster_expr_mat) <- "numeric"
rownames(cluster_expr_mat) <- rownames(log_expr_mat)

if (has_gsva) {
  cat("Scoring pathway activity with GSVA on cluster-level expression...\n")
  gsva_param <- GSVA::gsvaParam(cluster_expr_mat, hallmark_list, kcdf = "Gaussian")
  pathway_cluster_mat <- GSVA::gsva(gsva_param, verbose = FALSE)
  pathway_scoring_method <- "GSVA"
} else {
  cat("GSVA is unavailable; using hallmark gene-set mean row-z-score fallback...\n")
  pathway_cluster_mat <- do.call(rbind, lapply(names(hallmark_list), function(gs_name) {
    genes <- intersect(hallmark_list[[gs_name]], rownames(cluster_expr_mat))
    if (length(genes) == 0) {
      return(rep(NA_real_, length(present_clusters)))
    }
    gene_z <- zscore_rows(cluster_expr_mat[genes, present_clusters, drop = FALSE])
    colMeans(gene_z, na.rm = TRUE)
  }))
  rownames(pathway_cluster_mat) <- names(hallmark_list)
  colnames(pathway_cluster_mat) <- present_clusters
  pathway_scoring_method <- "hallmark_mean_row_zscore_fallback"
}

pathway_mean_mat <- pathway_cluster_mat
if (is.vector(pathway_mean_mat)) {
  pathway_mean_mat <- matrix(
    pathway_mean_mat,
    ncol = 1,
    dimnames = list(rownames(pathway_cluster_mat), present_clusters)
  )
}

pathway_dt <- as.data.table(pathway_mean_mat, keep.rownames = "hallmark_pathway")
pathway_dt[, pathway_label := pretty_hallmark(hallmark_pathway)]
setcolorder(pathway_dt, c("hallmark_pathway", "pathway_label", present_clusters))
fwrite(pathway_dt, file.path(tables_dir, "04_pathway_activity_by_cluster_mean.csv"))

pathway_long <- melt(
  pathway_dt,
  id.vars = c("hallmark_pathway", "pathway_label"),
  variable.name = "cluster",
  value.name = "mean_activity"
)
pathway_long[, cluster := factor(cluster, levels = present_clusters)]
pathway_long[, mean_diff_vs_others := {
  cur <- cluster[1]
  mean_activity - mean(pathway_long[hallmark_pathway == .BY$hallmark_pathway & cluster != cur, mean_activity], na.rm = TRUE)
}, by = .(hallmark_pathway, cluster)]

pathway_focus <- copy(pathway_long)[cluster %in% c("E13", "E17", "E18", "E19")]
pathway_focus_summary <- pathway_focus[, .(
  max_mean_activity = max(mean_activity, na.rm = TRUE),
  max_mean_diff = max(mean_diff_vs_others, na.rm = TRUE)
), by = .(hallmark_pathway, pathway_label)]
setorder(pathway_focus_summary, -max_mean_diff, -max_mean_activity)
selected_pathways <- unique(c(
  pathway_focus_summary$hallmark_pathway[1:min(10, nrow(pathway_focus_summary))],
  hallmark_keep
))
selected_pathways <- selected_pathways[selected_pathways %in% rownames(pathway_mean_mat)]

pathway_heatmap_mat <- pathway_mean_mat[selected_pathways, present_clusters, drop = FALSE]
pathway_heatmap_z <- zscore_rows(pathway_heatmap_mat)
heatmap_dt <- as.data.table(pathway_heatmap_z, keep.rownames = "hallmark_pathway")
heatmap_dt <- merge(
  heatmap_dt,
  unique(pathway_dt[, .(hallmark_pathway, pathway_label)]),
  by = "hallmark_pathway",
  all.x = TRUE
)
heatmap_long <- melt(
  heatmap_dt,
  id.vars = c("hallmark_pathway", "pathway_label"),
  variable.name = "cluster",
  value.name = "z_score"
)
row_order <- rev(unique(heatmap_dt$pathway_label))
heatmap_long[, pathway_label := factor(pathway_label, levels = row_order)]
heatmap_long[, cluster := factor(cluster, levels = present_clusters)]

pathway_plot <- ggplot(heatmap_long, aes(x = cluster, y = pathway_label, fill = z_score)) +
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
    title = "E  Pathway activity across progression-related epithelial states"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave(
  file.path(figures_dir, "Figure5_pathway_summary_main.pdf"),
  pathway_plot,
  width = 7.6,
  height = max(5.5, 0.32 * length(row_order) + 1.8),
  dpi = 300
)
ggsave(
  file.path(figures_dir, "Figure5_pathway_summary_main.png"),
  pathway_plot,
  width = 7.6,
  height = max(5.5, 0.32 * length(row_order) + 1.8),
  dpi = 300
)

cat("Reading key regulon annotations...\n")
regulon_cor_dt <- fread(regulon_cor_path)
regulon_reg_dt <- clean_regulon_table(regulon_reg_path)
setnames(regulon_reg_dt, old = c("TF", "TargetGenes"), new = c("TF", "TargetGenes"), skip_absent = TRUE)
regulon_reg_dt <- regulon_reg_dt[TF %in% gsub("\\(\\+\\)$", "", key_regulons)]
regulon_reg_dt[, regulon := paste0(TF, "(+)")]
regulon_reg_dt[, target_gene_vector := lapply(TargetGenes, extract_target_genes)]
regulon_reg_dt[, target_gene_n := lengths(target_gene_vector)]
regulon_reg_dt <- regulon_reg_dt[target_gene_n >= 5]

hallmark_simple <- split(hallmark_dt$gene_symbol, hallmark_dt$gs_name)
hallmark_simple <- lapply(hallmark_simple, unique)
all_genes <- unique(rownames(log_expr_mat))

cat("Running fgsea on key regulon target genes...\n")
enrichment_list <- lapply(seq_len(nrow(regulon_reg_dt)), function(i) {
  regulon_name <- regulon_reg_dt$regulon[i]
  targets <- intersect(regulon_reg_dt$target_gene_vector[[i]], all_genes)
  if (length(targets) < 5) {
    return(NULL)
  }

  ranks <- setNames(ifelse(all_genes %in% targets, 1, 0), all_genes)
  fg <- suppressWarnings(
    fgsea(
      pathways = hallmark_simple,
      stats = ranks,
      minSize = 5,
      maxSize = 500,
      nperm = 1000
    )
  )
  if (nrow(fg) == 0) {
    return(NULL)
  }
  fg <- as.data.table(fg)
  fg[, regulon := regulon_name]
  fg[, target_gene_n := length(targets)]
  fg[]
})
enrichment_dt <- rbindlist(enrichment_list, fill = TRUE)
if (nrow(enrichment_dt) == 0) {
  warning("No regulon target enrichment results were generated.")
  enrichment_dt <- data.table()
}

if (nrow(enrichment_dt) > 0) {
  enrichment_dt[, pathway_label := pretty_hallmark(pathway)]
  enrichment_dt <- merge(
    enrichment_dt,
    unique(regulon_cor_dt[, .(
      regulon, dominant_cluster, dominant_mean_auc,
      rho_raw_cnv_burden, fdr_raw_cnv_burden,
      rho_malignant_program_score, fdr_malignant_program_score,
      shared_significant_positive
    )]),
    by = "regulon",
    all.x = TRUE
  )
  setorder(enrichment_dt, regulon, padj, -NES)
  fwrite(enrichment_dt, file.path(tables_dir, "04_regulon_target_hallmark_enrichment.csv"))

  top_enrichment_dt <- enrichment_dt[padj < 0.05][order(regulon, padj, -NES), head(.SD, 8), by = regulon]
  fwrite(top_enrichment_dt, file.path(tables_dir, "04_regulon_target_hallmark_enrichment_top_hits.csv"))
} else {
  fwrite(enrichment_dt, file.path(tables_dir, "04_regulon_target_hallmark_enrichment.csv"))
}

saveRDS(
  list(
    pathway_cluster_matrix = pathway_cluster_mat,
    pathway_cluster_mean = pathway_mean_mat,
    pathway_activity_table = pathway_dt,
    key_regulon_target_enrichment = enrichment_dt,
    pathway_scoring_method = pathway_scoring_method
  ),
  file.path(rds_dir, "03_figure5_integrated_object.rds")
)

cat("Figure 5 pathway activity and enrichment completed.\n")
cat("Main pathway heatmap:", file.path(figures_dir, "Figure5_pathway_summary_main.pdf"), "\n")
cat("Pathway table:", file.path(tables_dir, "04_pathway_activity_by_cluster_mean.csv"), "\n")
cat("Regulon target enrichment table:", file.path(tables_dir, "04_regulon_target_hallmark_enrichment.csv"), "\n")

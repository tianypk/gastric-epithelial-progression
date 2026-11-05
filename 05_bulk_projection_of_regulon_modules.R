project_root <- "/Users/huahaijiao/Desktop/02_study/0206_artical/Hp_gastric_cancer"
result_root <- file.path(
  project_root,
  "gastric_cancer", "results", "Figure5_regulon_pathway_bulk"
)
tables_dir <- file.path(result_root, "tables")
figures_dir <- file.path(result_root, "figures")
manuscript_dir <- file.path(result_root, "manuscript")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(survival)
})

normalize_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

scale_by_gene <- function(expr_mat) {
  expr_mat <- as.matrix(expr_mat)
  z <- t(scale(t(expr_mat)))
  z[!is.finite(z)] <- NA_real_
  z
}

calc_module_score <- function(expr_mat, genes) {
  expr_mat <- as.matrix(expr_mat)
  rownames(expr_mat) <- normalize_gene(rownames(expr_mat))
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
  expr_z <- scale_by_gene(expr_mat)
  genes_use <- intersect(normalize_gene(genes), rownames(expr_z))
  if (length(genes_use) == 0) {
    return(data.table(
      ID = colnames(expr_z),
      module_score = NA_real_,
      n_genes_used = 0L
    ))
  }
  data.table(
    ID = colnames(expr_z),
    module_score = colMeans(expr_z[genes_use, , drop = FALSE], na.rm = TRUE),
    n_genes_used = length(genes_use)
  )
}

make_km_group <- function(x) {
  med <- median(x, na.rm = TRUE)
  fifelse(x >= med, "High", "Low")
}

fit_cox_summary <- function(dt, score_col = "module_score") {
  sub <- dt[is.finite(OS_time) & !is.na(OS_event) & !is.na(get(score_col))]
  if (nrow(sub) < 30 || length(unique(sub$OS_event)) < 2) {
    return(data.table(
      n = nrow(sub), events = sum(sub$OS_event, na.rm = TRUE),
      hr = NA_real_, lower95 = NA_real_, upper95 = NA_real_, pvalue = NA_real_
    ))
  }
  sub[, score_z := as.numeric(scale(get(score_col)))]
  fit <- coxph(Surv(OS_time, OS_event) ~ score_z, data = sub)
  sm <- summary(fit)
  data.table(
    n = nrow(sub),
    events = sum(sub$OS_event, na.rm = TRUE),
    hr = unname(sm$coefficients[1, "exp(coef)"]),
    lower95 = unname(sm$conf.int[1, "lower .95"]),
    upper95 = unname(sm$conf.int[1, "upper .95"]),
    pvalue = unname(sm$coefficients[1, "Pr(>|z|)"])
  )
}

safe_kruskal <- function(score, group) {
  ok <- is.finite(score) & !is.na(group) & group != ""
  score <- score[ok]
  group <- group[ok]
  if (length(score) < 10 || length(unique(group)) < 2) {
    return(NA_real_)
  }
  tryCatch(kruskal.test(score ~ group)$p.value, error = function(e) NA_real_)
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 10) {
    return(data.table(rho = NA_real_, pvalue = NA_real_, n = length(x)))
  }
  test <- tryCatch(
    suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(test)) {
    return(data.table(rho = NA_real_, pvalue = NA_real_, n = length(x)))
  }
  data.table(rho = unname(test$estimate), pvalue = test$p.value, n = length(x))
}

make_km_plot <- function(dt, dataset_label, score_label) {
  sub <- dt[is.finite(OS_time) & !is.na(OS_event) & !is.na(module_score)]
  if (nrow(sub) < 10 || length(unique(sub$OS_event)) < 2) {
    return(
      ggplot() +
        theme_void(base_size = 11) +
        labs(title = paste0(dataset_label, " - ", score_label, "\nInsufficient survival data"))
    )
  }
  sub[, score_group := make_km_group(module_score)]
  fit <- survfit(Surv(OS_time, OS_event) ~ score_group, data = sub)
  if (requireNamespace("survminer", quietly = TRUE)) {
    return(
      survminer::ggsurvplot(
        fit,
        data = sub,
        risk.table = FALSE,
        pval = TRUE,
        conf.int = FALSE,
        palette = c("#4F6FAE", "#BB4B39"),
        title = paste0(dataset_label, " - ", score_label),
        legend.title = score_label,
        legend.labs = c("High", "Low"),
        xlab = "Months",
        ylab = "Overall survival probability",
        ggtheme = theme_bw(base_size = 11) + theme(panel.grid = element_blank())
      )$plot
    )
  }

  fit_summary <- summary(fit)
  curve_dt <- data.table(
    time = fit_summary$time,
    survival = fit_summary$surv,
    strata = if (is.null(fit_summary$strata)) {
      "All"
    } else {
      sub("^score_group=", "", as.character(fit_summary$strata))
    }
  )
  logrank <- survdiff(Surv(OS_time, OS_event) ~ score_group, data = sub)
  logrank_p <- 1 - pchisq(logrank$chisq, df = max(1, length(logrank$n) - 1))

  ggplot(curve_dt, aes(x = time, y = survival, color = strata)) +
    geom_step(linewidth = 0.7) +
    scale_color_manual(values = c(High = "#4F6FAE", Low = "#BB4B39"), drop = FALSE) +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw(base_size = 11) +
    theme(panel.grid = element_blank()) +
    labs(
      title = paste0(
        dataset_label, " - ", score_label,
        " (log-rank p=", format.pval(logrank_p, digits = 2, eps = 1e-3), ")"
      ),
      x = "Months",
      y = "Overall survival probability",
      color = score_label
    )
}

make_boxplot <- function(dt, xvar, title_text) {
  sub <- dt[!is.na(get(xvar)) & get(xvar) != ""]
  ggplot(sub, aes(x = get(xvar), y = module_score, fill = get(xvar))) +
    geom_boxplot(outlier.size = 0.25, width = 0.72) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.text.x = element_text(angle = 30, hjust = 1)
    ) +
    labs(x = NULL, y = "Module score", title = title_text)
}

read_expr_rdata <- function(fp, object_name = NULL) {
  env <- new.env(parent = emptyenv())
  load(fp, envir = env)
  if (is.null(object_name)) {
    object_name <- ls(env)[1]
  }
  as.matrix(get(object_name, envir = env))
}

read_sample_table <- function(dataset_name) {
  file_dataset <- if (dataset_name == "TCGA") "TCGA_STAD" else dataset_name
  dt <- fread(file.path(
    project_root,
    "gastric_cancer", "results", "progression_score_main47", "tables",
    paste0("progression_score_47g_", file_dataset, "_sample_table.csv")
  ))
  dt[, dataset := dataset_name]
  dt
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

pyscenic_run_tag <- "pyscenic_out_cluster_aware_server_0121_20261001"
pyscenic_prefix <- "gastric_progression_cluster_aware_0121"
requested_key_regulons <- c("MAX(+)", "NR1H4(+)", "HNF4G(+)")
key_regulons <- requested_key_regulons
regulon_reg_path <- file.path(result_root, pyscenic_run_tag, paste0(pyscenic_prefix, ".reg.tsv"))
regulon_cor_path <- file.path(result_root, "tables", "02_regulon_correlation_with_cnv_program.csv")

regulon_reg_dt <- clean_regulon_table(regulon_reg_path)
regulon_reg_dt <- regulon_reg_dt[TF %in% gsub("\\(\\+\\)$", "", key_regulons)]
regulon_reg_dt[, regulon := paste0(TF, "(+)")]
regulon_reg_dt[, target_gene_vector := lapply(TargetGenes, extract_target_genes)]
regulon_reg_dt[, target_gene_n := lengths(target_gene_vector)]
regulon_reg_dt <- regulon_reg_dt[target_gene_n >= 5]
regulon_reg_dt <- regulon_reg_dt[, .(
  target_gene_vector = list(sort(unique(unlist(target_gene_vector))))
), by = regulon]
regulon_reg_dt[, target_gene_n := lengths(target_gene_vector)]
key_regulons <- requested_key_regulons[requested_key_regulons %in% regulon_reg_dt$regulon]
missing_key_regulons <- setdiff(requested_key_regulons, key_regulons)
if (length(missing_key_regulons) > 0) {
  warning("Skipping key regulons absent from the cluster-aware pySCENIC result: ",
          paste(missing_key_regulons, collapse = ", "))
}
if (length(key_regulons) == 0) {
  stop("None of the requested key regulons are present in the cluster-aware pySCENIC result.")
}

regulon_cor_dt <- fread(regulon_cor_path)
regulon_meta <- unique(regulon_cor_dt[
  regulon %in% key_regulons,
  .(regulon, dominant_cluster, rho_raw_cnv_burden, rho_malignant_program_score)
])

expr_list <- list(
  TCGA = read_expr_rdata(file.path(project_root, "TCGA_GC_mRNAexp.RData"), "TCGA_GC_mRNAexp"),
  ACRG = read_expr_rdata(file.path(project_root, "GSE62254_mRNAexp.RData"), "GSE62254_mRNAexp"),
  GSE15459 = read_expr_rdata(file.path(project_root, "GSE15459_mRNAexp.RData"), "GSE15459_mRNAexp"),
  GSE84437 = read_expr_rdata(file.path(project_root, "GSE84437_mRNAexp.RData"), "GSE84437_mRNAexp")
)

clin_list <- list(
  TCGA = read_sample_table("TCGA"),
  ACRG = read_sample_table("ACRG"),
  GSE15459 = read_sample_table("GSE15459"),
  GSE84437 = read_sample_table("GSE84437")
)

score_tables <- list()
dataset_summary_list <- list()
cox_list <- list()
assoc_list <- list()

for (dataset_name in names(expr_list)) {
  expr_mat <- expr_list[[dataset_name]]
  clin_dt <- copy(clin_list[[dataset_name]])

  dataset_scores <- rbindlist(lapply(seq_len(nrow(regulon_reg_dt)), function(i) {
    regulon_name <- regulon_reg_dt$regulon[i]
    target_genes <- regulon_reg_dt$target_gene_vector[[i]]
    score_dt <- calc_module_score(expr_mat, target_genes)
    score_dt[, regulon := regulon_name]
    score_dt[, dataset := dataset_name]
    score_dt
  }), fill = TRUE)

  merged_dt <- merge(dataset_scores, clin_dt, by = c("ID", "dataset"), all.x = TRUE)
  merged_dt <- merge(merged_dt, regulon_meta, by = "regulon", all.x = TRUE)
  if (!"stage_simple" %in% colnames(merged_dt)) {
    merged_dt[, stage_simple := NA_character_]
  }
  merged_dt[, module_group := make_km_group(module_score), by = regulon]

  dataset_summary <- merged_dt[, .(
    n_samples = uniqueN(ID),
    mean_module_score = mean(module_score, na.rm = TRUE),
    sd_module_score = sd(module_score, na.rm = TRUE),
    n_genes_used = unique(n_genes_used)[1]
  ), by = .(dataset, regulon, dominant_cluster, rho_raw_cnv_burden, rho_malignant_program_score)]

  cox_dt <- merged_dt[, cbind(
    data.table(dataset = dataset_name, regulon = regulon[1]),
    fit_cox_summary(.SD)
  ), by = regulon]

  assoc_dt <- merged_dt[, {
    lauren_p <- safe_kruskal(module_score, lauren_simple)
    acrg_p <- if ("acrg_subtype" %in% names(.SD)) safe_kruskal(module_score, acrg_subtype) else NA_real_
    stage_rank <- suppressWarnings(as.numeric(factor(stage_simple)))
    stage_cor <- safe_spearman(stage_rank, module_score)
    prog_cor <- safe_spearman(progression_score_47g, module_score)
    data.table(
      lauren_kruskal_p = lauren_p,
      acrg_subtype_kruskal_p = acrg_p,
      stage_spearman_rho = stage_cor$rho,
      stage_spearman_p = stage_cor$pvalue,
      progression_score_spearman_rho = prog_cor$rho,
      progression_score_spearman_p = prog_cor$pvalue,
      progression_score_n = prog_cor$n
    )
  }, by = regulon]
  assoc_dt[, dataset := dataset_name]

  fwrite(
    merged_dt,
    file.path(tables_dir, paste0("05_bulk_projection_", dataset_name, "_sample_table.csv"))
  )

  score_tables[[dataset_name]] <- merged_dt
  dataset_summary_list[[dataset_name]] <- dataset_summary
  cox_list[[dataset_name]] <- cox_dt
  assoc_list[[dataset_name]] <- assoc_dt
}

all_scores_dt <- rbindlist(score_tables, fill = TRUE)
summary_dt <- rbindlist(dataset_summary_list, fill = TRUE)
cox_dt <- rbindlist(cox_list, fill = TRUE)
assoc_dt <- rbindlist(assoc_list, fill = TRUE)

fwrite(all_scores_dt, file.path(tables_dir, "05_bulk_projection_all_samples_combined.csv"))
fwrite(summary_dt, file.path(tables_dir, "05_bulk_projection_dataset_summary.csv"))
fwrite(cox_dt, file.path(tables_dir, "05_bulk_projection_cox_summary.csv"))
fwrite(assoc_dt, file.path(tables_dir, "05_bulk_projection_clinic_association_summary.csv"))

heatmap_dt <- all_scores_dt[, .(
  mean_module_score = mean(module_score, na.rm = TRUE)
), by = .(dataset, regulon)]
heatmap_dt[, dataset := factor(dataset, levels = c("TCGA", "ACRG", "GSE15459", "GSE84437"))]
heatmap_dt[, regulon := factor(regulon, levels = key_regulons)]

p_heatmap <- ggplot(heatmap_dt, aes(x = dataset, y = regulon, fill = mean_module_score)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank()) +
  labs(x = NULL, y = NULL, fill = "Mean\nscore", title = "Bulk projection of key regulon target modules")

tcga_assoc <- unique(all_scores_dt[dataset == "TCGA", .(
  ID, regulon, module_score, lauren_simple, stage_simple
)])
p_tcga_lauren <- ggplot(
  tcga_assoc[!is.na(lauren_simple) & lauren_simple %in% c("Intestinal", "Diffuse", "Mixed")],
  aes(x = lauren_simple, y = module_score, fill = lauren_simple)
) +
  geom_boxplot(outlier.size = 0.25, width = 0.72) +
  facet_wrap(~regulon, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(x = NULL, y = "Module score", title = "TCGA Lauren subtype association")

km_specs <- data.table(
  dataset = c("TCGA", "ACRG", "ACRG"),
  regulon = c("MAX(+)", "NR1H4(+)", "HNF4G(+)")
)
km_specs <- km_specs[regulon %in% unique(all_scores_dt$regulon)]
km_plots <- lapply(seq_len(nrow(km_specs)), function(i) {
  make_km_plot(
    all_scores_dt[dataset == km_specs$dataset[i] & regulon == km_specs$regulon[i]],
    km_specs$dataset[i],
    km_specs$regulon[i]
  )
})
km_panel <- if (length(km_plots) > 0) {
  patchwork::wrap_plots(km_plots, ncol = 2)
} else {
  ggplot() + theme_void(base_size = 11) + labs(title = "No available key regulons for KM plots")
}

figure5_bulk <- p_heatmap / p_tcga_lauren / km_panel +
  plot_annotation(tag_levels = "A")

ggsave(
  file.path(figures_dir, "Figure5_bulk_projection_main.png"),
  figure5_bulk,
  width = 13,
  height = 16,
  dpi = 320
)
ggsave(
  file.path(figures_dir, "Figure5_bulk_projection_main.pdf"),
  figure5_bulk,
  width = 13,
  height = 16,
  dpi = 320
)
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(file.path(figures_dir, "Figure5_bulk_projection_main.svg"), width = 13, height = 16)
  print(figure5_bulk)
  dev.off()
} else {
  ggsave(
    file.path(figures_dir, "Figure5_bulk_projection_main.svg"),
    figure5_bulk,
    device = grDevices::svg,
    width = 13,
    height = 16,
    units = "in",
    dpi = 320,
    limitsize = FALSE
  )
  message("svglite is unavailable; used the base R SVG device for bulk projection.")
}

writeLines(
  c(
    "# Figure 5 bulk projection note",
    "",
    "## Key regulon target modules",
    paste0("- Regulons: ", paste(key_regulons, collapse = ", ")),
    "- Module score definition: mean gene-level z score of regulon target genes within each bulk cohort.",
    "",
    "## Bulk cohorts",
    "- TCGA-STAD",
    "- ACRG / GSE62254",
    "- GSE15459",
    "- GSE84437",
    "",
    "## Outputs",
    "- Sample-level bulk projection tables for each cohort",
    "- Combined bulk projection table",
    "- Cox summary by dataset and regulon",
    "- Main bulk projection figure"
  ),
  file.path(manuscript_dir, "Figure5_bulk_projection_note.md")
)

cat("Figure 5 bulk projection completed.\n")
cat("Main figure:", file.path(figures_dir, "Figure5_bulk_projection_main.pdf"), "\n")
cat("Combined table:", file.path(tables_dir, "05_bulk_projection_all_samples_combined.csv"), "\n")

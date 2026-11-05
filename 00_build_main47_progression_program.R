project_root <- "/Users/huahaijiao/Desktop/02_study/0206_artical/Hp_gastric_cancer"
final_root <- file.path(project_root, "gastric_cancer", "results", "最终投稿版")
main47_root <- file.path(project_root, "gastric_cancer", "results", "progression_score_main47")
tables_dir <- file.path(main47_root, "tables")
figures_dir <- file.path(main47_root, "figures")
manuscript_dir <- file.path(main47_root, "manuscript")
final_tables_dir <- file.path(final_root, "tables", "main")
final_figures_dir <- file.path(final_root, "figures", "main")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(final_figures_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(readxl)
  library(patchwork)
  library(metafor)
  library(svglite)
})

normalize_gene <- function(x) {
  toupper(trimws(as.character(x)))
}

clean_chr <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | toupper(x) %in% c("NA", "N/A", "NULL")] <- NA_character_
  x
}

scale_by_gene <- function(expr_mat) {
  expr_mat <- as.matrix(expr_mat)
  z <- t(scale(t(expr_mat)))
  z[!is.finite(z)] <- NA_real_
  z
}

read_expr_rdata <- function(path, object_name = NULL) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (is.null(object_name)) {
    object_name <- ls(env)[1]
  }
  as.matrix(get(object_name, envir = env))
}

safe_read_csv <- function(path) {
  for (enc in c("UTF-8", "UTF-8-BOM", "GB18030", "GBK", "latin1")) {
    out <- try(fread(path, encoding = enc), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(as.data.table(out))
    }
  }
  stop("Failed to read: ", path)
}

calc_progression_score <- function(expr_mat, up_genes, down_genes) {
  expr_mat <- as.matrix(expr_mat)
  rownames(expr_mat) <- normalize_gene(rownames(expr_mat))
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), , drop = FALSE]
  zmat <- scale_by_gene(expr_mat)
  up_use <- intersect(up_genes, rownames(zmat))
  down_use <- intersect(down_genes, rownames(zmat))
  if (length(up_use) == 0L || length(down_use) == 0L) {
    stop("Both up and down gene sets must have at least one gene in the matrix.")
  }
  data.table(
    ID = colnames(zmat),
    progression_score_47g =
      colMeans(zmat[up_use, , drop = FALSE], na.rm = TRUE) -
      colMeans(zmat[down_use, , drop = FALSE], na.rm = TRUE),
    n_up_47g = length(up_use),
    n_down_47g = length(down_use),
    n_total_47g = length(up_use) + length(down_use)
  )
}

load_acrg_clin <- function() {
  clin <- as.data.table(read_excel(
    file.path(project_root, "clin_GSE66229_mRNAexp原始数据.xlsx"),
    sheet = "Sheet1"
  ))
  setnames(clin, "GEO_ID", "ID")
  clin[, dataset := "ACRG"]
  clin[, lauren_simple := fifelse(
    grepl("intestinal", clean_chr(Lauren), ignore.case = TRUE), "Intestinal",
    fifelse(
      grepl("diffuse", clean_chr(Lauren), ignore.case = TRUE), "Diffuse",
      fifelse(grepl("mixed", clean_chr(Lauren), ignore.case = TRUE), "Mixed", NA_character_)
    )
  )]
  clin[, acrg_subtype := clean_chr(`ACRG.sub`)]
  clin[, recurrence_group := fifelse(
    suppressWarnings(as.numeric(Recur)) == 1, "Recurrence",
    fifelse(suppressWarnings(as.numeric(Recur)) == 0, "Non-recurrence", NA_character_)
  )]
  clin[, recurrence_time_months := suppressWarnings(as.numeric(`DFS.m`))]
  clin[, recurrence_event := suppressWarnings(as.numeric(Recur))]
  clin[, OS_time := suppressWarnings(as.numeric(`OS.m`))]
  clin[, OS_event := suppressWarnings(as.numeric(Death))]
  clin[]
}

load_gse26253_clin <- function() {
  clin <- safe_read_csv(file.path(project_root, "clin_GSE26253_mRNAexp.csv"))
  setnames(clin, names(clin)[1], "ID")
  clin[, dataset := "GSE26253"]
  clin[, lauren_simple := fifelse(
    grepl("intestinal", clean_chr(Lauren), ignore.case = TRUE), "Intestinal",
    fifelse(
      grepl("diffuse", clean_chr(Lauren), ignore.case = TRUE), "Diffuse",
      fifelse(grepl("mixed", clean_chr(Lauren), ignore.case = TRUE), "Mixed", NA_character_)
    )
  )]
  clin[, recurrence_group := fifelse(
    suppressWarnings(as.numeric(`Recurrence (1=yes, 0=no)`)) == 1, "Recurrence",
    fifelse(
      suppressWarnings(as.numeric(`Recurrence (1=yes, 0=no)`)) == 0,
      "Non-recurrence",
      NA_character_
    )
  )]
  clin[, recurrence_time_months := suppressWarnings(as.numeric(`RFS.m`))]
  clin[, recurrence_event := suppressWarnings(as.numeric(`Recurrence (1=yes, 0=no)`))]
  clin[, OS_time := suppressWarnings(as.numeric(`OS.m`))]
  clin[, OS_event := suppressWarnings(as.numeric(`Death (1=yes, 0=no)`))]
  clin[]
}

load_gse26899_clin <- function() {
  clin <- safe_read_csv(file.path(project_root, "clin_GSE26899_mRNAexp.csv"))
  setnames(clin, names(clin)[1], "ID")
  clin <- clin[tolower(clean_chr(group)) == "gastric cancer"]
  clin[, dataset := "GSE26899"]
  clin[, lauren_simple := fifelse(
    grepl("intestinal", clean_chr(Lauren), ignore.case = TRUE), "Intestinal",
    fifelse(
      grepl("diffuse", clean_chr(Lauren), ignore.case = TRUE), "Diffuse",
      fifelse(grepl("mixed", clean_chr(Lauren), ignore.case = TRUE), "Mixed", NA_character_)
    )
  )]
  clin[, recurrence_group := fifelse(
    suppressWarnings(as.numeric(`Recurrence (1=yes, 0=no)`)) == 1, "Recurrence",
    fifelse(
      suppressWarnings(as.numeric(`Recurrence (1=yes, 0=no)`)) == 0,
      "Non-recurrence",
      NA_character_
    )
  )]
  clin[, recurrence_time_months := suppressWarnings(as.numeric(`RFS.m`))]
  clin[, recurrence_event := suppressWarnings(as.numeric(`Recurrence (1=yes, 0=no)`))]
  clin[, OS_time := suppressWarnings(as.numeric(`OS.m`))]
  clin[, OS_event := suppressWarnings(as.numeric(`Death (1=yes, 0=no)`))]
  clin[]
}

load_tcga_clin <- function() {
  clin <- safe_read_csv(file.path(project_root, "clin_TCGA_add_lauren.csv"))
  setnames(clin, names(clin)[1], "ID")
  clin[, dataset := "TCGA-STAD"]
  clin[, lauren_simple := fifelse(
    grepl("INTEST", toupper(clean_chr(luanren_type))), "Intestinal",
    fifelse(
      grepl("DIFF", toupper(clean_chr(luanren_type))), "Diffuse",
      fifelse(grepl("MIX", toupper(clean_chr(luanren_type))), "Mixed", NA_character_)
    )
  )]
  clin[, recurrence_group := fifelse(
    suppressWarnings(as.numeric(DFS_STATUS)) == 1, "Recurrence",
    fifelse(suppressWarnings(as.numeric(DFS_STATUS)) == 0, "Non-recurrence", NA_character_)
  )]
  clin[, recurrence_time_months := suppressWarnings(as.numeric(DFS_MONTHS))]
  clin[, recurrence_event := suppressWarnings(as.numeric(DFS_STATUS))]
  clin[, OS_time := suppressWarnings(as.numeric(OS_MONTHS))]
  clin[is.na(OS_time) & !is.na(DAYS_TO_LAST_FOLLOWUP),
       OS_time := suppressWarnings(as.numeric(DAYS_TO_LAST_FOLLOWUP)) / 30.4375]
  clin[, OS_event := suppressWarnings(as.numeric(OS_STATUS))]
  clin[]
}

load_gse15459_clin <- function() {
  clin <- safe_read_csv(file.path(project_root, "clin_GSE15459_mRNAexp.csv"))
  setnames(clin, names(clin)[1], "ID")
  clin[, dataset := "GSE15459"]
  clin[, lauren_simple := fifelse(
    grepl("intestinal", clean_chr(Laurenclassification), ignore.case = TRUE), "Intestinal",
    fifelse(
      grepl("diffuse", clean_chr(Laurenclassification), ignore.case = TRUE), "Diffuse",
      fifelse(grepl("mixed", clean_chr(Laurenclassification), ignore.case = TRUE), "Mixed", NA_character_)
    )
  )]
  clin[, recurrence_group := NA_character_]
  clin[, recurrence_time_months := NA_real_]
  clin[, recurrence_event := NA_real_]
  clin[, OS_time := suppressWarnings(as.numeric(OS_Months))]
  clin[, OS_event := suppressWarnings(as.numeric(OS_status))]
  clin[]
}

load_gse84437_clin <- function() {
  clin <- safe_read_csv(file.path(project_root, "clin_GSE84437_mRNAexp.csv"))
  setnames(clin, names(clin)[1], "ID")
  clin[, dataset := "GSE84437"]
  clin[, lauren_simple := NA_character_]
  clin[, recurrence_group := NA_character_]
  clin[, recurrence_time_months := NA_real_]
  clin[, recurrence_event := NA_real_]
  clin[, OS_time := suppressWarnings(as.numeric(OS_month))]
  clin[, OS_event := suppressWarnings(as.numeric(OS_status))]
  clin[]
}

hedges_g <- function(group_a, group_b) {
  group_a <- group_a[is.finite(group_a)]
  group_b <- group_b[is.finite(group_b)]
  n_a <- length(group_a)
  n_b <- length(group_b)
  if (n_a < 3L || n_b < 3L) {
    return(list(est = NA_real_, se = NA_real_, var = NA_real_))
  }
  s_a <- sd(group_a)
  s_b <- sd(group_b)
  pooled_sd <- sqrt(((n_a - 1) * s_a^2 + (n_b - 1) * s_b^2) / (n_a + n_b - 2))
  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    return(list(est = NA_real_, se = NA_real_, var = NA_real_))
  }
  d <- (mean(group_a) - mean(group_b)) / pooled_sd
  correction <- 1 - 3 / (4 * (n_a + n_b) - 9)
  g <- correction * d
  var_g <- (n_a + n_b) / (n_a * n_b) + g^2 / (2 * (n_a + n_b - 2))
  list(est = g, se = sqrt(var_g), var = var_g)
}

make_effect_row <- function(dt, dataset_name, analysis_name, group_col, group_a, group_b) {
  sub <- dt[
    get(group_col) %in% c(group_a, group_b) &
      is.finite(progression_score_47g)
  ]
  x_a <- sub[get(group_col) == group_a, progression_score_47g]
  x_b <- sub[get(group_col) == group_b, progression_score_47g]
  eff <- hedges_g(x_a, x_b)
  if (!is.finite(eff$est)) {
    return(NULL)
  }
  p_value <- tryCatch(
    wilcox.test(x_a, x_b, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  data.table(
    dataset = dataset_name,
    analysis = analysis_name,
    group_a = group_a,
    group_b = group_b,
    n_group_a = length(x_a),
    n_group_b = length(x_b),
    mean_group_a = mean(x_a),
    mean_group_b = mean(x_b),
    median_group_a = median(x_a),
    median_group_b = median(x_b),
    effect_size = eff$est,
    effect_direction = ifelse(
      analysis_name == "Lauren intestinal vs diffuse",
      "positive = intestinal higher than diffuse",
      "positive = non-recurrence higher than recurrence"
    ),
    se = eff$se,
    var = eff$var,
    ci_low = eff$est - 1.96 * eff$se,
    ci_high = eff$est + 1.96 * eff$se,
    p_value = p_value
  )
}

fit_meta <- function(dt) {
  dt <- dt[is.finite(effect_size) & is.finite(se) & se > 0]
  fit <- rma.uni(yi = effect_size, sei = se, data = dt, method = "REML")
  data.table(
    pooled_effect = as.numeric(fit$b[1, 1]),
    ci_low = fit$ci.lb,
    ci_high = fit$ci.ub,
    p_value = fit$pval,
    tau2 = fit$tau2,
    i2 = fit$I2,
    qe_p = fit$QEp,
    included_datasets = paste(dt$dataset, collapse = ";"),
    n_datasets = nrow(dt)
  )
}

score_def <- fread(file.path(
  project_root,
  "gastric_cancer", "results", "Figure1_integrated_candidate_selection_figure",
  "r_repro_main47", "tables", "03_final_cross_platform_candidates_main47.csv"
))
score_def[, gene := normalize_gene(gene)]
up_genes <- score_def[direction == "up_gain", gene]
down_genes <- score_def[direction == "down_loss", gene]

expr_list <- list(
  ACRG = read_expr_rdata(file.path(project_root, "GSE62254_mRNAexp.RData"), "GSE62254_mRNAexp"),
  GSE26253 = read_expr_rdata(file.path(project_root, "GSE26253_mRNAexp.RData"), "GSE26253_mRNAexp"),
  GSE26899 = read_expr_rdata(file.path(project_root, "GSE26899_mRNAexp.RData"), "GSE26899_mRNAexp"),
  `TCGA-STAD` = read_expr_rdata(file.path(project_root, "TCGA_GC_mRNAexp.RData"), "TCGA_GC_mRNAexp"),
  GSE15459 = read_expr_rdata(file.path(project_root, "GSE15459_mRNAexp.RData"), "GSE15459_mRNAexp"),
  GSE84437 = read_expr_rdata(file.path(project_root, "GSE84437_mRNAexp.RData"), "GSE84437_mRNAexp")
)

clin_list <- list(
  ACRG = load_acrg_clin(),
  GSE26253 = load_gse26253_clin(),
  GSE26899 = load_gse26899_clin(),
  `TCGA-STAD` = load_tcga_clin(),
  GSE15459 = load_gse15459_clin(),
  GSE84437 = load_gse84437_clin()
)

score_tables <- list()
dataset_summary <- list()
for (dataset_name in names(expr_list)) {
  score_dt <- calc_progression_score(expr_list[[dataset_name]], up_genes, down_genes)
  score_dt[, dataset := dataset_name]
  clin_dt <- copy(clin_list[[dataset_name]])
  merged <- merge(score_dt, clin_dt, by = c("ID", "dataset"), all.x = TRUE)
  merged[, score_z := as.numeric(scale(progression_score_47g))]
  merged[, score_group := fifelse(
    progression_score_47g >= median(progression_score_47g, na.rm = TRUE),
    "High", "Low"
  )]
  score_tables[[dataset_name]] <- merged
  fwrite(
    merged,
    file.path(tables_dir, paste0("progression_score_47g_", gsub("-", "_", dataset_name), "_sample_table.csv"))
  )
  dataset_summary[[dataset_name]] <- data.table(
    dataset = dataset_name,
    n_samples_expr = ncol(expr_list[[dataset_name]]),
    n_samples_clin = nrow(clin_dt),
    n_samples_merged = nrow(merged),
    n_lauren_usable = sum(merged$lauren_simple %in% c("Intestinal", "Diffuse"), na.rm = TRUE),
    n_recurrence_usable = sum(merged$recurrence_group %in% c("Non-recurrence", "Recurrence"), na.rm = TRUE),
    n_survival_usable = sum(is.finite(merged$OS_time) & !is.na(merged$OS_event)),
    n_up_present = unique(merged$n_up_47g)[1],
    n_down_present = unique(merged$n_down_47g)[1],
    n_total_present = unique(merged$n_total_47g)[1]
  )
}

all_scores <- rbindlist(score_tables, fill = TRUE)
summary_dt <- rbindlist(dataset_summary, fill = TRUE)
fwrite(summary_dt, file.path(tables_dir, "progression_score_47g_dataset_summary.csv"))
fwrite(all_scores, file.path(tables_dir, "progression_score_47g_all_samples_combined.csv"))
fwrite(score_def, file.path(tables_dir, "progression_score_47g_definition.csv"))

lauren_datasets <- c("ACRG", "GSE26253", "GSE26899", "TCGA-STAD")
lauren_dt <- rbindlist(lapply(lauren_datasets, function(ds) {
  make_effect_row(
    score_tables[[ds]],
    ds,
    "Lauren intestinal vs diffuse",
    "lauren_simple",
    "Intestinal",
    "Diffuse"
  )
}), fill = TRUE)
lauren_dt[, p_value_fdr := p.adjust(p_value, method = "BH")]
lauren_meta <- cbind(
  data.table(
    analysis = "Lauren intestinal vs diffuse",
    meta_set = "all_Lauren_cohorts"
  ),
  fit_meta(lauren_dt)
)

recurrence_datasets <- c("ACRG", "GSE26253", "GSE26899", "TCGA-STAD")
recurrence_dt <- rbindlist(lapply(recurrence_datasets, function(ds) {
  make_effect_row(
    score_tables[[ds]],
    ds,
    ifelse(ds == "TCGA-STAD",
           "DFS proxy non-recurrence vs recurrence",
           "Recurrence non-recurrence vs recurrence"),
    "recurrence_group",
    "Non-recurrence",
    "Recurrence"
  )
}), fill = TRUE)
recurrence_dt[, p_value_fdr := p.adjust(p_value, method = "BH")]

meta_primary <- cbind(
  data.table(
    meta_set = "primary_only",
    analysis = "Recurrence non-recurrence vs recurrence"
  ),
  fit_meta(recurrence_dt[dataset %in% c("ACRG", "GSE26253")])
)
meta_primary_proxy <- cbind(
  data.table(
    meta_set = "primary_plus_proxy",
    analysis = "Recurrence non-recurrence vs recurrence"
  ),
  fit_meta(recurrence_dt[dataset %in% c("ACRG", "GSE26253", "TCGA-STAD")])
)

fwrite(lauren_dt, file.path(tables_dir, "Table_3_lauren_effect_sizes_47g.csv"))
fwrite(recurrence_dt, file.path(tables_dir, "Table_4_recurrence_effect_sizes_47g.csv"))
fwrite(rbind(meta_primary, meta_primary_proxy, fill = TRUE),
       file.path(tables_dir, "Table_5_recurrence_meta_47g.csv"))
fwrite(rbind(lauren_meta, meta_primary, meta_primary_proxy, fill = TRUE),
       file.path(tables_dir, "Table_6_cross_cohort_meta_summary_47g.csv"))

file.copy(
  file.path(tables_dir, "Table_3_lauren_effect_sizes_47g.csv"),
  file.path(final_tables_dir, "Table_3_lauren_effect_sizes.csv"),
  overwrite = TRUE
)
file.copy(
  file.path(tables_dir, "Table_4_recurrence_effect_sizes_47g.csv"),
  file.path(final_tables_dir, "Table_4_recurrence_effect_sizes.csv"),
  overwrite = TRUE
)
file.copy(
  file.path(tables_dir, "Table_5_recurrence_meta_47g.csv"),
  file.path(final_tables_dir, "Table_5_recurrence_meta_primary_plus_proxy.csv"),
  overwrite = TRUE
)
file.copy(
  file.path(tables_dir, "Table_6_cross_cohort_meta_summary_47g.csv"),
  file.path(final_tables_dir, "Table_6_cross_cohort_meta_summary.csv"),
  overwrite = TRUE
)

acrg <- score_tables[["ACRG"]]
acrg[, subtype_clean := factor(
  acrg_subtype,
  levels = c("MSI", "TP53positive", "TP53neg", "EMT")
)]
subtype_p <- kruskal.test(
  progression_score_47g ~ subtype_clean,
  data = acrg[!is.na(subtype_clean)]
)$p.value

program_sets <- list(
  `47-gene up-gain` = up_genes,
  `47-gene down-loss` = down_genes,
  `Regenerative/Malignant` = c(
    "OLFM4", "REG4", "REG1A", "TSPAN8", "CLDN4",
    "MDK", "EPCAM", "CLDN7", "MUC13", "PRAP1"
  ),
  `Intestinalized IM` = c(
    "PHGR1", "FABP1", "FABP2", "RBP2", "ANPEP",
    "LGALS4", "CDH17", "VIL1", "TFF3"
  ),
  `Gastric Homeostasis` = c(
    "GKN1", "GKN2", "TFF1", "TFF2", "MUC5AC",
    "MUC6", "PGC", "CLDN18", "GAST"
  ),
  `EMT-like` = c("VIM", "COL1A1", "COL1A2", "FN1", "SPARC", "ZEB1")
)

acrg_expr <- expr_list[["ACRG"]]
rownames(acrg_expr) <- normalize_gene(rownames(acrg_expr))
acrg_expr <- acrg_expr[!duplicated(rownames(acrg_expr)), , drop = FALSE]
acrg_z <- scale_by_gene(acrg_expr)
common_ids <- intersect(colnames(acrg_z), acrg$ID)
acrg_cor <- acrg[match(common_ids, ID)]
acrg_cor[, progression_score_47g := score_tables[["ACRG"]][
  match(common_ids, ID), progression_score_47g
]]
for (program_name in names(program_sets)) {
  genes <- intersect(normalize_gene(program_sets[[program_name]]), rownames(acrg_z))
  acrg_cor[, (program_name) := colMeans(acrg_z[genes, common_ids, drop = FALSE], na.rm = TRUE)]
}
cor_dt <- rbindlist(lapply(names(program_sets), function(program_name) {
  data.table(
    feature = program_name,
    spearman_rho = suppressWarnings(cor(
      acrg_cor$progression_score_47g,
      acrg_cor[[program_name]],
      method = "spearman",
      use = "pairwise.complete.obs"
    ))
  )
}), fill = TRUE)
fwrite(cor_dt, file.path(tables_dir, "figure3_acrg_47g_program_correlations.csv"))

p_subtype <- ggplot(
  acrg[!is.na(subtype_clean)],
  aes(subtype_clean, progression_score_47g, fill = subtype_clean)
) +
  geom_violin(trim = FALSE, alpha = 0.9, color = "white", linewidth = 0.3) +
  geom_boxplot(width = 0.18, outlier.size = 0.2, fill = "white") +
  scale_fill_manual(values = c(
    "MSI" = "#4C8C61",
    "TP53positive" = "#D98C2B",
    "TP53neg" = "#7A7A7A",
    "EMT" = "#BB4B39"
  )) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    axis.text.x = element_text(angle = 25, hjust = 1)
  ) +
  labs(
    x = NULL,
    y = "47-gene progression score",
    title = "ACRG molecular subtype",
    subtitle = sprintf("Kruskal-Wallis P = %.3g", subtype_p)
  )

make_forest_panel <- function(dt, pooled_row, title_text, xlab_text, proxy_label = FALSE) {
  plot_dt <- copy(dt)
  plot_dt[, display_dataset := dataset]
  if (proxy_label) {
    plot_dt[dataset == "TCGA-STAD", display_dataset := "TCGA-STAD (DFS proxy)"]
  }
  pooled_dt <- data.table(
    dataset = "Pooled",
    display_dataset = "Pooled",
    n_group_a = sum(plot_dt$n_group_a),
    n_group_b = sum(plot_dt$n_group_b),
    effect_size = pooled_row$pooled_effect,
    ci_low = pooled_row$ci_low,
    ci_high = pooled_row$ci_high
  )
  plot_dt <- rbind(plot_dt, pooled_dt, fill = TRUE)
  plot_dt[, label := ifelse(
    dataset == "Pooled",
    sprintf("Pooled (n=%d vs %d)", n_group_a, n_group_b),
    sprintf("%s (n=%d vs %d)", display_dataset, n_group_a, n_group_b)
  )]
  plot_dt[, label := factor(label, levels = rev(label))]
  plot_dt[, pooled := dataset == "Pooled"]
  ggplot(plot_dt, aes(x = effect_size, y = label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "#7A7A7A", linewidth = 0.4) +
    geom_segment(
      aes(x = ci_low, xend = ci_high, yend = label, color = pooled),
      linewidth = 0.9
    ) +
    geom_point(
      aes(fill = pooled, shape = pooled),
      size = 3.2,
      color = "black",
      stroke = 0.3
    ) +
    scale_color_manual(values = c("FALSE" = "#4F6FAE", "TRUE" = "#BB4B39"), guide = "none") +
    scale_fill_manual(values = c("FALSE" = "#4F6FAE", "TRUE" = "#BB4B39"), guide = "none") +
    scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 23), guide = "none") +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), axis.title.y = element_blank()) +
    labs(title = title_text, x = xlab_text)
}

p_lauren <- make_forest_panel(
  lauren_dt,
  lauren_meta,
  "Cross-cohort Lauren consistency",
  "Hedges' g (positive = intestinal higher)"
)

lauren_meta_for_plot <- lauren_meta
p_lauren <- p_lauren +
  labs(subtitle = sprintf(
    "Pooled g = %.2f (95%% CI %.2f to %.2f), P = %.3g, I2 = %.1f%%",
    lauren_meta_for_plot$pooled_effect,
    lauren_meta_for_plot$ci_low,
    lauren_meta_for_plot$ci_high,
    lauren_meta_for_plot$p_value,
    lauren_meta_for_plot$i2
  ))

p_corr <- ggplot(
  cor_dt,
  aes(
    1,
    factor(feature, levels = c(
      "47-gene up-gain",
      "47-gene down-loss",
      "Regenerative/Malignant",
      "Intestinalized IM",
      "Gastric Homeostasis",
      "EMT-like"
    )),
    fill = spearman_rho
  )
) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 3.6) +
  scale_fill_gradient2(low = "#3B6AA0", mid = "white", high = "#BB4B39", midpoint = 0) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank()
  ) +
  labs(title = "Correlation with biological programs", fill = "Spearman rho")

recurrence_plot_dt <- recurrence_dt[
  dataset %in% c("ACRG", "GSE26253", "TCGA-STAD"),
  .(
    dataset,
    n_group_a,
    n_group_b,
    effect_size,
    ci_low,
    ci_high
  )
]
p_recurrence <- make_forest_panel(
  recurrence_plot_dt,
  meta_primary_proxy,
  "Recurrence endpoint support",
  "Hedges' g (positive = non-recurrence higher)",
  proxy_label = TRUE
) +
  labs(subtitle = sprintf(
    "Primary + TCGA DFS proxy: g = %.2f (95%% CI %.2f to %.2f), P = %.3g, I2 = %.1f%%",
    meta_primary_proxy$pooled_effect,
    meta_primary_proxy$ci_low,
    meta_primary_proxy$ci_high,
    meta_primary_proxy$p_value,
    meta_primary_proxy$i2
  ))

figure3_main <- (p_subtype | p_lauren) / (p_corr | p_recurrence) +
  plot_layout(widths = c(1, 1.35), heights = c(1, 1)) +
  plot_annotation(
    title = "Figure 3. 47-gene progression program in Lauren-oriented bulk cohorts",
    tag_levels = "A"
  )

save_plot <- function(plot_obj, path_stem, width, height) {
  ggsave(paste0(path_stem, ".pdf"), plot_obj, width = width, height = height)
  ggsave(paste0(path_stem, ".png"), plot_obj, width = width, height = height, dpi = 320)
  svglite::svglite(paste0(path_stem, ".svg"), width = width, height = height)
  print(plot_obj)
  dev.off()
}

save_plot(
  figure3_main,
  file.path(figures_dir, "Figure3_main_47g_progression_program"),
  15.5,
  11.5
)
file.copy(
  file.path(figures_dir, "Figure3_main_47g_progression_program.pdf"),
  file.path(final_figures_dir, "Figure3_main.pdf"),
  overwrite = TRUE
)
file.copy(
  file.path(figures_dir, "Figure3_main_47g_progression_program.png"),
  file.path(final_figures_dir, "Figure3_main.png"),
  overwrite = TRUE
)
file.copy(
  file.path(figures_dir, "Figure3_main_47g_progression_program.svg"),
  file.path(final_figures_dir, "Figure3_main.svg"),
  overwrite = TRUE
)

writeLines(
  c(
    "# 47-gene progression program analysis",
    "",
    "The 47-gene score is the primary progression-program readout.",
    "It is calculated as the mean row-wise z-score of the 38 up-gain genes minus the mean row-wise z-score of the 9 down-loss genes within each cohort.",
    "Positive Lauren effect sizes are defined as intestinal minus diffuse.",
    "The 12-gene panel is retained only as a compact auxiliary readout and is not used for the primary Figure 3 score."
  ),
  file.path(manuscript_dir, "main47_progression_program_note.md")
)

cat("Main 47-gene progression program analysis completed.\n")
cat("Lauren rows:", nrow(lauren_dt), "\n")
cat("Recurrence rows:", nrow(recurrence_dt), "\n")
cat("47-gene score tables:", tables_dir, "\n")

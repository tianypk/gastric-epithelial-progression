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
tables_dir <- file.path(result_root, "tables")
rds_dir <- file.path(result_root, "rds")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

safe_wilcox <- function(score, group) {
  ok <- is.finite(score) & !is.na(group) & group != ""
  score <- score[ok]
  group <- group[ok]
  if (length(score) < 8 || length(unique(group)) != 2) {
    return(NA_real_)
  }
  tryCatch(
    wilcox.test(score ~ group, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
}

safe_logistic <- function(dt, score_col, group_col, positive_label) {
  sub <- dt[is.finite(get(score_col)) & !is.na(get(group_col))]
  sub <- sub[get(group_col) %in% c(positive_label, setdiff(unique(get(group_col)), positive_label))]
  if (nrow(sub) < 10 || length(unique(sub[[group_col]])) != 2) {
    return(data.table(glm_or = NA_real_, glm_lower95 = NA_real_, glm_upper95 = NA_real_, glm_pvalue = NA_real_))
  }
  sub[, y := as.integer(get(group_col) == positive_label)]
  sub[, score_z := as.numeric(scale(get(score_col)))]
  fit <- tryCatch(glm(y ~ score_z, data = sub, family = binomial()), error = function(e) NULL)
  if (is.null(fit)) {
    return(data.table(glm_or = NA_real_, glm_lower95 = NA_real_, glm_upper95 = NA_real_, glm_pvalue = NA_real_))
  }
  sm <- summary(fit)$coefficients
  ci <- suppressWarnings(tryCatch(confint.default(fit), error = function(e) NULL))
  if (is.null(ci)) {
    return(data.table(glm_or = NA_real_, glm_lower95 = NA_real_, glm_upper95 = NA_real_, glm_pvalue = NA_real_))
  }
  data.table(
    glm_or = exp(unname(sm["score_z", "Estimate"])),
    glm_lower95 = exp(unname(ci["score_z", 1])),
    glm_upper95 = exp(unname(ci["score_z", 2])),
    glm_pvalue = unname(sm["score_z", "Pr(>|z|)"])
  )
}

required_files <- c(
  file.path(rds_dir, "02_prjeb25780_scored_dataset.rds")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required files:\n", paste(missing_files, collapse = "\n"))
}

scored_bundle <- readRDS(file.path(rds_dir, "02_prjeb25780_scored_dataset.rds"))
score_dt <- as.data.table(scored_bundle$scored_data)

score_cols <- c(
  "progression_score_47g",
  "MAX_module_score",
  "NR1H4_module_score",
  "HNF4G_module_score"
)
group_configs <- list(
  response = list(column = "response_group", positive = "responder", group1 = "responder", group2 = "non_responder"),
  pd_l1_status = list(column = "pd_l1_status", positive = "positive", group1 = "positive", group2 = "negative"),
  msi_status_pentaplex = list(column = "msi_status_pentaplex", positive = "MSI", group1 = "MSI", group2 = "MSS"),
  ebv_status = list(column = "ebv_status", positive = "EBV_positive", group1 = "EBV_positive", group2 = "EBV_not_tcga_ebv_subtype")
)

summary_list <- list()

for (analysis_name in names(group_configs)) {
  cfg <- group_configs[[analysis_name]]
  for (score_col in score_cols) {
    sub <- copy(score_dt)
    sub <- sub[!is.na(get(cfg$column)) & !is.na(get(score_col))]
    sub <- sub[get(cfg$column) %in% c(cfg$group1, cfg$group2)]
    if (nrow(sub) == 0) {
      next
    }
    group1_vals <- sub[get(cfg$column) == cfg$group1][[score_col]]
    group2_vals <- sub[get(cfg$column) == cfg$group2][[score_col]]
    stat_dt <- data.table(
      analysis_group = analysis_name,
      score_name = score_col,
      comparison = paste(cfg$group1, "vs", cfg$group2),
      n_total = nrow(sub),
      group1_label = cfg$group1,
      group2_label = cfg$group2,
      group1_n = length(group1_vals),
      group2_n = length(group2_vals),
      group1_median = median(group1_vals, na.rm = TRUE),
      group2_median = median(group2_vals, na.rm = TRUE),
      delta_median = median(group1_vals, na.rm = TRUE) - median(group2_vals, na.rm = TRUE),
      wilcox_pvalue = safe_wilcox(sub[[score_col]], sub[[cfg$column]])
    )
    glm_dt <- if (analysis_name == "response") {
      safe_logistic(sub, score_col, cfg$column, cfg$positive)
    } else {
      data.table(glm_or = NA_real_, glm_lower95 = NA_real_, glm_upper95 = NA_real_, glm_pvalue = NA_real_)
    }
    stat_dt <- cbind(stat_dt, glm_dt)
    stat_dt[, direction_label := fifelse(delta_median > 0, cfg$group1, cfg$group2)]
    summary_list[[length(summary_list) + 1L]] <- stat_dt
  }
}

association_dt <- rbindlist(summary_list, fill = TRUE)
response_dt <- association_dt[analysis_group == "response"]
biomarker_dt <- association_dt[analysis_group != "response"]

plotting_long_dt <- melt(
  score_dt[, c(
    "sample_id", "patient_id", "response_group",
    "pd_l1_status", "msi_status_pentaplex", "ebv_status",
    score_cols
  ), with = FALSE],
  id.vars = c("sample_id", "patient_id", "response_group", "pd_l1_status", "msi_status_pentaplex", "ebv_status"),
  variable.name = "score_name",
  value.name = "score_value"
)

endpoint_readiness_dt <- data.table(
  endpoint_name = c("response", "pd_l1_status", "msi_status_pentaplex", "ebv_status", "pfs", "os"),
  n_non_missing = c(
    score_dt[, sum(!is.na(response_group))],
    score_dt[, sum(!is.na(pd_l1_status) & pd_l1_status != "N.A.")],
    score_dt[, sum(!is.na(msi_status_pentaplex))],
    score_dt[, sum(!is.na(ebv_status))],
    score_dt[, sum(!is.na(pfs_time))],
    score_dt[, sum(!is.na(os_time))]
  ),
  n_missing = nrow(score_dt) - c(
    score_dt[, sum(!is.na(response_group))],
    score_dt[, sum(!is.na(pd_l1_status) & pd_l1_status != "N.A.")],
    score_dt[, sum(!is.na(msi_status_pentaplex))],
    score_dt[, sum(!is.na(ebv_status))],
    score_dt[, sum(!is.na(pfs_time))],
    score_dt[, sum(!is.na(os_time))]
  ),
  recommended_use = c(
    "primary",
    "biomarker_context",
    "biomarker_context",
    "biomarker_context",
    "placeholder_only",
    "placeholder_only"
  ),
  panel_d_strategy = c(
    "OR_forest",
    "context_only",
    "context_only",
    "context_only",
    "not_ready_for_main_figure",
    "not_ready_for_main_figure"
  )
)

fwrite(response_dt, file.path(tables_dir, "03_prjeb25780_response_association_summary.csv"))
fwrite(biomarker_dt, file.path(tables_dir, "03_prjeb25780_biomarker_context_summary.csv"))
fwrite(plotting_long_dt, file.path(tables_dir, "03_prjeb25780_plotting_long_scores.csv"))
fwrite(endpoint_readiness_dt, file.path(tables_dir, "04_prjeb25780_paneld_endpoint_readiness.csv"))

cat("Figure 6 step 03 completed.\n")
cat("Response association rows:", nrow(response_dt), "\n")

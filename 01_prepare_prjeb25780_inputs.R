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
input_root <- file.path(project_root, "immunotherapy_datasets", "PRJEB25780")
tables_dir <- file.path(result_root, "tables")
rds_dir <- file.path(result_root, "rds")
logs_dir <- file.path(result_root, "logs")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)

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

input_files <- list(
  metadata_master = file.path(input_root, "PRJEB25780_paper_ready_metadata_master.csv"),
  metadata_tumor45 = file.path(input_root, "PRJEB25780_paper_ready_metadata_tumor45.csv"),
  clinical = file.path(input_root, "STAD-PRJEB25780_anti-PD-1.Response.clinical.tsv"),
  expression = file.path(input_root, "STAD-PRJEB25780_anti-PD-1.Response.expression.tsv")
)

missing_files <- names(input_files)[!file.exists(unlist(input_files))]
if (length(missing_files) > 0) {
  stop("Missing input files: ", paste(missing_files, collapse = ", "))
}

meta_master <- fread(input_files$metadata_master, fill = TRUE)
meta_tumor45 <- fread(input_files$metadata_tumor45, fill = TRUE)
clin_dt <- read_tsv_with_optional_index(input_files$clinical, "sample_id")
expr_header <- read_tsv_with_optional_index(input_files$expression, "GENE_SYMBOL")[0]
expr_sample_ids <- setdiff(colnames(expr_header), "GENE_SYMBOL")

manifest_dt <- copy(meta_tumor45)
manifest_dt[, sample_id := run_accession]
manifest_dt <- merge(
  manifest_dt,
  clin_dt[, .(
    sample_id,
    patient_name,
    tiger_response = response,
    tiger_response_NR = response_NR
  )],
  by = "sample_id",
  all.x = TRUE,
  suffixes = c("", "_clinical")
)

manifest_dt[, in_expression_matrix := sample_id %in% expr_sample_ids]
manifest_dt[, pfs_available := !is.na(pfs_time)]
manifest_dt[, os_available := !is.na(os_time)]
manifest_dt[, pd_l1_available := !is.na(pd_l1_status) & pd_l1_status != "N.A."]
manifest_dt[, msi_available := !is.na(msi_status_pentaplex)]
manifest_dt[, ebv_available := !is.na(ebv_status)]
manifest_dt[, response_available := !is.na(response_group)]

setcolorder(manifest_dt, c(
  "sample_id", "patient_id", "sample_alias", "study_patient_label",
  "sample_accession", "experiment_accession", "run_accession",
  "response", "response_group", "tiger_response", "tiger_response_NR",
  "pd_l1_status", "pd_l1_cps", "msi_status_pentaplex", "msi_wes",
  "ebv_status", "tcga_subtype", "immune_signature", "pathology",
  "pfs_time", "pfs_event", "os_time", "os_event",
  "response_available", "pd_l1_available", "msi_available", "ebv_available",
  "pfs_available", "os_available", "in_expression_matrix",
  "recommended_analysis_set", "annotation_confidence", "response_confidence",
  "msi_confidence", "ebv_confidence", "pd_l1_confidence",
  "sample_identity_source", "response_source", "notes"
))

cohort_overview_dt <- rbindlist(list(
  data.table(
    variable = c("tumor_rna_samples", "response", "pd_l1_status", "msi_status", "ebv_status", "pfs", "os"),
    n_known = c(
      manifest_dt[, sum(in_expression_matrix)],
      manifest_dt[, sum(response_available)],
      manifest_dt[, sum(pd_l1_available)],
      manifest_dt[, sum(msi_available)],
      manifest_dt[, sum(ebv_available)],
      manifest_dt[, sum(pfs_available)],
      manifest_dt[, sum(os_available)]
    ),
    n_total = nrow(manifest_dt),
    analysis_role = c(
      "primary_input",
      "primary_endpoint",
      "biomarker_context",
      "biomarker_context",
      "biomarker_context",
      "endpoint_placeholder",
      "endpoint_placeholder"
    )
  ),
  manifest_dt[, .(
    variable = "response_group_breakdown",
    level = response_group,
    n_known = .N,
    n_total = nrow(manifest_dt),
    analysis_role = "panel_b"
  ), by = response_group][, response_group := NULL],
  manifest_dt[, .(
    variable = "pd_l1_breakdown",
    level = pd_l1_status,
    n_known = .N,
    n_total = nrow(manifest_dt),
    analysis_role = "panel_c"
  ), by = pd_l1_status][, pd_l1_status := NULL]
), fill = TRUE)

panel_plan_dt <- data.table(
  panel = c("A", "B", "C", "D"),
  title_cn = c(
    "队列与数据可用性概览",
    "response 关联",
    "PD-L1 / MSI / EBV 生物标志物上下文",
    "response-oriented OR 汇总"
  ),
  main_input = c(
    "01_prjeb25780_cohort_overview.csv",
    "02_prjeb25780_projected_scores.csv",
    "03_prjeb25780_biomarker_context_summary.csv",
    "03_prjeb25780_response_association_summary.csv"
  ),
  current_status = c(
    "可直接生成",
    "待 score projection",
    "待 score projection",
    "待 score projection"
  ),
  notes = c(
    "用于交代 response/PD-L1/MSI/EBV/PFS/OS 覆盖率",
    "主回答 progression-regulon axis 与疗效关系",
    "用于说明与经典免疫标志物的关系",
    "当前公开资料缺乏 PFS/OS，因此先用 OR forest 代替 survival"
  )
)

bundle <- list(
  manifest = manifest_dt,
  cohort_overview = cohort_overview_dt,
  expr_sample_ids = expr_sample_ids,
  input_files = input_files
)

fwrite(panel_plan_dt, file.path(tables_dir, "00_figure6_panel_plan.csv"))
fwrite(manifest_dt, file.path(tables_dir, "01_prjeb25780_cohort_manifest.csv"))
fwrite(cohort_overview_dt, file.path(tables_dir, "01_prjeb25780_cohort_overview.csv"))
saveRDS(bundle, file.path(rds_dir, "01_prjeb25780_input_bundle.rds"))

cat("Figure 6 step 01 completed.\n")
cat("Manifest rows:", nrow(manifest_dt), "\n")

# 10_advanced_clinical_omics_analysis.R
# Extended analyses following a workflow like the provided figure:
# - Univariate Cox screening
# - LASSO-Cox signature and risk score
# - KM + time-dependent ROC + optional nomogram
# - Immune infiltration and immune-cell proportion
# - TMB and drug-sensitivity association analyses

source("pipeline/00_config.R")

prepare_survival_df <- function(df, time_col = "OS_time", event_col = "OS_event", feature_cols = NULL) {
  keep <- c(time_col, event_col, feature_cols)
  keep <- keep[keep %in% colnames(df)]
  out <- df[, keep, drop = FALSE]
  out <- out[stats::complete.cases(out), , drop = FALSE]
  out
}

run_univariate_cox <- function(df, time_col = "OS_time", event_col = "OS_event", feature_cols) {
  rows <- lapply(feature_cols, function(f) {
    if (!(f %in% colnames(df))) return(NULL)
    dat <- df[, c(time_col, event_col, f), drop = FALSE]
    dat <- dat[stats::complete.cases(dat), , drop = FALSE]
    if (nrow(dat) < 20) return(NULL)
    fit <- survival::coxph(stats::as.formula(paste0("survival::Surv(", time_col, ",", event_col, ") ~ `", f, "`")), data = dat)
    tbl <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
    tbl$feature <- f
    tbl
  })
  dplyr::bind_rows(rows)
}

run_lasso_cox <- function(df, time_col = "OS_time", event_col = "OS_event", feature_cols, nfolds = 5) {
  x <- as.matrix(df[, feature_cols, drop = FALSE])
  y <- survival::Surv(df[[time_col]], df[[event_col]])
  fit <- glmnet::cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = nfolds)
  coef_mat <- as.matrix(glmnet::coef.glmnet(fit, s = "lambda.1se"))
  nz <- which(coef_mat[, 1] != 0)
  selected <- rownames(coef_mat)[nz]

  risk <- as.numeric(stats::predict(fit, newx = x, s = "lambda.1se", type = "link"))
  list(fit = fit, selected = selected, coef = coef_mat[nz, , drop = FALSE], risk_score = risk)
}

evaluate_survival_signature <- function(df, risk_score, time_col = "OS_time", event_col = "OS_event", out_prefix = "lasso_cox") {
  out_df <- df
  out_df$risk_score <- risk_score
  cut <- stats::median(out_df$risk_score, na.rm = TRUE)
  out_df$risk_group <- factor(ifelse(out_df$risk_score >= cut, "High", "Low"))

  # KM curve
  if (requireNamespace("survminer", quietly = TRUE)) {
    sfit <- survival::survfit(survival::Surv(out_df[[time_col]], out_df[[event_col]]) ~ risk_group, data = out_df)
    p <- survminer::ggsurvplot(sfit, data = out_df, risk.table = TRUE, pval = TRUE, conf.int = FALSE, title = "Kaplan-Meier by LASSO-Cox risk group")
    ggplot2::ggsave(file.path(config$output$report_dir, paste0(out_prefix, "_km.png")), p$plot, width = 7, height = 5, dpi = 300)
  }

  # Time-dependent ROC
  if (requireNamespace("timeROC", quietly = TRUE)) {
    roc <- timeROC::timeROC(T = out_df[[time_col]],
                            delta = out_df[[event_col]],
                            marker = out_df$risk_score,
                            cause = 1,
                            times = c(12, 36, 60),
                            iid = TRUE)
    roc_tbl <- data.frame(time = roc$times, AUC = roc$AUC)
    data.table::fwrite(roc_tbl, file.path(config$output$report_dir, paste0(out_prefix, "_timeROC_auc.csv")))
  }

  out_df
}

build_nomogram_if_available <- function(df, time_col = "OS_time", event_col = "OS_event", covariates = c("risk_group", "age", "stage"), out_file = "nomogram.png") {
  if (!requireNamespace("rms", quietly = TRUE)) {
    warning("rms package not installed; nomogram skipped")
    return(NULL)
  }

  covariates <- covariates[covariates %in% colnames(df)]
  if (length(covariates) == 0) return(NULL)
  dat <- df[, c(time_col, event_col, covariates), drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  dd <- rms::datadist(dat)
  options(datadist = "dd")
  f <- stats::as.formula(paste0("survival::Surv(", time_col, ",", event_col, ") ~ ", paste0("`", covariates, "`", collapse = " + ")))
  fit <- rms::cph(f, data = dat, x = TRUE, y = TRUE, surv = TRUE)
  surv_fun <- rms::Survival(fit)
  nom <- rms::nomogram(fit,
                       fun = list(function(lp) surv_fun(12, lp),
                                  function(lp) surv_fun(36, lp),
                                  function(lp) surv_fun(60, lp)),
                       funlabel = c("1-year survival", "3-year survival", "5-year survival"),
                       lp = TRUE)

  png(file.path(config$output$report_dir, out_file), width = 1200, height = 900)
  plot(nom)
  dev.off()
  fit
}

run_immune_infiltration <- function(expr_mat, method = "xcell") {
  if (!requireNamespace("immunedeconv", quietly = TRUE)) {
    warning("immunedeconv not installed; immune infiltration skipped")
    return(NULL)
  }
  immunedeconv::deconvolute(expr_mat, method = method)
}

analyze_tmb_association <- function(df, tmb_col = "TMB", group_col = "risk_group") {
  if (!(tmb_col %in% colnames(df)) || !(group_col %in% colnames(df))) return(NULL)
  w <- stats::wilcox.test(df[[tmb_col]] ~ df[[group_col]])
  out <- data.frame(group = levels(df[[group_col]]), p_value = w$p.value)
  data.table::fwrite(out, file.path(config$output$report_dir, "tmb_riskgroup_wilcox.csv"))
  out
}

run_drug_sensitivity_stub <- function(expr_mat, model = "oncopredict") {
  # Placeholder to keep workflow reproducible while package/data requirements vary.
  # Replace with real calls to oncoPredict/pRRophetic according to local environment.
  msg <- paste0("Drug sensitivity step placeholder. Recommended model: ", model,
                ". Provide training pharmacogenomics reference and run prediction in your environment.")
  writeLines(msg, con = file.path(config$output$report_dir, "drug_sensitivity_note.txt"))
  invisible(msg)
}

run_extended_workflow <- function() {
  all_df <- readRDS(file.path(config$output$report_dir, "merged_subgroup_scored_data.rds"))

  # Identify candidate features from model artifact
  best <- readRDS(file.path(config$output$models_dir, "best_model_artifact.rds"))
  feat <- best$features[best$features %in% colnames(all_df)]

  surv_df <- prepare_survival_df(all_df, feature_cols = feat)
  uni <- run_univariate_cox(surv_df, feature_cols = feat)
  data.table::fwrite(uni, file.path(config$output$report_dir, "univariate_cox_signature.csv"))

  lasso <- run_lasso_cox(surv_df, feature_cols = feat)
  coef_tbl <- data.frame(feature = rownames(lasso$coef), coef = as.numeric(lasso$coef[, 1]))
  data.table::fwrite(coef_tbl, file.path(config$output$report_dir, "lasso_cox_coefficients.csv"))

  scored <- evaluate_survival_signature(surv_df, lasso$risk_score)
  saveRDS(scored, file.path(config$output$report_dir, "lasso_cox_scored_data.rds"))

  # Optional nomogram
  build_nomogram_if_available(scored)

  # Optional immune infiltration (requires expression matrix as genes x samples)
  expr <- t(as.matrix(all_df[, feat, drop = FALSE]))
  imm <- run_immune_infiltration(expr)
  if (!is.null(imm)) data.table::fwrite(as.data.frame(imm), file.path(config$output$report_dir, "immune_infiltration_scores.csv"))

  analyze_tmb_association(scored)
  run_drug_sensitivity_stub(expr)
}

if (interactive()) {
  run_extended_workflow()
}

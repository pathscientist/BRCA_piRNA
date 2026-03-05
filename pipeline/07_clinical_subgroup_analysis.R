# 07_clinical_subgroup_analysis.R
# Merge BRCA + yyfbatch1 + yyfbatch2, compute model score, do subgroup ROC,
# compare tumor vs normal T-score by subgroup, and perform Cox + forest plot
# using binary score variable (not continuous).

source("pipeline/00_config.R")

safe_stage_binary <- function(stage_vec) {
  x <- as.character(stage_vec)
  out <- ifelse(grepl("I$|II$|stage i|stage ii", tolower(x)), "Early", "Late")
  factor(out)
}

safe_age_binary <- function(age_vec) {
  med <- stats::median(age_vec, na.rm = TRUE)
  factor(ifelse(age_vec >= med, "HighAge", "LowAge"))
}

safe_subtype_binary <- function(subtype_vec) {
  x <- as.character(subtype_vec)
  top <- names(sort(table(x), decreasing = TRUE))[1]
  factor(ifelse(x == top, top, "Other"))
}

predict_score <- function(model_artifact, df) {
  feats <- model_artifact$features
  feats <- feats[feats %in% colnames(df)]
  probs <- predict(model_artifact$fit, newdata = df[, feats, drop = FALSE], type = "prob")[[config$input$positive_class]]
  score <- stats::qlogis(pmin(pmax(probs, 1e-6), 1 - 1e-6))
  list(prob = probs, t_score = as.numeric(score))
}

subgroup_roc <- function(df, subgroup_col, label_col, prob_col, out_prefix) {
  groups <- na.omit(unique(df[[subgroup_col]]))
  roc_rows <- list()
  curve <- list()

  for (g in groups) {
    tmp <- df[df[[subgroup_col]] == g, , drop = FALSE]
    if (length(unique(tmp[[label_col]])) < 2) next
    r <- pROC::roc(tmp[[label_col]], tmp[[prob_col]], levels = rev(levels(tmp[[label_col]])), quiet = TRUE)
    roc_rows[[as.character(g)]] <- data.frame(subgroup = as.character(g), AUROC = as.numeric(pROC::auc(r)), n = nrow(tmp))
    curve[[as.character(g)]] <- data.frame(subgroup = as.character(g), fpr = 1 - r$specificities, tpr = r$sensitivities)
  }

  roc_tbl <- dplyr::bind_rows(roc_rows)
  curve_tbl <- dplyr::bind_rows(curve)

  p <- ggplot2::ggplot(curve_tbl, ggplot2::aes(x = fpr, y = tpr, color = subgroup)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("Subgroup ROC: ", subgroup_col), x = "FPR", y = "TPR")

  ggplot2::ggsave(file.path(config$output$report_dir, paste0(out_prefix, "_", subgroup_col, "_roc.png")), p, width = 7, height = 5, dpi = 300)
  data.table::fwrite(roc_tbl, file.path(config$output$report_dir, paste0(out_prefix, "_", subgroup_col, "_auc.csv")))
}

plot_tscore_box <- function(df, subgroup_col, label_col, score_col, out_prefix) {
  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = subgroup_col, y = score_col, fill = label_col)) +
    ggplot2::geom_boxplot(outlier.alpha = 0.4) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("T-score by ", subgroup_col, " and label"), x = subgroup_col, y = "T-score")

  ggplot2::ggsave(file.path(config$output$report_dir, paste0(out_prefix, "_", subgroup_col, "_tscore_box.png")), p, width = 8, height = 5, dpi = 300)
}

run_cox_binary_score <- function(df, score_bin_col, time_col = "OS_time", event_col = "OS_event") {
  needed <- c(score_bin_col, time_col, event_col, "age_bin", "stage_bin", "subtype_bin")
  needed <- needed[needed %in% colnames(df)]
  cox_df <- df[, needed, drop = FALSE]
  cox_df <- cox_df[stats::complete.cases(cox_df), , drop = FALSE]
  if (nrow(cox_df) < 20) return(NULL)

  f <- stats::as.formula(paste0("survival::Surv(", time_col, ",", event_col, ") ~ ", paste(setdiff(needed, c(time_col, event_col)), collapse = " + ")))
  fit <- survival::coxph(f, data = cox_df)
  tbl <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)

  p <- ggplot2::ggplot(tbl, ggplot2::aes(x = reorder(term, estimate), y = estimate, ymin = conf.low, ymax = conf.high)) +
    ggplot2::geom_pointrange() +
    ggplot2::geom_hline(yintercept = 1, linetype = 2, color = "red") +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Cox regression forest plot (binary score)", x = "Covariate", y = "Hazard Ratio")

  data.table::fwrite(tbl, file.path(config$output$report_dir, "cox_binary_score_results.csv"))
  ggplot2::ggsave(file.path(config$output$report_dir, "cox_binary_score_forest.png"), p, width = 8, height = 5, dpi = 300)
  saveRDS(fit, file.path(config$output$report_dir, "cox_binary_score_model.rds"))
}

if (interactive()) {
  all_df <- readRDS(file.path(config$output$dir, "all_batch_corrected.rds"))
  all_df <- all_df[all_df[[config$input$batch_col]] %in% c("BRCA", "yyfbatch1", "yyfbatch2"), , drop = FALSE]
  all_df[[config$input$label_col]] <- factorize_label(all_df[[config$input$label_col]], config$input$positive_class, config$input$negative_class)

  best_model <- readRDS(file.path(config$output$models_dir, "best_model_artifact.rds"))
  sc <- predict_score(best_model, all_df)
  all_df$model_prob <- sc$prob
  all_df$T_score <- sc$t_score
  all_df$score_bin <- factor(ifelse(all_df$model_prob >= best_model$threshold$threshold, "HighScore", "LowScore"))

  if ("age" %in% colnames(all_df)) all_df$age_bin <- safe_age_binary(all_df$age)
  if ("stage" %in% colnames(all_df)) all_df$stage_bin <- safe_stage_binary(all_df$stage)
  if ("subtype" %in% colnames(all_df)) all_df$subtype_bin <- safe_subtype_binary(all_df$subtype)

  subgroup_roc(all_df, "age_bin", config$input$label_col, "model_prob", "merged")
  subgroup_roc(all_df, "stage_bin", config$input$label_col, "model_prob", "merged")
  subgroup_roc(all_df, "subtype_bin", config$input$label_col, "model_prob", "merged")

  plot_tscore_box(all_df, "age_bin", config$input$label_col, "T_score", "merged")
  plot_tscore_box(all_df, "stage_bin", config$input$label_col, "T_score", "merged")
  plot_tscore_box(all_df, "subtype_bin", config$input$label_col, "T_score", "merged")

  run_cox_binary_score(all_df, score_bin_col = "score_bin", time_col = "OS_time", event_col = "OS_event")
  saveRDS(all_df, file.path(config$output$report_dir, "merged_subgroup_scored_data.rds"))
}

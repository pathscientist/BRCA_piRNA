# 06_reports_plots.R
# Generate summary report artifacts and key plots.

source("pipeline/00_config.R")

plot_roc_pr <- function(detail_obj, prefix) {
  roc_obj <- pROC::roc(response = detail_obj$truth, predictor = detail_obj$prob, levels = rev(levels(detail_obj$truth)), quiet = TRUE)
  roc_df <- data.frame(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)

  p_roc <- ggplot2::ggplot(roc_df, ggplot2::aes(x = fpr, y = tpr)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
    ggplot2::labs(title = paste("ROC -", prefix), x = "False Positive Rate", y = "True Positive Rate") +
    ggplot2::theme_minimal()

  ggplot2::ggsave(filename = file.path(config$output$report_dir, paste0(prefix, "_roc.png")), plot = p_roc, width = 6, height = 5, dpi = 300)
}

plot_model_ranking <- function(summary_tbl) {
  p <- ggplot2::ggplot(summary_tbl, ggplot2::aes(x = reorder(paste(feature_set, model, sep = "|"), AUROC), y = AUROC, fill = n_features)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Feature Set | Model", y = "Cross-validated AUROC", title = "Internal model ranking")

  ggplot2::ggsave(file.path(config$output$report_dir, "internal_model_ranking.png"), p, width = 8, height = 10, dpi = 300)
}

write_summary <- function() {
  best <- readRDS(file.path(config$output$models_dir, "best_model_row.rds"))
  val <- data.table::fread(file.path(config$output$report_dir, "validation_metrics.csv"), data.table = FALSE)

  lines <- c(
    "# piRNA Diagnosis Pipeline Summary",
    "",
    paste("Best feature set:", best$feature_set),
    paste("Best model:", best$model),
    paste("Features used:", best$n_features),
    paste("Tuned threshold:", round(best$threshold, 4)),
    "",
    "## Validation metrics",
    capture.output(print(val)),
    "",
    "## sessionInfo()",
    capture.output(sessionInfo())
  )

  writeLines(lines, con = file.path(config$output$report_dir, "final_summary_report.md"))
}

if (interactive()) {
  summary_tbl <- data.table::fread(file.path(config$output$models_dir, "internal_cv_summary.csv"), data.table = FALSE)
  holdout_detail <- readRDS(file.path(config$output$report_dir, "holdout_validation_detail.rds"))

  plot_model_ranking(summary_tbl)
  plot_roc_pr(holdout_detail, "internal_holdout")

  ext_path <- file.path(config$output$report_dir, "external_validation_detail.rds")
  if (file.exists(ext_path)) {
    ext_detail <- readRDS(ext_path)
    plot_roc_pr(ext_detail, "external")
  }

  write_summary()
}

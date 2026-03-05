# 04_model_selection.R
# Select best internal model with AUROC primary, AUPRC secondary, simplicity tie-break.

source("pipeline/00_config.R")

pick_best_model <- function(summary_tbl, delta_auc = config$feature$near_best_delta_auc) {
  if (nrow(summary_tbl) == 0) stop("No model results available.")

  best_auc <- max(summary_tbl$AUROC, na.rm = TRUE)
  candidates <- summary_tbl[summary_tbl$AUROC >= (best_auc - delta_auc), , drop = FALSE]

  candidates <- candidates[order(-candidates$AUROC, -candidates$AUPRC, candidates$n_features), , drop = FALSE]
  candidates[1, , drop = FALSE]
}

load_model_artifact <- function(feature_set, model) {
  key <- paste(feature_set, model, sep = "__")
  path <- file.path(config$output$models_dir, paste0(key, ".rds"))
  readRDS(path)
}

save_selected_model <- function(best_row) {
  mdl <- load_model_artifact(best_row$feature_set, best_row$model)
  saveRDS(best_row, file.path(config$output$models_dir, "best_model_row.rds"))
  saveRDS(mdl, file.path(config$output$models_dir, "best_model_artifact.rds"))
  data.table::fwrite(best_row, file.path(config$output$models_dir, "best_model_row.csv"))
}

if (interactive()) {
  summary_tbl <- data.table::fread(file.path(config$output$models_dir, "internal_cv_summary.csv"), data.table = FALSE)
  best <- pick_best_model(summary_tbl)
  save_selected_model(best)
}

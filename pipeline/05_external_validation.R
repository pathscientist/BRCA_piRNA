# 05_external_validation.R
# Lock preprocessing + features + model + threshold; run held-out and external validation.

source("pipeline/00_config.R")

predict_with_artifact <- function(model_artifact, df, sample_id_col, label_col) {
  feats <- model_artifact$features
  feats <- feats[feats %in% colnames(df)]
  if (length(feats) < 2) stop("Too few features found in validation data")

  x <- df[, feats, drop = FALSE]
  y <- df[[label_col]]
  probs <- predict(model_artifact$fit, newdata = x, type = "prob")[[config$input$positive_class]]
  list(truth = y, prob = probs)
}

validate_dataset <- function(model_artifact, df, dataset_name) {
  pred <- predict_with_artifact(model_artifact, df, config$input$sample_id_col, config$input$label_col)
  m <- calc_binary_metrics(
    truth = pred$truth,
    prob = pred$prob,
    positive = config$input$positive_class,
    threshold = model_artifact$threshold$threshold
  )

  row <- dplyr::mutate(m$metrics, dataset = dataset_name, threshold = model_artifact$threshold$threshold)
  list(row = row, confusion = m$confusion, prob = pred$prob, truth = pred$truth)
}

if (interactive()) {
  best_model <- readRDS(file.path(config$output$models_dir, "best_model_artifact.rds"))
  test_df <- readRDS(file.path(config$output$dir, "test_data.rds"))
  external_df <- readRDS(file.path(config$output$dir, "external_data.rds"))

  test_val <- validate_dataset(best_model, test_df, "internal_holdout")
  out_rows <- list(test_val$row)

  if (!is.null(external_df) && nrow(external_df) > 0) {
    ext_val <- validate_dataset(best_model, external_df, "external")
    out_rows <- append(out_rows, list(ext_val$row))
    saveRDS(ext_val, file.path(config$output$report_dir, "external_validation_detail.rds"))
  }

  metrics_tbl <- dplyr::bind_rows(out_rows)
  data.table::fwrite(metrics_tbl, file.path(config$output$report_dir, "validation_metrics.csv"))
  saveRDS(test_val, file.path(config$output$report_dir, "holdout_validation_detail.rds"))
}

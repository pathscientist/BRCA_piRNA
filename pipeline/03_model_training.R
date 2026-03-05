# 03_model_training.R
# Train and tune multiple classifiers over candidate feature sets (internal CV only).

source("pipeline/00_config.R")

prepare_xy <- function(df, sample_id_col, label_col, features) {
  x <- df[, features, drop = FALSE]
  y <- df[[label_col]]
  list(x = x, y = y)
}

train_control <- function() {
  caret::trainControl(
    method = "repeatedcv",
    number = config$split$outer_folds,
    repeats = config$split$outer_repeats,
    classProbs = TRUE,
    summaryFunction = caret::twoClassSummary,
    savePredictions = "final"
  )
}

train_one_model <- function(x, y, model_name) {
  ctrl <- train_control()

  if (model_name == "logistic") {
    return(caret::train(x = x, y = y, method = "glm", family = binomial(), trControl = ctrl, metric = "ROC"))
  }

  if (model_name == "elastic_net") {
    grid <- expand.grid(alpha = seq(0, 1, by = 0.25), lambda = 10 ^ seq(-4, -1, length.out = 10))
    return(caret::train(x = x, y = y, method = "glmnet", trControl = ctrl, tuneGrid = grid, preProcess = c("center", "scale"), metric = "ROC"))
  }

  if (model_name == "rf") {
    grid <- expand.grid(mtry = unique(pmax(1, floor(c(sqrt(ncol(x)), ncol(x) / 3, ncol(x) / 2)))))
    return(caret::train(x = x, y = y, method = "rf", trControl = ctrl, tuneGrid = grid, ntree = 500, metric = "ROC"))
  }

  if (model_name == "xgb") {
    grid <- expand.grid(
      nrounds = c(100, 200),
      max_depth = c(2, 4, 6),
      eta = c(0.03, 0.1),
      gamma = c(0, 1),
      colsample_bytree = c(0.6, 0.8),
      min_child_weight = c(1, 5),
      subsample = c(0.7, 0.9)
    )
    return(caret::train(x = x, y = y, method = "xgbTree", trControl = ctrl, tuneGrid = grid, metric = "ROC", verbosity = 0))
  }

  if (model_name == "svm_rbf") {
    grid <- expand.grid(C = 2 ^ (-1:3), sigma = 2 ^ (-7:-3))
    return(caret::train(x = x, y = y, method = "svmRadial", trControl = ctrl, tuneGrid = grid, preProcess = c("center", "scale"), metric = "ROC"))
  }

  stop("Unknown model: ", model_name)
}

run_model_grid <- function(train_df, sample_id_col, label_col, candidate_sets) {
  model_names <- c("logistic", "elastic_net", "rf", "xgb", "svm_rbf")
  results <- list()
  summary_rows <- list()

  for (set_name in names(candidate_sets)) {
    feats <- unique(candidate_sets[[set_name]])
    feats <- feats[feats %in% colnames(train_df)]
    feats <- feats[seq_len(min(length(feats), config$feature$max_features))]
    if (length(feats) < 2) next

    xy <- prepare_xy(train_df, sample_id_col, label_col, feats)

    for (mdl in model_names) {
      key <- paste(set_name, mdl, sep = "__")
      message("Training: ", key)
      fit <- train_one_model(xy$x, xy$y, mdl)
      pred <- fit$pred
      pos <- config$input$positive_class
      prob <- pred[[pos]]
      truth <- pred$obs

      th <- optimize_threshold(truth, prob, metric = "youden", positive = pos)
      m <- calc_binary_metrics(truth, prob, positive = pos, threshold = th$threshold)

      results[[key]] <- list(
        feature_set = set_name,
        model = mdl,
        features = feats,
        fit = fit,
        threshold = th,
        cv_metrics = m$metrics
      )

      summary_rows[[key]] <- dplyr::mutate(m$metrics, feature_set = set_name, model = mdl, n_features = length(feats), threshold = th$threshold)
      saveRDS(results[[key]], file.path(config$output$models_dir, paste0(key, ".rds")))
    }
  }

  summary_tbl <- dplyr::bind_rows(summary_rows)
  data.table::fwrite(summary_tbl, file.path(config$output$models_dir, "internal_cv_summary.csv"))
  list(results = results, summary = summary_tbl)
}

if (interactive()) {
  train_df <- readRDS(file.path(config$output$dir, "train_data.rds"))
  candidate_sets <- readRDS(file.path(config$output$fs_dir, "candidate_sets.rds"))
  out <- run_model_grid(train_df, config$input$sample_id_col, config$input$label_col, candidate_sets)
  saveRDS(out, file.path(config$output$models_dir, "all_models_internal.rds"))
}

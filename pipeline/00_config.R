# 00_config.R
# Global configuration and helper utilities for miRNA TPM diagnosis pipeline.

set.seed(2026)

config <- list(
  input = list(
    matrix_path = "<PATH_TO_MATRIX>",
    labels_path = "<PATH_TO_LABELS>",
    sample_id_col = "<SAMPLE_ID_COL>",
    label_col = "<LABEL_COL>",
    positive_class = "Cancer",
    negative_class = "Normal",
    batch_col = "batch",
    cohort_col = "cohort",
    validation_batch = "yyfbatch1", # choose yyfbatch1 or yyfbatch2
    optional_covariates = c("age", "stage", "subtype", "OS_time", "OS_event")
  ),
  split = list(
    train_prop = 0.8,
    outer_folds = 5,
    outer_repeats = 2,
    inner_folds = 5
  ),
  balance = list(
    brca_name = "BRCA1",
    add_excess_fraction = 0.4
  ),
  feature = list(
    top_n = c(5, 8, 9),
    min_selected_methods = 3,
    near_best_delta_auc = 0.01,
    max_features = 9
  ),
  output = list(
    dir = "artifacts",
    fs_dir = "artifacts/feature_selection",
    models_dir = "artifacts/models",
    report_dir = "artifacts/reports"
  )
)

required_packages <- c(
  "data.table", "dplyr", "caret", "glmnet", "pROC", "PRROC",
  "randomForest", "xgboost", "e1071", "Boruta", "mRMRe",
  "ggplot2", "purrr", "tibble", "survival", "broom"
)

check_packages <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    message("Missing packages (install before full run): ", paste(missing_pkgs, collapse = ", "))
  }
}

ensure_dirs <- function(paths) {
  for (p in paths) {
    if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
}

factorize_label <- function(x, pos = "Cancer", neg = "Normal") {
  factor(x, levels = c(neg, pos))
}

calc_binary_metrics <- function(truth, prob, positive = "Cancer", threshold = 0.5) {
  pred <- factor(ifelse(prob >= threshold, positive, setdiff(levels(truth), positive)[1]), levels = levels(truth))
  cm <- caret::confusionMatrix(pred, truth, positive = positive)
  precision <- cm$byClass[["Pos Pred Value"]]
  recall <- cm$byClass[["Sensitivity"]]
  f1 <- ifelse((precision + recall) == 0, 0, 2 * precision * recall / (precision + recall))
  mcc_num <- (cm$table[2, 2] * cm$table[1, 1]) - (cm$table[1, 2] * cm$table[2, 1])
  mcc_den <- sqrt(prod(rowSums(cm$table)) * prod(colSums(cm$table)))
  mcc <- ifelse(mcc_den == 0, 0, mcc_num / mcc_den)

  roc_obj <- pROC::roc(response = truth, predictor = prob, levels = rev(levels(truth)), quiet = TRUE)
  pr_obj <- PRROC::pr.curve(scores.class0 = prob[truth == positive], scores.class1 = prob[truth != positive], curve = FALSE)

  list(
    confusion = cm,
    metrics = tibble::tibble(
      AUROC = as.numeric(pROC::auc(roc_obj)),
      AUPRC = as.numeric(pr_obj$auc.integral),
      Accuracy = as.numeric(cm$overall[["Accuracy"]]),
      Sensitivity = as.numeric(cm$byClass[["Sensitivity"]]),
      Specificity = as.numeric(cm$byClass[["Specificity"]]),
      F1 = as.numeric(f1),
      MCC = as.numeric(mcc)
    )
  )
}

optimize_threshold <- function(truth, prob, metric = c("youden", "f1"), positive = "Cancer") {
  metric <- match.arg(metric)
  thresholds <- sort(unique(prob))
  best_score <- -Inf
  best_t <- 0.5

  for (t in thresholds) {
    pred <- factor(ifelse(prob >= t, positive, setdiff(levels(truth), positive)[1]), levels = levels(truth))
    cm <- caret::confusionMatrix(pred, truth, positive = positive)
    sens <- cm$byClass[["Sensitivity"]]
    spec <- cm$byClass[["Specificity"]]
    precision <- cm$byClass[["Pos Pred Value"]]
    f1 <- ifelse((precision + sens) == 0, 0, 2 * precision * sens / (precision + sens))
    score <- if (metric == "youden") sens + spec - 1 else f1
    if (!is.na(score) && score > best_score) {
      best_score <- score
      best_t <- t
    }
  }

  list(threshold = best_t, score = best_score, objective = metric)
}

check_packages(required_packages)
ensure_dirs(unlist(config$output))

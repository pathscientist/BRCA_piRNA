# 02_feature_selection.R
# Multi-method feature selection (>=8) on training data only.

source("pipeline/00_config.R")

get_xy <- function(df, sample_id_col, label_col) {
  feature_cols <- setdiff(colnames(df), c(sample_id_col, label_col, config$input$optional_covariates))
  x <- as.data.frame(df[, feature_cols, drop = FALSE])
  y <- df[[label_col]]
  list(x = x, y = y, feature_cols = feature_cols)
}

rank_de <- function(x, y) {
  pvals <- vapply(x, function(col) {
    suppressWarnings(wilcox.test(col ~ y)$p.value)
  }, numeric(1))
  padj <- p.adjust(pvals, method = "fdr")
  tibble::tibble(feature = names(pvals), score = -log10(padj + 1e-12), selected = padj < 0.05)
}

rank_lasso <- function(x, y, alpha = 1) {
  xmat <- as.matrix(x)
  ybin <- ifelse(y == levels(y)[2], 1, 0)
  fit <- glmnet::cv.glmnet(xmat, ybin, family = "binomial", alpha = alpha, nfolds = 5)
  coefs <- glmnet::coef.glmnet(fit, s = "lambda.1se")
  vals <- as.matrix(coefs)
  feats <- rownames(vals)[-1]
  coef_vals <- abs(vals[-1, 1])
  tibble::tibble(feature = feats, score = coef_vals, selected = coef_vals > 0)
}

rank_mrmr <- function(x, y, top_k = 30) {
  if (!requireNamespace("mRMRe", quietly = TRUE)) {
    return(tibble::tibble(feature = colnames(x), score = NA_real_, selected = FALSE, method_warning = "mRMRe missing"))
  }
  dat <- data.frame(label = as.numeric(y == levels(y)[2]), x)
  dd <- mRMRe::mRMR.data(data = dat)
  fs <- mRMRe::mRMR.classic(data = dd, target_indices = 1, feature_count = min(top_k, ncol(x)))
  idx <- unlist(mRMRe::solutions(fs))
  feats <- colnames(dat)[idx]
  tibble::tibble(feature = colnames(x), score = as.numeric(colnames(x) %in% feats), selected = colnames(x) %in% feats)
}

rank_boruta <- function(x, y) {
  if (!requireNamespace("Boruta", quietly = TRUE)) {
    return(tibble::tibble(feature = colnames(x), score = NA_real_, selected = FALSE, method_warning = "Boruta missing"))
  }
  bor <- Boruta::Boruta(x = x, y = y, doTrace = 0, maxRuns = 100)
  stats <- Boruta::attStats(bor)
  tibble::tibble(feature = rownames(stats), score = stats$meanImp, selected = stats$decision %in% c("Confirmed", "Tentative"))
}

rank_rf <- function(x, y) {
  fit <- randomForest::randomForest(x = x, y = y, importance = TRUE, ntree = 500)
  imp <- randomForest::importance(fit)[, "MeanDecreaseGini"]
  tibble::tibble(feature = names(imp), score = as.numeric(imp), selected = imp > stats::median(imp, na.rm = TRUE))
}

rank_xgb <- function(x, y) {
  ybin <- ifelse(y == levels(y)[2], 1, 0)
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(x), label = ybin)
  fit <- xgboost::xgb.train(
    params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = 3, eta = 0.1),
    data = dtrain,
    nrounds = 100,
    verbose = 0
  )
  imp <- xgboost::xgb.importance(model = fit)
  out <- tibble::tibble(feature = colnames(x), score = 0)
  out$score[match(imp$Feature, out$feature)] <- imp$Gain
  tibble::add_column(out, selected = out$score > stats::median(out$score, na.rm = TRUE), .after = "score")
}

rank_svm_rfe <- function(x, y) {
  ctrl <- caret::rfeControl(functions = caret::caretFuncs, method = "cv", number = 5)
  sizes <- unique(pmin(c(5, 10, 20, 30), ncol(x)))
  fit <- caret::rfe(x = x, y = y, sizes = sizes, rfeControl = ctrl, method = "svmRadial")
  prof <- fit$variables
  agg <- prof %>% dplyr::group_by(var) %>% dplyr::summarise(score = max(Overall), .groups = "drop")
  tibble::tibble(feature = colnames(x), score = agg$score[match(colnames(x), agg$var)] %||% 0) |>
    dplyr::mutate(score = ifelse(is.na(score), 0, score), selected = feature %in% fit$optVariables)
}

rank_mutual_information <- function(x, y) {
  if (!requireNamespace("FSelectorRcpp", quietly = TRUE)) {
    return(tibble::tibble(feature = colnames(x), score = NA_real_, selected = FALSE, method_warning = "FSelectorRcpp missing"))
  }
  dat <- data.frame(y = y, x)
  mi <- FSelectorRcpp::information_gain(y ~ ., dat)
  tibble::tibble(feature = rownames(mi), score = mi$importance, selected = mi$importance > stats::median(mi$importance, na.rm = TRUE))
}

`%||%` <- function(a, b) if (length(a) == 0) b else a

run_feature_selection <- function(train_df, sample_id_col, label_col) {
  xy <- get_xy(train_df, sample_id_col, label_col)
  x <- xy$x
  y <- xy$y

  methods <- list(
    DE = rank_de(x, y),
    LASSO = rank_lasso(x, y, alpha = 1),
    ElasticNet = rank_lasso(x, y, alpha = 0.5),
    mRMR = rank_mrmr(x, y),
    Boruta = rank_boruta(x, y),
    RF = rank_rf(x, y),
    XGBoost = rank_xgb(x, y),
    SVM_RFE = rank_svm_rfe(x, y),
    MutualInfo = rank_mutual_information(x, y)
  )

  purrr::iwalk(methods, function(tbl, nm) {
    data.table::fwrite(tbl, file.path(config$output$fs_dir, paste0("rank_", nm, ".csv")))
  })

  selected_list <- lapply(methods, function(tbl) tbl$feature[which(tbl$selected %in% TRUE)])
  freq <- sort(table(unlist(selected_list)), decreasing = TRUE)

  intersection_set <- Reduce(intersect, selected_list)
  union_freq_set <- names(freq[freq >= config$feature$min_selected_methods])

  trim_set <- function(v) unique(v)[seq_len(min(length(unique(v)), config$feature$max_features))]
  intersection_set <- trim_set(intersection_set)
  union_freq_set <- trim_set(union_freq_set)

  topn_sets <- lapply(config$feature$top_n, function(k) names(freq)[seq_len(min(k, length(freq), config$feature$max_features))])
  names(topn_sets) <- paste0("top_", config$feature$top_n)

  saveRDS(selected_list, file.path(config$output$fs_dir, "selected_by_method.rds"))
  saveRDS(freq, file.path(config$output$fs_dir, "selection_frequency.rds"))

  list(
    methods = methods,
    selected_list = selected_list,
    candidate_sets = c(
      list(intersection_all = intersection_set, freq_union = union_freq_set),
      topn_sets
    )
  )
}

if (interactive()) {
  train_df <- readRDS(file.path(config$output$dir, "train_data.rds"))
  fs <- run_feature_selection(train_df, config$input$sample_id_col, config$input$label_col)
  saveRDS(fs$candidate_sets, file.path(config$output$fs_dir, "candidate_sets.rds"))
}

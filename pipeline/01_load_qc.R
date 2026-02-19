# 01_load_qc.R
# Load TPM matrix + labels, perform global batch correction, define locked validation batch,
# and optionally rebalance BRCA1 class ratio in training data.

source("pipeline/00_config.R")

read_expression_and_labels <- function(matrix_path, labels_path, sample_id_col, label_col) {
  expr <- data.table::fread(matrix_path, data.table = FALSE)
  labels <- data.table::fread(labels_path, data.table = FALSE)

  if (!(sample_id_col %in% colnames(labels))) stop("sample_id_col not found in labels")
  if (!(label_col %in% colnames(labels))) stop("label_col not found in labels")

  if (sample_id_col %in% colnames(expr)) {
    expr_df <- expr
  } else {
    rownames(expr) <- expr[[1]]
    expr_df <- expr[, -1, drop = FALSE]
    expr_df[[sample_id_col]] <- rownames(expr)
    expr_df <- expr_df[, c(sample_id_col, setdiff(colnames(expr_df), sample_id_col)), drop = FALSE]
  }

  merged <- dplyr::inner_join(labels, expr_df, by = sample_id_col)
  merged[[label_col]] <- factorize_label(merged[[label_col]], config$input$positive_class, config$input$negative_class)
  stopifnot(all(!is.na(merged[[label_col]])))
  merged
}

get_feature_cols <- function(df, sample_id_col, label_col) {
  setdiff(colnames(df), c(sample_id_col, label_col, config$input$optional_covariates, config$input$batch_col, config$input$cohort_col))
}

remove_batch_effect_all <- function(df, sample_id_col, label_col, batch_col) {
  if (!(batch_col %in% colnames(df))) {
    warning("batch_col not found; skipping batch correction")
    return(df)
  }

  feature_cols <- get_feature_cols(df, sample_id_col, label_col)
  x <- as.matrix(df[, feature_cols, drop = FALSE])
  batch <- as.factor(df[[batch_col]])

  corrected <- NULL
  if (requireNamespace("sva", quietly = TRUE)) {
    corrected <- t(sva::ComBat(dat = t(x), batch = batch, par.prior = TRUE, prior.plots = FALSE))
  } else if (requireNamespace("limma", quietly = TRUE)) {
    corrected <- limma::removeBatchEffect(x, batch = batch)
  } else {
    warning("Neither sva nor limma available; batch correction skipped")
    return(df)
  }

  df[, feature_cols] <- as.data.frame(corrected)
  df
}

rebalance_brca_training <- function(train_df, cohort_col, label_col, brca_name, add_excess_fraction = 0.4) {
  if (!(cohort_col %in% colnames(train_df))) return(train_df)

  brca <- train_df[train_df[[cohort_col]] == brca_name, , drop = FALSE]
  other <- train_df[train_df[[cohort_col]] != brca_name, , drop = FALSE]
  if (nrow(brca) == 0) return(train_df)

  tumor <- brca[brca[[label_col]] == config$input$positive_class, , drop = FALSE]
  normal <- brca[brca[[label_col]] == config$input$negative_class, , drop = FALSE]

  if (nrow(tumor) == 0 || nrow(normal) == 0) return(train_df)

  n_pairs <- min(nrow(tumor), nrow(normal))
  major_class <- if (nrow(tumor) > nrow(normal)) config$input$positive_class else config$input$negative_class

  major_df <- brca[brca[[label_col]] == major_class, , drop = FALSE]
  minor_df <- brca[brca[[label_col]] != major_class, , drop = FALSE]

  major_pair_idx <- sample(seq_len(nrow(major_df)), n_pairs)
  minor_pair_idx <- sample(seq_len(nrow(minor_df)), n_pairs)
  major_pair <- major_df[major_pair_idx, , drop = FALSE]
  minor_pair <- minor_df[minor_pair_idx, , drop = FALSE]

  extra_n <- floor((nrow(major_df) - n_pairs) * add_excess_fraction)
  extra_df <- major_df[setdiff(seq_len(nrow(major_df)), major_pair_idx), , drop = FALSE]
  if (extra_n > 0 && nrow(extra_df) > 0) {
    extra_n <- min(extra_n, nrow(extra_df))
    extra_df <- extra_df[sample(seq_len(nrow(extra_df)), extra_n), , drop = FALSE]
  } else {
    extra_df <- extra_df[0, , drop = FALSE]
  }

  balanced_brca <- rbind(major_pair, minor_pair, extra_df)
  out <- rbind(other, balanced_brca)
  out[sample(seq_len(nrow(out))), , drop = FALSE]
}

split_train_and_validation_batch <- function(df, batch_col, validation_batch) {
  if (!(batch_col %in% colnames(df))) stop("batch_col missing, cannot split validation batch")
  ext <- df[df[[batch_col]] == validation_batch, , drop = FALSE]
  dev <- df[df[[batch_col]] != validation_batch, , drop = FALSE]
  list(dev = dev, external = ext)
}

split_dev_train_test <- function(dev_df, label_col) {
  idx <- caret::createDataPartition(dev_df[[label_col]], p = config$split$train_prop, list = FALSE)
  list(train = dev_df[idx, , drop = FALSE], test = dev_df[-idx, , drop = FALSE])
}

save_split_objects <- function(split_obj) {
  saveRDS(split_obj$train, file.path(config$output$dir, "train_data.rds"))
  saveRDS(split_obj$test, file.path(config$output$dir, "test_data.rds"))
  saveRDS(split_obj$external, file.path(config$output$dir, "external_data.rds"))
  saveRDS(split_obj$all_corrected, file.path(config$output$dir, "all_batch_corrected.rds"))
}

if (interactive()) {
  dat <- read_expression_and_labels(
    config$input$matrix_path,
    config$input$labels_path,
    config$input$sample_id_col,
    config$input$label_col
  )
  dat <- remove_batch_effect_all(dat, config$input$sample_id_col, config$input$label_col, config$input$batch_col)

  split1 <- split_train_and_validation_batch(dat, config$input$batch_col, config$input$validation_batch)
  split2 <- split_dev_train_test(split1$dev, config$input$label_col)

  train_bal <- rebalance_brca_training(
    split2$train,
    config$input$cohort_col,
    config$input$label_col,
    config$balance$brca_name,
    config$balance$add_excess_fraction
  )

  save_split_objects(list(train = train_bal, test = split2$test, external = split1$external, all_corrected = dat))
}

################################################################################
#                                                                              #
#   Breast Cancer miRNA Diagnostic Pipeline                                    #
#   Complete Feature Selection → Model Building → Independent Validation       #
#                                                                              #
#   Phase 1: 8 feature selection methods → consensus miRNA set                 #
#   Phase 2: 6 ML models with hyperparameter tuning → best model              #
#   Phase 3: Independent validation with robustness checks                     #
#                                                                              #
################################################################################

start_time <- Sys.time()

# ==============================================================================
# PHASE 0: SETUP & PREPROCESSING
# ==============================================================================

# --- 0.1 Install/load packages ------------------------------------------------
required_packages <- c(
  "caret", "randomForest", "e1071", "glmnet", "xgboost", "nnet",
  "pROC", "limma", "Boruta", "UpSetR", "pheatmap",
  "ggplot2", "dplyr", "tidyr", "gridExtra", "viridis", "scales",
  "smotefamily", "praznik", "PRROC"
)

bioc_packages <- c("limma")

# Install CRAN packages if missing
for (pkg in setdiff(required_packages, bioc_packages)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, quiet = TRUE)
  }
}

# Install Bioconductor packages if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", quiet = TRUE)
}
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# Load all
suppressPackageStartupMessages({
  library(caret)
  library(randomForest)
  library(e1071)
  library(glmnet)
  library(xgboost)
  library(nnet)
  library(pROC)
  library(limma)
  library(Boruta)
  library(UpSetR)
  library(pheatmap)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(viridis)
  library(scales)
  library(smotefamily)
  library(praznik)
  library(PRROC)
})

cat("All packages loaded successfully.\n")

# --- 0.2 Global settings -----------------------------------------------------
SEED <- 2024
COLOR_PALETTE <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928")

theme_set(theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey40"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ))

# Create output directories
dir.create("results/feature_selection", recursive = TRUE, showWarnings = FALSE)
dir.create("results/models", recursive = TRUE, showWarnings = FALSE)
dir.create("results/validation", recursive = TRUE, showWarnings = FALSE)

# --- 0.3 Load data ------------------------------------------------------------
cat("\n===== LOADING DATA =====\n")

# >>> EDIT THIS LINE to point to your CSV file <<<
# mirna_data <- read.csv("your_mirna_tpm_data.csv", row.names = 1, check.names = FALSE)
# Expected: rows = samples, columns = miRNAs + "label" column

# For demonstration, check if mirna_data already exists in environment
if (!exists("mirna_data")) {
  stop("Please load your data as 'mirna_data' before running this script.\n",
       "  Example: mirna_data <- read.csv('your_file.csv', row.names = 1)\n",
       "  Rows = samples, Columns = miRNAs + 'label' column")
}

mirna_data$label <- as.factor(mirna_data$label)

cat("Raw data dimensions:", nrow(mirna_data), "samples x", ncol(mirna_data) - 1, "miRNAs\n")
cat("Class distribution:\n")
print(table(mirna_data$label))

class_ratio <- max(table(mirna_data$label)) / min(table(mirna_data$label))
cat("Class ratio:", round(class_ratio, 2), "\n")
if (class_ratio > 2) {
  cat("WARNING: Classes are imbalanced (ratio > 2:1). SMOTE will be applied during modeling.\n")
  IMBALANCED <- TRUE
} else {
  IMBALANCED <- FALSE
}

# --- 0.4 Preprocessing -------------------------------------------------------
cat("\n===== PREPROCESSING =====\n")

# Separate label from features
labels <- mirna_data$label
expr_mat <- mirna_data[, colnames(mirna_data) != "label"]

# Ensure all expression columns are numeric
expr_mat <- as.data.frame(lapply(expr_mat, as.numeric))
rownames(expr_mat) <- rownames(mirna_data)

# Log2 transform TPM: log2(TPM + 1)
cat("Applying log2(TPM + 1) transformation...\n")
expr_mat <- log2(expr_mat + 1)

# Remove miRNAs with zero expression in >80% of samples
zero_frac <- colMeans(expr_mat == 0)
keep_nonzero <- zero_frac <= 0.80
expr_mat <- expr_mat[, keep_nonzero]
cat("Removed", sum(!keep_nonzero), "miRNAs with >80% zero expression.",
    ncol(expr_mat), "remaining.\n")

# Remove near-zero variance features
nzv <- nearZeroVar(expr_mat)
if (length(nzv) > 0) {
  expr_mat <- expr_mat[, -nzv]
  cat("Removed", length(nzv), "near-zero variance miRNAs.", ncol(expr_mat), "remaining.\n")
} else {
  cat("No near-zero variance miRNAs found.\n")
}

# Handle missing values
if (any(is.na(expr_mat))) {
  cat("Imputing missing values with column medians...\n")
  for (j in seq_len(ncol(expr_mat))) {
    na_idx <- is.na(expr_mat[, j])
    if (any(na_idx)) {
      expr_mat[na_idx, j] <- median(expr_mat[, j], na.rm = TRUE)
    }
  }
}

# Recombine
mirna_clean <- cbind(expr_mat, label = labels)

cat("Preprocessed data:", nrow(mirna_clean), "samples x", ncol(expr_mat), "miRNAs\n")

# --- 0.5 Train / Validation Split (70/30, stratified) ------------------------
cat("\n===== SPLITTING DATA =====\n")

set.seed(SEED)
train_idx <- createDataPartition(mirna_clean$label, p = 0.7, list = FALSE)
discovery_data <- mirna_clean[train_idx, ]
validation_data <- mirna_clean[-train_idx, ]

disc_labels <- discovery_data$label
disc_features <- discovery_data[, colnames(discovery_data) != "label"]

val_labels <- validation_data$label
val_features <- validation_data[, colnames(validation_data) != "label"]

cat("Discovery set:", nrow(discovery_data), "samples\n")
cat("Validation set:", nrow(validation_data), "samples (LOCKED until Phase 3)\n")
cat("Discovery class distribution:\n")
print(table(disc_labels))

# Compute and save preprocessing parameters (center + scale) from discovery set
preproc_params <- preProcess(disc_features, method = c("center", "scale"))
disc_features_scaled <- predict(preproc_params, disc_features)

cat("\nPreprocessing complete. Starting feature selection.\n")


# ==============================================================================
# PHASE 1: COMPREHENSIVE FEATURE SELECTION (8 METHODS)
# ==============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  PHASE 1: FEATURE SELECTION (8 METHODS)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Storage for selected features from each method
fs_results <- list()

# --- 1.1 Differential Expression Analysis (limma) ----------------------------
cat("--- 1.1 limma Differential Expression ---\n")
tryCatch({
  design <- model.matrix(~ disc_labels)
  fit <- lmFit(t(disc_features), design)
  fit <- eBayes(fit)
  limma_res <- topTable(fit, coef = 2, number = Inf, sort.by = "p")
  limma_res$miRNA <- rownames(limma_res)

  # Selection criteria: adj.P.Val < 0.01 and |logFC| > 1
  fs_limma <- limma_res$miRNA[limma_res$adj.P.Val < 0.01 & abs(limma_res$logFC) > 1]
  fs_results$limma <- fs_limma
  cat("  limma selected:", length(fs_limma), "miRNAs\n")

  # Volcano plot
  limma_res$sig <- ifelse(limma_res$adj.P.Val < 0.01 & abs(limma_res$logFC) > 1, "Significant", "NS")
  p_volcano <- ggplot(limma_res, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Significant" = "#E31A1C", "NS" = "grey60")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40") +
    labs(title = "Volcano Plot (limma)", x = "log2 Fold Change", y = "-log10(adj. P-value)") +
    theme(legend.title = element_blank())
  ggsave("results/feature_selection/volcano_plot.png", p_volcano, width = 8, height = 6, dpi = 300)
  print(p_volcano)

}, error = function(e) {
  cat("  WARNING: limma failed:", conditionMessage(e), "\n")
  fs_results$limma <<- character(0)
})

# --- 1.2 Wilcoxon Rank-Sum Test -----------------------------------------------
cat("\n--- 1.2 Wilcoxon Rank-Sum Test ---\n")
tryCatch({
  cancer_idx <- disc_labels == levels(disc_labels)[1]
  normal_idx <- !cancer_idx

  wilcox_pvals <- sapply(colnames(disc_features), function(mirna) {
    wilcox.test(disc_features[cancer_idx, mirna],
                disc_features[normal_idx, mirna])$p.value
  })

  wilcox_padj <- p.adjust(wilcox_pvals, method = "BH")

  # Median fold change
  median_cancer <- apply(disc_features[cancer_idx, ], 2, median)
  median_normal <- apply(disc_features[normal_idx, ], 2, median)
  median_fc <- abs(median_cancer - median_normal)  # on log2 scale, difference = log2 FC

  fs_wilcoxon <- names(wilcox_padj)[wilcox_padj < 0.01 & median_fc > log2(1.5)]
  fs_results$wilcoxon <- fs_wilcoxon
  cat("  Wilcoxon selected:", length(fs_wilcoxon), "miRNAs\n")

}, error = function(e) {
  cat("  WARNING: Wilcoxon test failed:", conditionMessage(e), "\n")
  fs_results$wilcoxon <<- character(0)
})

# --- 1.3 Random Forest Variable Importance -----------------------------------
cat("\n--- 1.3 Random Forest Importance ---\n")
tryCatch({
  set.seed(SEED)
  rf_fs <- randomForest(x = disc_features, y = disc_labels,
                        ntree = 1000, importance = TRUE)

  rf_imp <- importance(rf_fs, type = 1)  # MeanDecreaseAccuracy
  rf_imp_sorted <- sort(rf_imp[, 1], decreasing = TRUE)

  # Select: importance > mean + 1*SD, or top 50
  threshold_imp <- mean(rf_imp[, 1]) + sd(rf_imp[, 1])
  fs_rf_imp <- names(rf_imp_sorted[rf_imp_sorted > threshold_imp])
  if (length(fs_rf_imp) > 50) fs_rf_imp <- names(rf_imp_sorted[1:50])

  fs_results$rf_importance <- fs_rf_imp
  cat("  RF importance selected:", length(fs_rf_imp), "miRNAs\n")

  # Plot top 30
  top30 <- data.frame(
    miRNA = names(rf_imp_sorted[1:min(30, length(rf_imp_sorted))]),
    Importance = rf_imp_sorted[1:min(30, length(rf_imp_sorted))]
  )
  top30$miRNA <- factor(top30$miRNA, levels = rev(top30$miRNA))

  p_rf_imp <- ggplot(top30, aes(x = miRNA, y = Importance)) +
    geom_col(fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    labs(title = "Top 30 miRNAs by RF Importance",
         x = "", y = "Mean Decrease Accuracy")
  ggsave("results/feature_selection/rf_importance.png", p_rf_imp, width = 8, height = 8, dpi = 300)
  print(p_rf_imp)

}, error = function(e) {
  cat("  WARNING: RF importance failed:", conditionMessage(e), "\n")
  fs_results$rf_importance <<- character(0)
})

# --- 1.4 LASSO Regression (L1) -----------------------------------------------
cat("\n--- 1.4 LASSO (L1 Regularization) ---\n")
tryCatch({
  set.seed(SEED)
  x_mat <- as.matrix(disc_features)
  y_bin <- as.numeric(disc_labels) - 1  # 0/1

  cv_lasso <- cv.glmnet(x_mat, y_bin, family = "binomial", alpha = 1,
                         nfolds = 10, type.measure = "auc")

  lasso_coef <- coef(cv_lasso, s = "lambda.1se")
  fs_lasso <- rownames(lasso_coef)[which(lasso_coef[, 1] != 0)]
  fs_lasso <- setdiff(fs_lasso, "(Intercept)")

  fs_results$lasso <- fs_lasso
  cat("  LASSO selected:", length(fs_lasso), "miRNAs (at lambda.1se)\n")

  # Coefficient path plot
  png("results/feature_selection/lasso_cv_plot.png", width = 8, height = 6,
      units = "in", res = 300)
  plot(cv_lasso, main = "LASSO Cross-Validation")
  dev.off()
  plot(cv_lasso, main = "LASSO Cross-Validation")

}, error = function(e) {
  cat("  WARNING: LASSO failed:", conditionMessage(e), "\n")
  fs_results$lasso <<- character(0)
})

# --- 1.5 Elastic Net (L1+L2) -------------------------------------------------
cat("\n--- 1.5 Elastic Net ---\n")
tryCatch({
  set.seed(SEED)
  x_mat <- as.matrix(disc_features)
  y_bin <- as.numeric(disc_labels) - 1

  alpha_grid <- seq(0.1, 0.9, by = 0.1)
  best_auc <- 0
  best_alpha <- 0.5
  best_enet_model <- NULL

  for (a in alpha_grid) {
    cv_enet <- cv.glmnet(x_mat, y_bin, family = "binomial", alpha = a,
                          nfolds = 10, type.measure = "auc")
    max_auc <- max(cv_enet$cvm)
    if (max_auc > best_auc) {
      best_auc <- max_auc
      best_alpha <- a
      best_enet_model <- cv_enet
    }
  }

  cat("  Best alpha:", best_alpha, "with CV AUC:", round(best_auc, 4), "\n")

  enet_coef <- coef(best_enet_model, s = "lambda.1se")
  fs_elasticnet <- rownames(enet_coef)[which(enet_coef[, 1] != 0)]
  fs_elasticnet <- setdiff(fs_elasticnet, "(Intercept)")

  fs_results$elasticnet <- fs_elasticnet
  cat("  Elastic Net selected:", length(fs_elasticnet), "miRNAs\n")

}, error = function(e) {
  cat("  WARNING: Elastic Net failed:", conditionMessage(e), "\n")
  fs_results$elasticnet <<- character(0)
})

# --- 1.6 Recursive Feature Elimination (RFE) ---------------------------------
cat("\n--- 1.6 Recursive Feature Elimination ---\n")
tryCatch({
  set.seed(SEED)

  rfe_ctrl <- rfeControl(
    functions = rfFuncs,
    method = "repeatedcv",
    number = 10,
    repeats = 3,
    verbose = FALSE
  )

  rfe_sizes <- c(5, 10, 15, 20, 30, 50, 75, 100)
  rfe_sizes <- rfe_sizes[rfe_sizes <= ncol(disc_features)]

  rfe_result <- rfe(
    x = disc_features,
    y = disc_labels,
    sizes = rfe_sizes,
    rfeControl = rfe_ctrl
  )

  fs_rfe <- predictors(rfe_result)
  fs_results$rfe <- fs_rfe
  cat("  RFE selected:", length(fs_rfe), "miRNAs (optimal subset)\n")

  # Plot RFE performance
  rfe_perf <- data.frame(
    Variables = rfe_result$results$Variables,
    Accuracy = rfe_result$results$Accuracy
  )
  p_rfe <- ggplot(rfe_perf, aes(x = Variables, y = Accuracy)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(color = "steelblue", size = 3) +
    labs(title = "RFE: Features vs CV Accuracy",
         x = "Number of Features", y = "CV Accuracy")
  ggsave("results/feature_selection/rfe_performance.png", p_rfe, width = 8, height = 6, dpi = 300)
  print(p_rfe)

}, error = function(e) {
  cat("  WARNING: RFE failed:", conditionMessage(e), "\n")
  fs_results$rfe <<- character(0)
})

# --- 1.7 Boruta Feature Selection ---------------------------------------------
cat("\n--- 1.7 Boruta Feature Selection ---\n")
tryCatch({
  set.seed(SEED)
  boruta_result <- Boruta(x = disc_features, y = disc_labels,
                          maxRuns = 300, doTrace = 0)

  fs_boruta <- getSelectedAttributes(boruta_result, withTentative = FALSE)
  fs_results$boruta <- fs_boruta
  cat("  Boruta selected:", length(fs_boruta), "confirmed miRNAs\n")

  # Boruta plot
  png("results/feature_selection/boruta_plot.png", width = 10, height = 6,
      units = "in", res = 300)
  plot(boruta_result, xlab = "", las = 2, cex.axis = 0.6,
       main = "Boruta Feature Selection")
  dev.off()
  plot(boruta_result, xlab = "", las = 2, cex.axis = 0.6,
       main = "Boruta Feature Selection")

}, error = function(e) {
  cat("  WARNING: Boruta failed:", conditionMessage(e), "\n")
  fs_results$boruta <<- character(0)
})

# --- 1.8 mRMR (Minimum Redundancy Maximum Relevance) -------------------------
cat("\n--- 1.8 mRMR (praznik) ---\n")
tryCatch({
  # praznik::MRMR expects a data frame with the target as a factor column
  mrmr_result <- MRMR(disc_features, disc_labels, k = min(50, ncol(disc_features)))

  fs_mrmr <- colnames(disc_features)[mrmr_result$selection]
  fs_results$mrmr <- fs_mrmr
  cat("  mRMR selected:", length(fs_mrmr), "miRNAs\n")

}, error = function(e) {
  cat("  WARNING: mRMR failed:", conditionMessage(e), "\n")
  fs_results$mrmr <<- character(0)
})

# --- 1.9 Consensus Feature Set -----------------------------------------------
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  CONSENSUS FEATURE SET CONSTRUCTION\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# Remove any methods that failed (empty results)
active_methods <- names(fs_results)[sapply(fs_results, length) > 0]
cat("Active methods:", length(active_methods), "/", length(fs_results), "\n")
cat("Methods:", paste(active_methods, collapse = ", "), "\n\n")

# Frequency table: how many methods selected each miRNA
all_selected <- unique(unlist(fs_results[active_methods]))
freq_table <- data.frame(
  miRNA = all_selected,
  count = sapply(all_selected, function(m) {
    sum(sapply(fs_results[active_methods], function(fs) m %in% fs))
  }),
  stringsAsFactors = FALSE
)
freq_table$methods <- sapply(freq_table$miRNA, function(m) {
  paste(names(fs_results[active_methods])[sapply(fs_results[active_methods], function(fs) m %in% fs)],
        collapse = ", ")
})
freq_table <- freq_table[order(-freq_table$count), ]
rownames(freq_table) <- NULL

cat("Top miRNAs by selection frequency:\n")
print(head(freq_table, 20))

# UpSet plot
upset_list <- fs_results[active_methods]
# Convert to binary matrix for UpSetR
all_mirnas <- unique(unlist(upset_list))
upset_matrix <- data.frame(
  miRNA = all_mirnas,
  stringsAsFactors = FALSE
)
for (method in names(upset_list)) {
  upset_matrix[[method]] <- as.integer(all_mirnas %in% upset_list[[method]])
}

png("results/feature_selection/upset_plot.png", width = 12, height = 8,
    units = "in", res = 300)
upset(upset_matrix[, -1], nsets = length(active_methods),
      order.by = "freq", main.bar.color = "steelblue",
      sets.bar.color = COLOR_PALETTE[seq_along(active_methods)],
      text.scale = 1.3)
dev.off()

upset(upset_matrix[, -1], nsets = length(active_methods),
      order.by = "freq", main.bar.color = "steelblue",
      sets.bar.color = COLOR_PALETTE[seq_along(active_methods)],
      text.scale = 1.3)

# STRATEGY A: Strict Intersection
n_methods <- length(active_methods)
threshold_A <- ceiling(n_methods * 0.75)

consensus_A <- freq_table$miRNA[freq_table$count >= threshold_A]
cat("\nSTRATEGY A (Strict, >=", threshold_A, "/", n_methods, "methods):",
    length(consensus_A), "miRNAs\n")

# Adjust threshold if too few or too many
if (length(consensus_A) < 5 && threshold_A > 2) {
  threshold_A <- threshold_A - 1
  consensus_A <- freq_table$miRNA[freq_table$count >= threshold_A]
  cat("  Relaxed to >=", threshold_A, "methods:", length(consensus_A), "miRNAs\n")
}
if (length(consensus_A) > 50) {
  threshold_A <- threshold_A + 1
  consensus_A <- freq_table$miRNA[freq_table$count >= threshold_A]
  cat("  Tightened to >=", threshold_A, "methods:", length(consensus_A), "miRNAs\n")
}

# STRATEGY B: Minimum Optimal Set (forward stepwise by frequency)
cat("\nSTRATEGY B (Minimum Optimal Set via forward stepwise)...\n")
ranked_mirnas <- freq_table$miRNA  # already sorted by count descending
best_auc_B <- 0
no_improve_count <- 0
selected_B <- c()
auc_trajectory <- c()

set.seed(SEED)
for (i in seq_along(ranked_mirnas)) {
  candidate <- c(selected_B, ranked_mirnas[i])

  if (length(candidate) == 1) {
    # With 1 feature, simple CV
    tmp_df <- data.frame(feat = disc_features[, candidate], label = disc_labels)
    ctrl_tmp <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                             summaryFunction = twoClassSummary, verboseIter = FALSE)
    set.seed(SEED)
    tmp_model <- tryCatch({
      train(label ~ ., data = tmp_df, method = "rf", trControl = ctrl_tmp,
            metric = "ROC", tuneGrid = data.frame(mtry = 1), ntree = 500)
    }, error = function(e) NULL)

    if (!is.null(tmp_model)) {
      current_auc <- max(tmp_model$results$ROC)
    } else {
      current_auc <- 0.5
    }
  } else {
    tmp_df <- data.frame(disc_features[, candidate], label = disc_labels)
    ctrl_tmp <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                             summaryFunction = twoClassSummary, verboseIter = FALSE)
    set.seed(SEED)
    mtry_val <- min(floor(sqrt(length(candidate))), length(candidate))
    tmp_model <- tryCatch({
      train(label ~ ., data = tmp_df, method = "rf", trControl = ctrl_tmp,
            metric = "ROC", tuneGrid = data.frame(mtry = mtry_val), ntree = 500)
    }, error = function(e) NULL)

    if (!is.null(tmp_model)) {
      current_auc <- max(tmp_model$results$ROC)
    } else {
      current_auc <- best_auc_B
    }
  }

  auc_trajectory <- c(auc_trajectory, current_auc)

  if (current_auc > best_auc_B + 0.005) {
    best_auc_B <- current_auc
    selected_B <- candidate
    no_improve_count <- 0
  } else {
    no_improve_count <- no_improve_count + 1
  }

  if (no_improve_count >= 3) {
    cat("  Stopped at", length(selected_B), "miRNAs (3 consecutive non-improvements)\n")
    break
  }

  # Safety: stop at 60 features max for forward selection
  if (i >= 60) break
}

consensus_B <- selected_B
cat("STRATEGY B selected:", length(consensus_B), "miRNAs with AUC:",
    round(best_auc_B, 4), "\n")

# Plot AUC trajectory
if (length(auc_trajectory) > 1) {
  traj_df <- data.frame(n_features = seq_along(auc_trajectory), AUC = auc_trajectory)
  p_traj <- ggplot(traj_df, aes(x = n_features, y = AUC)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(color = "steelblue", size = 2) +
    geom_vline(xintercept = length(consensus_B), linetype = "dashed", color = "red") +
    labs(title = "Forward Stepwise Feature Selection",
         subtitle = paste("Optimal:", length(consensus_B), "miRNAs"),
         x = "Number of Features", y = "5-Fold CV AUC")
  ggsave("results/feature_selection/stepwise_auc.png", p_traj, width = 8, height = 6, dpi = 300)
  print(p_traj)
}

# Choose the final set: prefer smaller set if AUC difference < 0.01
# Evaluate Strategy A AUC for comparison
if (length(consensus_A) > 0) {
  set.seed(SEED)
  tmp_df_A <- data.frame(disc_features[, consensus_A, drop = FALSE], label = disc_labels)
  ctrl_tmp <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                           summaryFunction = twoClassSummary, verboseIter = FALSE)
  mtry_A <- min(floor(sqrt(length(consensus_A))), length(consensus_A))
  model_A <- tryCatch({
    train(label ~ ., data = tmp_df_A, method = "rf", trControl = ctrl_tmp,
          metric = "ROC", tuneGrid = data.frame(mtry = mtry_A), ntree = 500)
  }, error = function(e) NULL)
  auc_A <- if (!is.null(model_A)) max(model_A$results$ROC) else 0
} else {
  auc_A <- 0
}

cat("\nComparison:\n")
cat("  Strategy A:", length(consensus_A), "miRNAs, AUC:", round(auc_A, 4), "\n")
cat("  Strategy B:", length(consensus_B), "miRNAs, AUC:", round(best_auc_B, 4), "\n")

# Decision logic
if (length(consensus_A) == 0 && length(consensus_B) == 0) {
  # Fallback: LASSO + top 20 RF importance
  cat("WARNING: Both strategies yielded empty sets. Using fallback.\n")
  fallback <- unique(c(fs_results$lasso, head(freq_table$miRNA, 20)))
  final_mirnas <- fallback
} else if (length(consensus_A) == 0) {
  final_mirnas <- consensus_B
} else if (length(consensus_B) == 0) {
  final_mirnas <- consensus_A
} else if (abs(auc_A - best_auc_B) < 0.01) {
  # AUC similar — prefer smaller set
  if (length(consensus_B) <= length(consensus_A)) {
    final_mirnas <- consensus_B
    cat("  Selected: Strategy B (smaller set, similar AUC)\n")
  } else {
    final_mirnas <- consensus_A
    cat("  Selected: Strategy A (smaller set, similar AUC)\n")
  }
} else if (best_auc_B > auc_A) {
  final_mirnas <- consensus_B
  cat("  Selected: Strategy B (higher AUC)\n")
} else {
  final_mirnas <- consensus_A
  cat("  Selected: Strategy A (higher AUC)\n")
}

cat("\n*** FINAL miRNA SIGNATURE:", length(final_mirnas), "miRNAs ***\n")
cat(paste(final_mirnas, collapse = ", "), "\n")

# Individual AUC for each selected miRNA
individual_auc <- sapply(final_mirnas, function(m) {
  roc_obj <- roc(disc_labels, disc_features[, m], quiet = TRUE)
  as.numeric(auc(roc_obj))
})

# Summary table
fs_summary <- data.frame(
  miRNA = final_mirnas,
  selection_count = freq_table$count[match(final_mirnas, freq_table$miRNA)],
  methods = freq_table$methods[match(final_mirnas, freq_table$miRNA)],
  individual_AUC = round(individual_auc, 4),
  stringsAsFactors = FALSE
)
fs_summary <- fs_summary[order(-fs_summary$selection_count, -fs_summary$individual_AUC), ]
write.csv(fs_summary, "results/feature_selection/fs_summary_table.csv", row.names = FALSE)
writeLines(final_mirnas, "results/feature_selection/consensus_mirnas.txt")

cat("\nFeature selection summary saved.\n")
print(fs_summary)

# Heatmap of consensus miRNAs
ann_col <- data.frame(Label = disc_labels)
rownames(ann_col) <- rownames(disc_features)
ann_colors <- list(Label = setNames(COLOR_PALETTE[1:2], levels(disc_labels)))

pheatmap(t(disc_features[, final_mirnas]),
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         clustering_method = "ward.D2",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = paste("Consensus miRNA Signature (n =", length(final_mirnas), ")"),
         filename = "results/feature_selection/heatmap_consensus.png",
         width = 10, height = 8)

# Also display on screen
pheatmap(t(disc_features[, final_mirnas]),
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         show_colnames = FALSE,
         clustering_method = "ward.D2",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = paste("Consensus miRNA Signature (n =", length(final_mirnas), ")"))


# ==============================================================================
# PHASE 2: BUILD THE BEST DIAGNOSIS MODEL
# ==============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  PHASE 2: MODEL BUILDING & HYPERPARAMETER TUNING\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# --- 2.0 Prepare discovery data with selected features -----------------------
train_x <- disc_features[, final_mirnas, drop = FALSE]
train_y <- disc_labels

cat("Training data:", nrow(train_x), "samples x", ncol(train_x), "features\n")

# --- 2.1 Define trainControl --------------------------------------------------
set.seed(SEED)
if (IMBALANCED) {
  cat("Using SMOTE sampling to handle class imbalance.\n")
  train_ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    sampling = "smote",
    verboseIter = FALSE
  )
} else {
  train_ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    repeats = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final",
    verboseIter = FALSE
  )
}

# Storage
all_models <- list()

# --- MODEL 1: Random Forest --------------------------------------------------
cat("\n--- Training Model 1: Random Forest ---\n")
set.seed(SEED)
mtry_max <- min(ncol(train_x), 20)
rf_grid <- expand.grid(mtry = seq(2, mtry_max, by = max(1, floor(mtry_max / 10))))

all_models$RF <- train(
  x = train_x, y = train_y,
  method = "rf",
  trControl = train_ctrl,
  tuneGrid = rf_grid,
  ntree = 1000,
  importance = TRUE,
  metric = "ROC"
)
cat("  Best mtry:", all_models$RF$bestTune$mtry,
    "| CV AUC:", round(max(all_models$RF$results$ROC), 4), "\n")

# --- MODEL 2: SVM (Radial) ---------------------------------------------------
cat("\n--- Training Model 2: SVM (Radial) ---\n")
set.seed(SEED)
svm_grid <- expand.grid(
  C = c(0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100),
  sigma = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5, 1)
)

all_models$SVM <- train(
  x = train_x, y = train_y,
  method = "svmRadial",
  trControl = train_ctrl,
  tuneGrid = svm_grid,
  preProcess = c("center", "scale"),
  metric = "ROC"
)
cat("  Best C:", all_models$SVM$bestTune$C,
    "sigma:", all_models$SVM$bestTune$sigma,
    "| CV AUC:", round(max(all_models$SVM$results$ROC), 4), "\n")

# --- MODEL 3: XGBoost --------------------------------------------------------
cat("\n--- Training Model 3: XGBoost ---\n")
set.seed(SEED)

# Use random search for XGBoost (full grid is too large)
train_ctrl_xgb <- train_ctrl
train_ctrl_xgb$search <- "random"

all_models$XGB <- train(
  x = train_x, y = train_y,
  method = "xgbTree",
  trControl = train_ctrl_xgb,
  tuneLength = 200,
  metric = "ROC",
  verbosity = 0
)
cat("  Best params — nrounds:", all_models$XGB$bestTune$nrounds,
    "depth:", all_models$XGB$bestTune$max_depth,
    "eta:", all_models$XGB$bestTune$eta,
    "| CV AUC:", round(max(all_models$XGB$results$ROC), 4), "\n")

# --- MODEL 4: Elastic Net Logistic Regression ---------------------------------
cat("\n--- Training Model 4: Elastic Net (glmnet) ---\n")
set.seed(SEED)
glm_grid <- expand.grid(
  alpha = seq(0, 1, by = 0.1),
  lambda = 10^seq(-4, 0, length.out = 50)
)

all_models$GLM <- train(
  x = train_x, y = train_y,
  method = "glmnet",
  trControl = train_ctrl,
  tuneGrid = glm_grid,
  preProcess = c("center", "scale"),
  metric = "ROC"
)
cat("  Best alpha:", all_models$GLM$bestTune$alpha,
    "lambda:", round(all_models$GLM$bestTune$lambda, 6),
    "| CV AUC:", round(max(all_models$GLM$results$ROC), 4), "\n")

# --- MODEL 5: k-Nearest Neighbors --------------------------------------------
cat("\n--- Training Model 5: KNN ---\n")
set.seed(SEED)
knn_grid <- expand.grid(k = seq(1, 41, by = 2))

all_models$KNN <- train(
  x = train_x, y = train_y,
  method = "knn",
  trControl = train_ctrl,
  tuneGrid = knn_grid,
  preProcess = c("center", "scale"),
  metric = "ROC"
)
cat("  Best k:", all_models$KNN$bestTune$k,
    "| CV AUC:", round(max(all_models$KNN$results$ROC), 4), "\n")

# --- MODEL 6: Neural Network -------------------------------------------------
cat("\n--- Training Model 6: Neural Network ---\n")
set.seed(SEED)
nnet_grid <- expand.grid(
  size = c(1, 3, 5, 7, 10),
  decay = c(0, 0.001, 0.01, 0.1, 1)
)

all_models$NNET <- train(
  x = train_x, y = train_y,
  method = "nnet",
  trControl = train_ctrl,
  tuneGrid = nnet_grid,
  preProcess = c("center", "scale"),
  MaxNWts = 5000,
  maxit = 500,
  trace = FALSE,
  metric = "ROC"
)
cat("  Best size:", all_models$NNET$bestTune$size,
    "decay:", all_models$NNET$bestTune$decay,
    "| CV AUC:", round(max(all_models$NNET$results$ROC), 4), "\n")

# --- 2.2 Model Comparison ----------------------------------------------------
cat("\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
cat("  MODEL COMPARISON\n")
cat(paste(rep("-", 70), collapse = ""), "\n\n")

cv_results <- resamples(all_models)
cv_summary <- summary(cv_results)

# Performance summary table
model_perf <- data.frame(
  Model = names(all_models),
  AUC_mean = cv_summary$statistics$ROC[, "Mean"],
  AUC_sd = cv_summary$statistics$ROC[, "SD"],
  Sens_mean = cv_summary$statistics$Sens[, "Mean"],
  Spec_mean = cv_summary$statistics$Spec[, "Mean"],
  stringsAsFactors = FALSE
)
model_perf <- model_perf[order(-model_perf$AUC_mean), ]
cat("Model comparison (sorted by mean AUC):\n")
print(model_perf)

write.csv(model_perf, "results/models/model_comparison.csv", row.names = FALSE)

# Statistical comparison
cv_diffs <- diff(cv_results)

# Boxplot of CV AUC
cv_long <- data.frame()
for (m in names(all_models)) {
  vals <- cv_results$values[[paste0(m, "~ROC")]]
  cv_long <- rbind(cv_long, data.frame(Model = m, AUC = vals))
}
cv_long$Model <- factor(cv_long$Model, levels = model_perf$Model)

p_cv_box <- ggplot(cv_long, aes(x = Model, y = AUC, fill = Model)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
  scale_fill_manual(values = COLOR_PALETTE) +
  labs(title = "Cross-Validation AUC Distribution",
       subtitle = "10-fold CV, repeated 5 times",
       x = "", y = "AUC") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1))
ggsave("results/models/cv_boxplot.png", p_cv_box, width = 8, height = 6, dpi = 300)
print(p_cv_box)

# Select the best model
best_model_name <- model_perf$Model[1]
cat("\n*** BEST MODEL:", best_model_name, "***\n")
cat("  Mean AUC:", round(model_perf$AUC_mean[1], 4),
    "(SD:", round(model_perf$AUC_sd[1], 4), ")\n")

best_model <- all_models[[best_model_name]]

# --- 2.3 Model Refinement (best model fine-tuning) ---------------------------
cat("\n--- 2.3 Refining Best Model ---\n")

# Fine-tune: narrow grid around best hyperparameters
# This is model-specific; we handle the top models
if (best_model_name == "RF") {
  best_mtry <- best_model$bestTune$mtry
  fine_grid <- expand.grid(mtry = seq(max(1, best_mtry - 3), best_mtry + 3, by = 1))
  set.seed(SEED)
  best_model_refined <- train(
    x = train_x, y = train_y, method = "rf",
    trControl = train_ctrl, tuneGrid = fine_grid,
    ntree = 1500, importance = TRUE, metric = "ROC"
  )
} else if (best_model_name == "SVM") {
  best_C <- best_model$bestTune$C
  best_sig <- best_model$bestTune$sigma
  fine_grid <- expand.grid(
    C = best_C * c(0.5, 0.75, 1, 1.5, 2),
    sigma = best_sig * c(0.5, 0.75, 1, 1.5, 2)
  )
  set.seed(SEED)
  best_model_refined <- train(
    x = train_x, y = train_y, method = "svmRadial",
    trControl = train_ctrl, tuneGrid = fine_grid,
    preProcess = c("center", "scale"), metric = "ROC"
  )
} else if (best_model_name == "XGB") {
  bt <- best_model$bestTune
  fine_grid <- expand.grid(
    nrounds = bt$nrounds * c(0.75, 1, 1.25),
    max_depth = c(max(1, bt$max_depth - 1), bt$max_depth, bt$max_depth + 1),
    eta = bt$eta * c(0.5, 1, 2),
    gamma = bt$gamma,
    colsample_bytree = bt$colsample_bytree,
    min_child_weight = bt$min_child_weight,
    subsample = bt$subsample
  )
  fine_grid$nrounds <- as.integer(round(fine_grid$nrounds))
  set.seed(SEED)
  best_model_refined <- train(
    x = train_x, y = train_y, method = "xgbTree",
    trControl = train_ctrl, tuneGrid = fine_grid,
    metric = "ROC", verbosity = 0
  )
} else if (best_model_name == "GLM") {
  best_a <- best_model$bestTune$alpha
  best_l <- best_model$bestTune$lambda
  fine_grid <- expand.grid(
    alpha = seq(max(0, best_a - 0.1), min(1, best_a + 0.1), by = 0.02),
    lambda = best_l * 10^seq(-1, 1, length.out = 20)
  )
  set.seed(SEED)
  best_model_refined <- train(
    x = train_x, y = train_y, method = "glmnet",
    trControl = train_ctrl, tuneGrid = fine_grid,
    preProcess = c("center", "scale"), metric = "ROC"
  )
} else if (best_model_name == "KNN") {
  best_k <- best_model$bestTune$k
  fine_grid <- expand.grid(k = seq(max(1, best_k - 6), best_k + 6, by = 1))
  set.seed(SEED)
  best_model_refined <- train(
    x = train_x, y = train_y, method = "knn",
    trControl = train_ctrl, tuneGrid = fine_grid,
    preProcess = c("center", "scale"), metric = "ROC"
  )
} else if (best_model_name == "NNET") {
  bt <- best_model$bestTune
  fine_grid <- expand.grid(
    size = c(max(1, bt$size - 1), bt$size, bt$size + 1, bt$size + 2),
    decay = bt$decay * c(0.1, 0.5, 1, 2, 10)
  )
  fine_grid <- fine_grid[fine_grid$decay >= 0, ]
  set.seed(SEED)
  best_model_refined <- train(
    x = train_x, y = train_y, method = "nnet",
    trControl = train_ctrl, tuneGrid = fine_grid,
    preProcess = c("center", "scale"),
    MaxNWts = 5000, maxit = 500, trace = FALSE, metric = "ROC"
  )
} else {
  best_model_refined <- best_model
}

refined_auc <- max(best_model_refined$results$ROC)
original_auc <- max(best_model$results$ROC)
cat("  Original AUC:", round(original_auc, 4),
    "→ Refined AUC:", round(refined_auc, 4), "\n")

if (refined_auc >= original_auc) {
  final_model <- best_model_refined
  cat("  Using refined model.\n")
} else {
  final_model <- best_model
  cat("  Refinement did not improve; keeping original.\n")
}

# Optimal probability threshold (Youden's J)
# Use CV predictions from the final model
cv_preds <- final_model$pred
# For repeated CV, average predictions per sample
cv_preds_agg <- cv_preds %>%
  group_by(rowIndex) %>%
  summarise(
    obs = first(obs),
    prob = mean(.data[[levels(train_y)[2]]])  # probability of second class
  )

roc_train <- roc(cv_preds_agg$obs, cv_preds_agg$prob, quiet = TRUE)

# Youden's J
youden_idx <- which.max(roc_train$sensitivities + roc_train$specificities - 1)
threshold_youden <- roc_train$thresholds[youden_idx]

# F1 threshold search
thresholds_search <- seq(0.1, 0.9, by = 0.01)
f1_scores <- sapply(thresholds_search, function(thr) {
  pred_class <- ifelse(cv_preds_agg$prob >= thr, levels(train_y)[2], levels(train_y)[1])
  pred_class <- factor(pred_class, levels = levels(train_y))
  cm <- confusionMatrix(pred_class, cv_preds_agg$obs, positive = levels(train_y)[2])
  f1 <- 2 * cm$byClass["Precision"] * cm$byClass["Recall"] /
    (cm$byClass["Precision"] + cm$byClass["Recall"])
  ifelse(is.na(f1), 0, f1)
})
threshold_f1 <- thresholds_search[which.max(f1_scores)]

cat("\nOptimal thresholds:\n")
cat("  Youden's J:", round(threshold_youden, 4), "\n")
cat("  Max F1:", round(threshold_f1, 4), "\n")
cat("  Recommended (diagnostic): Youden's J =", round(threshold_youden, 4), "\n")

optimal_threshold <- threshold_youden

# Calibration curve
cal_df <- data.frame(
  predicted = cv_preds_agg$prob,
  observed = as.numeric(cv_preds_agg$obs == levels(train_y)[2])
)
cal_df$bin <- cut(cal_df$predicted, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
cal_summary <- cal_df %>%
  group_by(bin) %>%
  summarise(
    mean_predicted = mean(predicted),
    mean_observed = mean(observed),
    n = n()
  ) %>%
  filter(!is.na(bin))

p_cal <- ggplot(cal_summary, aes(x = mean_predicted, y = mean_observed)) +
  geom_point(aes(size = n), color = "steelblue") +
  geom_line(color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  labs(title = "Calibration Curve (CV Predictions)",
       x = "Mean Predicted Probability", y = "Observed Proportion",
       size = "N samples") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1))
ggsave("results/models/calibration_curve.png", p_cal, width = 7, height = 7, dpi = 300)
print(p_cal)

# Save the final model and preprocessing params
saveRDS(final_model, "results/models/best_model.rds")
saveRDS(preproc_params, "results/models/preprocessing_params.rds")
saveRDS(final_mirnas, "results/models/selected_mirnas.rds")

cat("\nPhase 2 complete. Best model saved.\n")

# ROC curves for all models on CV predictions
p_roc_all <- ggplot() + theme_minimal(base_size = 13)
roc_plot_data <- data.frame()
auc_labels <- c()

for (i in seq_along(all_models)) {
  m_name <- names(all_models)[i]
  m_preds <- all_models[[m_name]]$pred
  m_agg <- m_preds %>%
    group_by(rowIndex) %>%
    summarise(obs = first(obs),
              prob = mean(.data[[levels(train_y)[2]]]))
  m_roc <- roc(m_agg$obs, m_agg$prob, quiet = TRUE)
  m_auc <- round(as.numeric(auc(m_roc)), 3)
  auc_labels <- c(auc_labels, paste0(m_name, " (AUC=", m_auc, ")"))

  roc_plot_data <- rbind(roc_plot_data, data.frame(
    FPR = 1 - m_roc$specificities,
    TPR = m_roc$sensitivities,
    Model = paste0(m_name, " (AUC=", m_auc, ")")
  ))
}

p_roc_all <- ggplot(roc_plot_data, aes(x = FPR, y = TPR, color = Model)) +
  geom_line(size = 1, alpha = 0.8) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = COLOR_PALETTE[seq_along(all_models)]) +
  labs(title = "ROC Curves — All Models (CV)",
       x = "False Positive Rate", y = "True Positive Rate",
       color = "Model") +
  coord_equal() +
  theme(legend.position = "right")
ggsave("results/models/roc_all_models.png", p_roc_all, width = 9, height = 7, dpi = 300)
print(p_roc_all)


# ==============================================================================
# PHASE 3: INDEPENDENT VALIDATION
# ==============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  PHASE 3: INDEPENDENT VALIDATION\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# --- 3.0 Prepare validation set ----------------------------------------------
cat("--- 3.0 Preparing validation data ---\n")

val_x <- val_features[, final_mirnas, drop = FALSE]
val_y <- val_labels

# Apply SAME preprocessing from training (center + scale with training parameters)
val_x_scaled <- predict(preproc_params, val_features[, final_mirnas, drop = FALSE])

cat("Validation set:", nrow(val_x), "samples x", ncol(val_x), "features\n")
cat("Validation class distribution:\n")
print(table(val_y))

# --- 3.1 Predict & Evaluate --------------------------------------------------
cat("\n--- 3.1 Predictions on Validation Set ---\n")

val_prob <- predict(final_model, val_x, type = "prob")
prob_positive <- val_prob[, levels(train_y)[2]]

# Apply optimal threshold
val_pred_opt <- ifelse(prob_positive >= optimal_threshold,
                       levels(train_y)[2], levels(train_y)[1])
val_pred_opt <- factor(val_pred_opt, levels = levels(train_y))

# Also default 0.5 threshold
val_pred_default <- predict(final_model, val_x)

# Confusion matrix (optimal threshold)
cm_val <- confusionMatrix(val_pred_opt, val_y, positive = levels(train_y)[2])
cat("\nValidation Results (threshold =", round(optimal_threshold, 4), "):\n")
print(cm_val)

# AUC with 95% CI (DeLong method)
roc_val <- roc(val_y, prob_positive, quiet = TRUE)
auc_val <- as.numeric(auc(roc_val))
ci_val <- ci.auc(roc_val, method = "delong")

cat("\nValidation AUC:", round(auc_val, 4),
    "(95% CI:", round(ci_val[1], 4), "-", round(ci_val[3], 4), ")\n")

# Full metrics
accuracy_val <- cm_val$overall["Accuracy"]
sens_val <- cm_val$byClass["Sensitivity"]
spec_val <- cm_val$byClass["Specificity"]
ppv_val <- cm_val$byClass["Pos Pred Value"]
npv_val <- cm_val$byClass["Neg Pred Value"]
f1_val <- 2 * ppv_val * sens_val / (ppv_val + sens_val)

val_metrics <- data.frame(
  Metric = c("AUC", "Accuracy", "Sensitivity", "Specificity", "PPV", "NPV", "F1"),
  Value = round(c(auc_val, accuracy_val, sens_val, spec_val, ppv_val, npv_val, f1_val), 4)
)
cat("\nValidation metrics summary:\n")
print(val_metrics)

write.csv(val_metrics, "results/validation/validation_results.csv", row.names = FALSE)

# Compare training CV vs validation performance
cv_auc_train <- round(max(final_model$results$ROC), 4)
cat("\nPerformance comparison:\n")
cat("  Training CV AUC:", cv_auc_train, "\n")
cat("  Validation AUC:", round(auc_val, 4), "\n")
if (cv_auc_train - auc_val > 0.05) {
  cat("  WARNING: AUC dropped >0.05 — potential overfitting!\n")
} else {
  cat("  Good generalization (AUC drop <=0.05)\n")
}

# ROC curve with CI band
roc_ci <- ci.se(roc_val, specificities = seq(0, 1, 0.01), conf.level = 0.95)
roc_ci_df <- data.frame(
  Specificity = as.numeric(rownames(roc_ci)),
  Sensitivity = roc_ci[, 2],
  lower = roc_ci[, 1],
  upper = roc_ci[, 3]
)

p_roc_val <- ggplot(roc_ci_df, aes(x = 1 - Specificity, y = Sensitivity)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  labs(title = "Validation ROC Curve with 95% CI",
       subtitle = paste0("AUC = ", round(auc_val, 3),
                         " (95% CI: ", round(ci_val[1], 3), "–", round(ci_val[3], 3), ")"),
       x = "False Positive Rate", y = "True Positive Rate") +
  coord_equal()
ggsave("results/validation/roc_validation.png", p_roc_val, width = 7, height = 7, dpi = 300)
print(p_roc_val)

# Confusion matrix heatmap
cm_data <- as.data.frame(cm_val$table)
p_cm <- ggplot(cm_data, aes(x = Prediction, y = Reference, fill = Freq)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = Freq), size = 8, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Confusion Matrix (Validation Set)",
       subtitle = paste("Threshold:", round(optimal_threshold, 3)),
       x = "Predicted", y = "Actual") +
  theme(legend.position = "none")
ggsave("results/validation/confusion_matrix.png", p_cm, width = 6, height = 5, dpi = 300)
print(p_cm)

# Precision-Recall curve
pr_obj <- pr.curve(scores.class0 = prob_positive[val_y == levels(train_y)[2]],
                   scores.class1 = prob_positive[val_y == levels(train_y)[1]],
                   curve = TRUE)
pr_auc_val <- pr_obj$auc.integral

pr_df <- data.frame(Recall = pr_obj$curve[, 1], Precision = pr_obj$curve[, 2])
p_pr <- ggplot(pr_df, aes(x = Recall, y = Precision)) +
  geom_line(color = "#E31A1C", size = 1.2) +
  labs(title = "Precision-Recall Curve (Validation)",
       subtitle = paste0("PR-AUC = ", round(pr_auc_val, 3)),
       x = "Recall (Sensitivity)", y = "Precision (PPV)") +
  ylim(0, 1)
ggsave("results/validation/pr_curve.png", p_pr, width = 7, height = 6, dpi = 300)
print(p_pr)

# --- 3.2 Robustness Checks ---------------------------------------------------
cat("\n--- 3.2 Robustness Checks ---\n")

# (a) Bootstrap validation — 1000 resamples
cat("Running bootstrap validation (1000 resamples)...\n")
set.seed(SEED)
n_boot <- 1000
boot_aucs <- numeric(n_boot)

for (b in seq_len(n_boot)) {
  boot_idx <- sample(seq_len(nrow(val_x)), replace = TRUE)
  boot_y <- val_y[boot_idx]

  # Skip if only one class in bootstrap sample
  if (length(unique(boot_y)) < 2) {
    boot_aucs[b] <- NA
    next
  }

  boot_prob <- prob_positive[boot_idx]
  boot_roc <- tryCatch(roc(boot_y, boot_prob, quiet = TRUE), error = function(e) NULL)
  boot_aucs[b] <- if (!is.null(boot_roc)) as.numeric(auc(boot_roc)) else NA
}

boot_aucs <- boot_aucs[!is.na(boot_aucs)]
boot_ci <- quantile(boot_aucs, probs = c(0.025, 0.975))

cat("  Bootstrap AUC — Median:", round(median(boot_aucs), 4),
    "95% CI: [", round(boot_ci[1], 4), ",", round(boot_ci[2], 4), "]\n")

p_boot <- ggplot(data.frame(AUC = boot_aucs), aes(x = AUC)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "white") +
  geom_vline(xintercept = auc_val, color = "#E31A1C", size = 1.2, linetype = "solid") +
  geom_vline(xintercept = boot_ci, color = "grey40", linetype = "dashed") +
  labs(title = "Bootstrap Distribution of Validation AUC",
       subtitle = paste0("Median = ", round(median(boot_aucs), 3),
                         ", 95% CI: [", round(boot_ci[1], 3), ", ", round(boot_ci[2], 3), "]"),
       x = "AUC", y = "Count")
ggsave("results/validation/bootstrap_ci.png", p_boot, width = 8, height = 6, dpi = 300)
print(p_boot)

# (b) Permutation test — 1000 permutations
cat("Running permutation test (1000 permutations)...\n")
set.seed(SEED)
n_perm <- 1000
perm_aucs <- numeric(n_perm)

for (p in seq_len(n_perm)) {
  perm_y <- sample(val_y)
  perm_roc <- tryCatch(roc(perm_y, prob_positive, quiet = TRUE), error = function(e) NULL)
  perm_aucs[p] <- if (!is.null(perm_roc)) as.numeric(auc(perm_roc)) else 0.5
}

perm_pval <- mean(perm_aucs >= auc_val)
cat("  Permutation p-value:", perm_pval, "\n")

p_perm <- ggplot(data.frame(AUC = perm_aucs), aes(x = AUC)) +
  geom_histogram(bins = 50, fill = "grey70", color = "white") +
  geom_vline(xintercept = auc_val, color = "#E31A1C", size = 1.2) +
  annotate("text", x = auc_val, y = Inf, label = paste("Observed AUC =", round(auc_val, 3)),
           vjust = 2, hjust = -0.1, color = "#E31A1C", fontface = "bold") +
  labs(title = "Permutation Test: Null Distribution of AUC",
       subtitle = paste0("p-value = ", perm_pval),
       x = "AUC (permuted labels)", y = "Count")
ggsave("results/validation/permutation_test.png", p_perm, width = 8, height = 6, dpi = 300)
print(p_perm)

# (c) Per-sample prediction confidence
pred_conf_df <- data.frame(
  Probability = prob_positive,
  True_Label = val_y,
  Confidence = ifelse(prob_positive > 0.6 | prob_positive < 0.4, "Confident", "Uncertain")
)

p_pred_dist <- ggplot(pred_conf_df, aes(x = Probability, fill = True_Label)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity", color = "white") +
  geom_vline(xintercept = c(0.4, 0.6), linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = COLOR_PALETTE[1:2]) +
  labs(title = "Prediction Probability Distribution",
       subtitle = paste("Uncertain zone (0.4–0.6):",
                        sum(pred_conf_df$Confidence == "Uncertain"), "samples"),
       x = "Predicted Probability", y = "Count",
       fill = "True Label")
ggsave("results/validation/prediction_distribution.png", p_pred_dist,
       width = 8, height = 6, dpi = 300)
print(p_pred_dist)

# --- 3.3 Clinical Interpretation Summary -------------------------------------
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  FINAL CLINICAL REPORT\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

report_lines <- c(
  paste(rep("=", 70), collapse = ""),
  "  BREAST CANCER miRNA DIAGNOSTIC MODEL — FINAL REPORT",
  paste(rep("=", 70), collapse = ""),
  "",
  paste("Diagnostic miRNA Signature:", length(final_mirnas), "miRNAs"),
  paste("  ", paste(final_mirnas, collapse = ", ")),
  "",
  paste("Best Model:", best_model_name),
  paste("Best Hyperparameters:"),
  paste("  ", paste(names(final_model$bestTune),
                    final_model$bestTune, sep = " = ", collapse = ", ")),
  "",
  paste("Discovery Performance (10x5 CV):"),
  paste("  AUC =", round(cv_auc_train, 4)),
  "",
  paste("Validation Performance:"),
  paste("  AUC =", round(auc_val, 4),
        "(95% CI:", round(ci_val[1], 4), "–", round(ci_val[3], 4), ")"),
  paste("  Accuracy =", round(accuracy_val, 4)),
  paste("  Sensitivity =", round(sens_val, 4)),
  paste("  Specificity =", round(spec_val, 4)),
  paste("  PPV =", round(ppv_val, 4)),
  paste("  NPV =", round(npv_val, 4)),
  paste("  F1 =", round(f1_val, 4)),
  paste("  PR-AUC =", round(pr_auc_val, 4)),
  "",
  paste("Optimal Threshold:", round(optimal_threshold, 4), "(Youden's J)"),
  paste("  At this threshold:"),
  paste("    Sensitivity:", round(sens_val * 100, 1), "%"),
  paste("    Specificity:", round(spec_val * 100, 1), "%"),
  "",
  paste("Bootstrap 95% CI for Validation AUC:",
        "[", round(boot_ci[1], 4), ",", round(boot_ci[2], 4), "]"),
  paste("Permutation test p-value:", perm_pval),
  "",
  paste(rep("=", 70), collapse = "")
)

cat(paste(report_lines, collapse = "\n"))
cat("\n")

writeLines(report_lines, "results/validation/final_report.txt")

# ==============================================================================
# SESSION INFO & RUNTIME
# ==============================================================================

end_time <- Sys.time()
runtime <- difftime(end_time, start_time, units = "mins")
cat("\nTotal runtime:", round(as.numeric(runtime), 2), "minutes\n")
cat("\nSession info:\n")
sessionInfo()

cat("\n*** PIPELINE COMPLETE ***\n")
cat("All results saved to results/ directory.\n")

################################################################################
#                                                                              #
#   03 — Model Training: 6 Classifiers + Comparison (PHASE 2)                 #
#                                                                              #
#   Train RF (primary) + 5 comparison models on train_df using consensus      #
#   piRNAs. 10-fold CV x5 repeats. Select best, fine-tune, find optimal       #
#   threshold.                                                                 #
#                                                                              #
#   Requires: results from 01 + 02                                            #
#   Produces: results/models/final_model.rds + comparison plots               #
#                                                                              #
################################################################################

start_time <- Sys.time()
cat("=================================================================\n")
cat("  PHASE 2: Model Training — 6 Candidate Classifiers\n")
cat("=================================================================\n\n")

# ==============================================================================
# 0. PACKAGES
# ==============================================================================
required_pkgs <- c(
  "caret", "randomForest", "e1071", "glmnet", "xgboost", "nnet",
  "pROC", "PRROC", "ggplot2", "dplyr", "cowplot"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(caret)
  library(randomForest)
  library(e1071)
  library(glmnet)
  library(xgboost)
  library(nnet)
  library(pROC)
  library(PRROC)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})
cat("All packages loaded.\n\n")

SEED <- 2024
set.seed(SEED)

my_theme <- theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")

# ==============================================================================
# 0.1 LOAD PHASE 0+1 OUTPUTS
# ==============================================================================
cat("Loading previous outputs...\n")
train_df   <- readRDS("results/models/train_df.rds")
top_feats  <- readRDS("results/models/final_features.rds")
gene_cols  <- readRDS("results/models/common_pirnas.rds")

cat("  Training samples:", nrow(train_df), "\n")
cat("  Signature piRNAs:", length(top_feats), "\n")
cat("  Features:", paste(top_feats, collapse = ", "), "\n\n")

# Prepare training data
train_data <- train_df[, c("Group", top_feats)]

# Check class imbalance
class_table <- table(train_data$Group)
imbalance_ratio <- max(class_table) / min(class_table)
cat("Class distribution:", class_table["Normal"], "Normal,",
    class_table["Tumor"], "Tumor\n")
cat("Imbalance ratio:", round(imbalance_ratio, 2), "\n")

use_sampling <- if (imbalance_ratio > 2) "down" else NULL
if (!is.null(use_sampling)) cat("Using downsampling due to imbalance > 2:1\n")

# ==============================================================================
# 1. COMMON TRAIN CONTROL
# ==============================================================================
fitControl <- trainControl(
  method          = "repeatedcv",
  number          = 10,
  repeats         = 5,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  sampling        = use_sampling
)

# ==============================================================================
# 2. TRAIN 6 CANDIDATE MODELS
# ==============================================================================
cat("\n========== Training 6 Models ==========\n")
model_list <- list()

# --- 2.1 Random Forest ---
cat("\n--- Model 1/6: Random Forest ---\n")
set.seed(SEED)
mtry_vals <- seq(2, min(length(top_feats), 10), by = 1)
model_list$RF <- train(
  Group ~ ., data = train_data,
  method = "rf", metric = "ROC", trControl = fitControl,
  tuneGrid = data.frame(mtry = mtry_vals),
  ntree = 1000, importance = TRUE
)
cat("  Best mtry:", model_list$RF$bestTune$mtry,
    "  CV AUC:", round(max(model_list$RF$results$ROC), 4), "\n")

# --- 2.2 SVM-RBF ---
cat("\n--- Model 2/6: SVM-RBF ---\n")
set.seed(SEED)
svm_grid <- expand.grid(
  C     = c(0.01, 0.1, 1, 10, 100),
  sigma = c(0.001, 0.01, 0.1, 0.5)
)
model_list$SVM <- train(
  Group ~ ., data = train_data,
  method = "svmRadial", metric = "ROC", trControl = fitControl,
  tuneGrid = svm_grid,
  preProcess = c("center", "scale")
)
cat("  Best C:", model_list$SVM$bestTune$C,
    " sigma:", model_list$SVM$bestTune$sigma,
    "  CV AUC:", round(max(model_list$SVM$results$ROC), 4), "\n")

# --- 2.3 XGBoost ---
cat("\n--- Model 3/6: XGBoost ---\n")
set.seed(SEED)
xgb_ctrl <- trainControl(
  method = "repeatedcv", number = 10, repeats = 5,
  classProbs = TRUE, summaryFunction = twoClassSummary,
  savePredictions = "final", sampling = use_sampling,
  search = "random"
)
model_list$XGB <- train(
  Group ~ ., data = train_data,
  method = "xgbTree", metric = "ROC", trControl = xgb_ctrl,
  tuneLength = 50,
  verbosity = 0
)
cat("  CV AUC:", round(max(model_list$XGB$results$ROC), 4), "\n")

# --- 2.4 Elastic Net ---
cat("\n--- Model 4/6: Elastic Net ---\n")
set.seed(SEED)
en_grid <- expand.grid(
  alpha  = seq(0, 1, by = 0.1),
  lambda = 10^seq(-4, 0, length = 30)
)
model_list$ElasticNet <- train(
  Group ~ ., data = train_data,
  method = "glmnet", metric = "ROC", trControl = fitControl,
  tuneGrid = en_grid,
  preProcess = c("center", "scale")
)
cat("  Best alpha:", model_list$ElasticNet$bestTune$alpha,
    " lambda:", format(model_list$ElasticNet$bestTune$lambda, scientific = TRUE),
    "  CV AUC:", round(max(model_list$ElasticNet$results$ROC), 4), "\n")

# --- 2.5 KNN ---
cat("\n--- Model 5/6: KNN ---\n")
set.seed(SEED)
knn_grid <- data.frame(k = seq(1, 31, by = 2))
model_list$KNN <- train(
  Group ~ ., data = train_data,
  method = "knn", metric = "ROC", trControl = fitControl,
  tuneGrid = knn_grid,
  preProcess = c("center", "scale")
)
cat("  Best k:", model_list$KNN$bestTune$k,
    "  CV AUC:", round(max(model_list$KNN$results$ROC), 4), "\n")

# --- 2.6 Neural Net ---
cat("\n--- Model 6/6: Neural Net ---\n")
set.seed(SEED)
nnet_grid <- expand.grid(
  size  = c(1, 3, 5, 7),
  decay = c(0, 0.001, 0.01, 0.1)
)
model_list$NNet <- train(
  Group ~ ., data = train_data,
  method = "nnet", metric = "ROC", trControl = fitControl,
  tuneGrid = nnet_grid,
  preProcess = c("center", "scale"),
  trace = FALSE, MaxNWts = 5000
)
cat("  Best size:", model_list$NNet$bestTune$size,
    " decay:", model_list$NNet$bestTune$decay,
    "  CV AUC:", round(max(model_list$NNet$results$ROC), 4), "\n")

# ==============================================================================
# 3. MODEL COMPARISON
# ==============================================================================
cat("\n========== Model Comparison ==========\n")

# 3a. Resamples comparison
resamp <- resamples(model_list)
resamp_summary <- summary(resamp)

cat("\nCV AUC Summary:\n")
print(resamp_summary$statistics$ROC)

# Comparison table
comp_table <- data.frame(
  Model   = names(model_list),
  AUC_Median = sapply(model_list, function(m) median(m$resample$ROC)),
  AUC_Mean   = sapply(model_list, function(m) mean(m$resample$ROC)),
  Sens_Mean  = sapply(model_list, function(m) mean(m$resample$Sens)),
  Spec_Mean  = sapply(model_list, function(m) mean(m$resample$Spec)),
  stringsAsFactors = FALSE
)
comp_table <- comp_table[order(-comp_table$AUC_Median), ]
cat("\nModel Comparison Table:\n")
print(comp_table)
write.csv(comp_table, "results/models/model_comparison_table.csv", row.names = FALSE)

# 3b. Boxplot of CV AUC
resamp_df <- data.frame(
  Model = rep(names(model_list), each = nrow(resamp$values) - 0),
  AUC   = unlist(lapply(names(model_list), function(n) model_list[[n]]$resample$ROC))
)
p_box <- ggplot(resamp_df, aes(x = reorder(Model, AUC, FUN = median), y = AUC, fill = Model)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1.5) +
  scale_fill_manual(values = c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928")) +
  coord_flip() +
  labs(title = "10x5 Repeated CV AUC — Model Comparison",
       x = "", y = "AUC (ROC)") +
  my_theme + theme(legend.position = "none")
ggsave("results/models/cv_boxplot.png", p_box, width = 8, height = 6, dpi = 300, bg = "white")
cat("  Saved: results/models/cv_boxplot.png\n")

# 3c. Statistical comparison
cat("\nStatistical differences (p-values):\n")
diffs <- diff(resamp)
print(summary(diffs)$table$ROC)

# 3d. Overlaid ROC curves from CV predictions
cat("\nGenerating overlaid ROC plot...\n")
roc_plot_data <- list()
for (nm in names(model_list)) {
  pred_cv <- model_list[[nm]]$pred
  bt <- model_list[[nm]]$bestTune
  for (col_nm in names(bt)) pred_cv <- pred_cv[pred_cv[[col_nm]] == bt[[col_nm]], ]
  roc_obj <- roc(pred_cv$obs, pred_cv$Tumor, levels = c("Normal", "Tumor"), quiet = TRUE)
  roc_plot_data[[nm]] <- data.frame(
    fpr   = 1 - roc_obj$specificities,
    tpr   = roc_obj$sensitivities,
    Model = sprintf("%s (AUC=%.3f)", nm, as.numeric(auc(roc_obj)))
  )
}
roc_all_df <- do.call(rbind, roc_plot_data)
roc_all_df$Model <- factor(roc_all_df$Model, levels = unique(roc_all_df$Model))

p_roc_all <- ggplot(roc_all_df, aes(fpr, tpr, color = Model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_path(linewidth = 1) +
  scale_color_manual(values = c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928")) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(title = "ROC Curves — All 6 Candidate Models (CV)",
       x = "False Positive Rate", y = "True Positive Rate", color = "") +
  my_theme +
  theme(legend.position = c(0.7, 0.25),
        legend.background = element_rect(fill = alpha("white", 0.9), color = "grey70"))
ggsave("results/models/roc_all_models.png", p_roc_all,
       width = 8, height = 7, dpi = 300, bg = "white")
cat("  Saved: results/models/roc_all_models.png\n")

# ==============================================================================
# 4. SELECT BEST MODEL
# ==============================================================================
cat("\n========== Selecting Best Model ==========\n")

best_name <- comp_table$Model[1]  # highest median AUC
best_model <- model_list[[best_name]]
cat("Best model:", best_name, "\n")
cat("CV AUC:", round(max(best_model$results$ROC), 4), "\n")
cat("Best hyperparameters:\n")
print(best_model$bestTune)

# ==============================================================================
# 5. FINE-TUNE BEST MODEL (narrow grid around best hyperparams)
# ==============================================================================
cat("\n========== Fine-Tuning Best Model ==========\n")

if (best_name == "RF") {
  best_mtry <- best_model$bestTune$mtry
  fine_grid <- data.frame(mtry = unique(pmax(1, c(best_mtry - 1, best_mtry, best_mtry + 1))))
  set.seed(SEED)
  final_model <- train(
    Group ~ ., data = train_data,
    method = "rf", metric = "ROC", trControl = fitControl,
    tuneGrid = fine_grid, ntree = 1000, importance = TRUE
  )
} else {
  final_model <- best_model
}

cat("Final model CV AUC:", round(max(final_model$results$ROC), 4), "\n")

# ==============================================================================
# 6. OPTIMAL THRESHOLD
# ==============================================================================
cat("\n========== Optimal Threshold ==========\n")

# Get CV predictions for the final model
pred_final <- final_model$pred
bt_final   <- final_model$bestTune
for (col_nm in names(bt_final)) pred_final <- pred_final[pred_final[[col_nm]] == bt_final[[col_nm]], ]

roc_final <- roc(pred_final$obs, pred_final$Tumor, levels = c("Normal", "Tumor"), quiet = TRUE)

# Youden's J
youden_coords <- coords(roc_final, "best", best.method = "youden",
                          ret = c("threshold", "sensitivity", "specificity"))
optimal_thr <- as.numeric(youden_coords$threshold[1])
cat("  Youden's J threshold:", round(optimal_thr, 4), "\n")
cat("    Sensitivity:", round(as.numeric(youden_coords$sensitivity[1]), 4), "\n")
cat("    Specificity:", round(as.numeric(youden_coords$specificity[1]), 4), "\n")

# High-sensitivity threshold (>= 95%)
all_coords <- coords(roc_final, x = "all",
                      ret = c("threshold", "sensitivity", "specificity"))
viable_95 <- all_coords[all_coords$sensitivity >= 0.95, ]
if (nrow(viable_95) > 0) {
  best_95 <- viable_95[which.max(viable_95$specificity), ]
  cat("  Screening threshold (Sens>=95%):", round(best_95$threshold, 4), "\n")
  cat("    Sensitivity:", round(best_95$sensitivity, 4), "\n")
  cat("    Specificity:", round(best_95$specificity, 4), "\n")
}

# ==============================================================================
# 7. CALIBRATION CURVE
# ==============================================================================
cat("\n========== Calibration Curve ==========\n")

cal_data <- data.frame(
  obs  = ifelse(pred_final$obs == "Tumor", 1, 0),
  pred = pred_final$Tumor
)
cal_data$bin <- cut(cal_data$pred, breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
cal_summary <- cal_data %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(pred),
    mean_obs  = mean(obs),
    n         = n(),
    .groups   = "drop"
  )

p_cal <- ggplot(cal_summary, aes(x = mean_pred, y = mean_obs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = n), color = "#E31A1C") +
  geom_line(color = "#E31A1C", linewidth = 0.8) +
  scale_size_continuous(range = c(2, 8)) +
  labs(title = "Calibration Curve", x = "Mean Predicted Probability",
       y = "Observed Proportion", size = "N") +
  my_theme
ggsave("results/models/calibration.png", p_cal, width = 8, height = 6, dpi = 300, bg = "white")
cat("  Saved: results/models/calibration.png\n")

# ==============================================================================
# 8. LEARNING CURVE
# ==============================================================================
cat("\n========== Learning Curve ==========\n")

fractions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
lc_results <- data.frame(fraction = fractions, n_samples = NA, AUC = NA)

for (i in seq_along(fractions)) {
  frac <- fractions[i]
  set.seed(SEED)
  if (frac < 1) {
    sub_idx <- createDataPartition(train_data$Group, p = frac, list = FALSE)
    sub_data <- train_data[sub_idx, ]
  } else {
    sub_data <- train_data
  }
  lc_results$n_samples[i] <- nrow(sub_data)

  lc_ctrl <- trainControl(
    method = "cv", number = 5,
    classProbs = TRUE, summaryFunction = twoClassSummary,
    savePredictions = "final", sampling = use_sampling
  )

  set.seed(SEED)
  lc_model <- train(
    Group ~ ., data = sub_data,
    method = "rf", metric = "ROC", trControl = lc_ctrl,
    tuneGrid = data.frame(mtry = final_model$bestTune$mtry),
    ntree = 500
  )
  lc_results$AUC[i] <- max(lc_model$results$ROC)
  cat(sprintf("  %.0f%%: n=%d, AUC=%.4f\n", frac * 100,
              lc_results$n_samples[i], lc_results$AUC[i]))
}

p_lc <- ggplot(lc_results, aes(x = n_samples, y = AUC)) +
  geom_line(linewidth = 1.2, color = "#1F78B4") +
  geom_point(size = 3, color = "#1F78B4") +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -1, size = 3.5) +
  ylim(c(min(lc_results$AUC) - 0.05, 1)) +
  labs(title = "Learning Curve", x = "Training Set Size", y = "CV AUC") +
  my_theme
ggsave("results/models/learning_curve.png", p_lc,
       width = 8, height = 6, dpi = 300, bg = "white")
cat("  Saved: results/models/learning_curve.png\n")

# ==============================================================================
# 9. SAVE FINAL OUTPUTS
# ==============================================================================
cat("\n========== Saving Outputs ==========\n")

saveRDS(final_model,  "results/models/final_model.rds")
saveRDS(top_feats,    "results/models/final_features.rds")
saveRDS(optimal_thr,  "results/models/optimal_threshold.rds")
saveRDS(model_list,   "results/models/all_candidate_models.rds")
saveRDS(comp_table,   "results/models/model_comparison.rds")

cat("  Saved: results/models/final_model.rds\n")
cat("  Saved: results/models/optimal_threshold.rds\n")
cat("  Saved: results/models/all_candidate_models.rds\n")

# CV AUC with CI
cv_auc <- as.numeric(auc(roc_final))
cv_ci  <- as.numeric(ci.auc(roc_final))
cat(sprintf("\nBest model: %s, CV AUC = %.4f (95%% CI: %.4f–%.4f)\n",
            best_name, cv_auc, cv_ci[1], cv_ci[3]))

# Check
if (cv_auc < 0.85) {
  cat("WARNING: CV AUC < 0.85 — model performance may be suboptimal.\n")
}

end_time <- Sys.time()
cat(sprintf("\n=== PHASE 2 complete. Runtime: %.1f minutes ===\n",
            as.numeric(difftime(end_time, start_time, units = "mins"))))
cat("=================================================================\n")
cat("  Next: run 04_independent_validation.R\n")
cat("=================================================================\n")

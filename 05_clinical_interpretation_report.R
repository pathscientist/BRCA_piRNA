################################################################################
#                                                                              #
#   05 — Clinical Interpretation & Final Report (PHASE 4)                      #
#                                                                              #
#   Signature characterisation table (Table 1), summary heatmap across all    #
#   7 cohorts, final text report, and session info.                           #
#                                                                              #
#   Requires: results from 01 + 02 + 03 + 04                                 #
#   Produces: results/final_signature_table.csv                               #
#             results/feature_selection/heatmap_all_cohorts.png               #
#             results/validation/final_report.txt                             #
#             results/session_info.txt                                        #
#                                                                              #
################################################################################

start_time <- Sys.time()
cat("=================================================================\n")
cat("  PHASE 4: Clinical Interpretation & Final Report\n")
cat("=================================================================\n\n")

# ==============================================================================
# 0. PACKAGES
# ==============================================================================
required_pkgs <- c(
  "ggplot2", "dplyr", "pheatmap", "viridis", "RColorBrewer",
  "pROC", "caret", "randomForest"
)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(viridis)
  library(RColorBrewer)
  library(pROC)
  library(caret)
  library(randomForest)
})
cat("All packages loaded.\n\n")

SEED <- 2024
set.seed(SEED)

# ==============================================================================
# 0.1 LOAD ALL PREVIOUS OUTPUTS
# ==============================================================================
cat("Loading all pipeline outputs...\n")

combat_df_all   <- readRDS("results/models/combat_df_all.rds")
train_df        <- readRDS("results/models/train_df.rds")
top_feats       <- readRDS("results/models/final_features.rds")
gene_cols       <- readRDS("results/models/common_pirnas.rds")
model           <- readRDS("results/models/final_model.rds")
optimal_thr     <- readRDS("results/models/optimal_threshold.rds")
val_details     <- readRDS("results/validation/validation_details.rds")

eval_v1    <- val_details$eval_v1
eval_v2    <- val_details$eval_v2
auc_cv     <- val_details$auc_cv
auc_cv_ci  <- val_details$auc_cv_ci

# Load feature selection summary if available
fs_table_path <- "results/feature_selection/fs_summary_table.csv"
if (file.exists(fs_table_path)) {
  fs_summary <- read.csv(fs_table_path, stringsAsFactors = FALSE)
} else {
  fs_summary <- NULL
}

cat("  Signature:", paste(top_feats, collapse = ", "), "\n")
cat("  Model loaded. Threshold:", round(optimal_thr, 4), "\n\n")

# ==============================================================================
# 1. SIGNATURE CHARACTERISATION TABLE (Table 1 of the paper)
# ==============================================================================
cat("========== 1. Signature Characterisation Table ==========\n")

x_train <- as.matrix(train_df[, gene_cols])
y_train <- train_df$Group

# limma DE results for the final signature
tryCatch({
  suppressPackageStartupMessages(library(limma))
  design <- model.matrix(~ y_train)
  fit <- lmFit(t(x_train), design)
  fit <- eBayes(fit)
  limma_res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
}, error = function(e) {
  cat("  limma not available; using Wilcoxon for DE stats.\n")
  limma_res <<- NULL
})

sig_table <- data.frame(piRNA_ID = top_feats, stringsAsFactors = FALSE)

# Selection frequency
if (!is.null(fs_summary)) {
  sig_table$N_methods <- fs_summary$N_methods[match(sig_table$piRNA_ID, fs_summary$piRNA)]
  sig_table$which_methods <- fs_summary$which_methods[match(sig_table$piRNA_ID, fs_summary$piRNA)]
} else {
  sig_table$N_methods <- NA
  sig_table$which_methods <- NA
}

# Direction (up/down in tumor)
med_tumor  <- apply(x_train[y_train == "Tumor", top_feats, drop = FALSE], 2, median)
med_normal <- apply(x_train[y_train == "Normal", top_feats, drop = FALSE], 2, median)
sig_table$Direction <- ifelse(med_tumor > med_normal, "Up", "Down")

# limma stats
if (exists("limma_res") && !is.null(limma_res)) {
  sig_table$log2FC <- limma_res[top_feats, "logFC"]
  sig_table$adj_P  <- limma_res[top_feats, "adj.P.Val"]
} else {
  sig_table$log2FC <- med_tumor - med_normal
  pvals <- sapply(top_feats, function(f) {
    suppressWarnings(wilcox.test(x_train[, f] ~ y_train)$p.value)
  })
  sig_table$adj_P <- p.adjust(pvals, method = "BH")
}

# Individual AUC (univariate, on discovery set, 5-fold CV)
set.seed(SEED)
folds <- createFolds(y_train, k = 5, returnTrain = TRUE)
sig_table$Individual_AUC <- sapply(top_feats, function(feat) {
  aucs <- numeric(length(folds))
  for (i in seq_along(folds)) {
    tr_idx <- folds[[i]]; te_idx <- setdiff(seq_along(y_train), tr_idx)
    set.seed(SEED)
    rf_tmp <- randomForest(x = x_train[tr_idx, feat, drop = FALSE], y = y_train[tr_idx],
                           ntree = 500)
    prob <- predict(rf_tmp, x_train[te_idx, feat, drop = FALSE], type = "prob")[, "Tumor"]
    y01 <- ifelse(y_train[te_idx] == "Tumor", 1, 0)
    if (length(unique(y01)) < 2) { aucs[i] <- 0.5; next }
    aucs[i] <- as.numeric(auc(roc(y01, prob, direction = "<", quiet = TRUE)))
  }
  mean(aucs)
})

# RF importance rank
set.seed(SEED)
rf_full <- randomForest(x = x_train[, top_feats, drop = FALSE], y = y_train,
                        ntree = 1000, importance = TRUE)
rf_imp <- importance(rf_full, type = 2)  # MeanDecreaseGini
sig_table$RF_importance_rank <- rank(-rf_imp[match(top_feats, rownames(rf_imp)), 1])

# Known associations
sig_table$Known_associations <- "see literature"

# Round for display
sig_table$log2FC         <- round(sig_table$log2FC, 3)
sig_table$adj_P          <- signif(sig_table$adj_P, 3)
sig_table$Individual_AUC <- round(sig_table$Individual_AUC, 3)

cat("\nSignature Characterisation Table:\n")
print(sig_table)

write.csv(sig_table, "results/final_signature_table.csv", row.names = FALSE)
cat("\n  Saved: results/final_signature_table.csv\n")

# ==============================================================================
# 2. SUMMARY HEATMAP (all 7 cohorts)
# ==============================================================================
cat("\n========== 2. Summary Heatmap ==========\n")

# Expression matrix: rows = piRNAs, columns = samples
heat_mat <- t(as.matrix(combat_df_all[, top_feats]))

# Annotation columns
ann_col <- data.frame(
  Group    = combat_df_all$Group,
  Cohort   = combat_df_all$Batch,
  Specimen = ifelse(combat_df_all$Batch == "yyfbatch2", "Plasma/Serum", "Tissue"),
  row.names = rownames(combat_df_all)
)

# Annotation colours
ann_colors <- list(
  Group    = c("Normal" = "#4393C3", "Tumor" = "#D6604D"),
  Cohort   = setNames(
    c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3")[
      seq_len(length(unique(combat_df_all$Batch)))],
    sort(unique(combat_df_all$Batch))
  ),
  Specimen = c("Tissue" = "#FFD700", "Plasma/Serum" = "#20B2AA")
)

# Cap extreme values for cleaner heatmap
heat_mat_capped <- heat_mat
heat_mat_capped[heat_mat_capped > 3]  <- 3
heat_mat_capped[heat_mat_capped < -3] <- -3

png("results/feature_selection/heatmap_all_cohorts.png",
    width = 12, height = 8, units = "in", res = 300)
pheatmap(
  heat_mat_capped,
  annotation_col    = ann_col,
  annotation_colors = ann_colors,
  show_colnames     = FALSE,
  clustering_method = "ward.D2",
  color             = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  main              = paste(length(top_feats), "piRNA Diagnostic Signature — All Cohorts"),
  fontsize          = 10,
  fontsize_row      = 9,
  border_color      = NA
)
dev.off()
cat("  Saved: results/feature_selection/heatmap_all_cohorts.png\n")

# ==============================================================================
# 3. FINAL REPORT
# ==============================================================================
cat("\n========== 3. Final Report ==========\n")

# Compute discovery CV sensitivity/specificity at optimal threshold
pred_cv <- model$pred
bt <- model$bestTune
for (nm in names(bt)) pred_cv <- pred_cv[pred_cv[[nm]] == bt[[nm]], ]
cv_pred_class <- factor(ifelse(pred_cv$Tumor >= optimal_thr, "Tumor", "Normal"),
                        levels = c("Normal", "Tumor"))
cv_cm <- confusionMatrix(cv_pred_class, pred_cv$obs, positive = "Tumor")
cv_sens <- as.numeric(cv_cm$byClass["Sensitivity"])
cv_spec <- as.numeric(cv_cm$byClass["Specificity"])

# Best model name
best_name <- model$method
best_params <- paste(names(model$bestTune), "=", model$bestTune, collapse = ", ")

report_lines <- c(
  "===================================================================",
  "  piRNA Diagnostic Signature — Final Report",
  "===================================================================",
  "",
  sprintf("  Signature: %s", paste(top_feats, collapse = ", ")),
  sprintf("  Number of piRNAs: %d", length(top_feats)),
  sprintf("  Model: %s with %s", best_name, best_params),
  sprintf("  Optimal threshold: %.4f (Youden's J)", optimal_thr),
  "",
  "  -----------------------------------------------------------",
  "  Discovery (5 public tissue cohorts, 10x5 CV):",
  sprintf("    AUC = %.4f (95%% CI: %.4f–%.4f)", auc_cv, auc_cv_ci[1], auc_cv_ci[3]),
  sprintf("    Sensitivity = %.1f%%, Specificity = %.1f%%", cv_sens * 100, cv_spec * 100),
  "",
  sprintf("  Tissue Validation (yyfbatch1, n=%d):", eval_v1$n),
  sprintf("    AUC = %.4f (95%% CI: %.4f–%.4f)", eval_v1$auc, eval_v1$auc_lower, eval_v1$auc_upper),
  sprintf("    Sensitivity = %.1f%%, Specificity = %.1f%%", eval_v1$sens * 100, eval_v1$spec * 100),
  "",
  sprintf("  Plasma Validation (yyfbatch2, n=%d):", eval_v2$n),
  sprintf("    AUC = %.4f (95%% CI: %.4f–%.4f)", eval_v2$auc, eval_v2$auc_lower, eval_v2$auc_upper),
  sprintf("    Sensitivity = %.1f%%, Specificity = %.1f%%", eval_v2$sens * 100, eval_v2$spec * 100),
  "",
  "  -----------------------------------------------------------",
  sprintf("  Tissue vs Plasma AUC gap: %.4f", eval_v1$auc - eval_v2$auc),
  sprintf("  Overfitting assessment:"),
  sprintf("    Delta(train-tissue) = %.4f", auc_cv - eval_v1$auc),
  sprintf("    Delta(train-plasma) = %.4f", auc_cv - eval_v2$auc),
  "===================================================================",
  "",
  sprintf("  Report generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
)

cat("\n")
writeLines(report_lines)
writeLines(report_lines, con = "results/validation/final_report.txt")
cat("\n  Saved: results/validation/final_report.txt\n")

# ==============================================================================
# 4. SESSION INFO
# ==============================================================================
cat("\n========== 4. Session Info ==========\n")

end_time <- Sys.time()
total_runtime <- as.numeric(difftime(end_time, start_time, units = "mins"))

session_lines <- c(
  "===================================================================",
  "  Session Information",
  "===================================================================",
  "",
  sprintf("  Phase 4 runtime: %.1f minutes", total_runtime),
  sprintf("  Report date: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  capture.output(sessionInfo())
)

writeLines(session_lines, con = "results/session_info.txt")
cat("  Saved: results/session_info.txt\n")

cat(sprintf("\n=== PHASE 4 complete. Runtime: %.1f minutes ===\n", total_runtime))
cat("=================================================================\n")
cat("  PIPELINE COMPLETE.\n")
cat("  \n")
cat("  Output files:\n")
cat("    results/final_signature_table.csv\n")
cat("    results/feature_selection/heatmap_all_cohorts.png\n")
cat("    results/validation/final_report.txt\n")
cat("    results/session_info.txt\n")
cat("=================================================================\n")

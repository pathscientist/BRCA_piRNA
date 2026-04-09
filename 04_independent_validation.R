################################################################################
#                                                                              #
#   04 — Independent Validation: Plasma ROC (PHASE 3)                         #
#                                                                              #
#   Evaluate final model on yyfbatch2 (independent plasma hold-out).          #
#   Also evaluate per-cohort (BRCA1, yyfbatch1) on training data to show      #
#   per-batch performance.                                                     #
#                                                                              #
#   Training:   BRCA1 (tissue) + yyfbatch1 (plasma)                           #
#   Validation: yyfbatch2 (plasma — NEVER seen during training)               #
#   Both yyfbatch1 and yyfbatch2 are plasma.                                   #
#                                                                              #
#   Requires: results from 01 + 02 + 03                                       #
#   Produces: results/validation/*                                             #
#                                                                              #
################################################################################

start_time <- Sys.time()
cat("=================================================================\n")
cat("  PHASE 3: Independent Validation — Plasma Hold-Out\n")
cat("=================================================================\n\n")

# ==============================================================================
# 0. PACKAGES
# ==============================================================================
required_pkgs <- c(
  "caret", "randomForest", "pROC", "PRROC",
  "ggplot2", "dplyr", "cowplot", "patchwork", "scales"
)
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
}
suppressPackageStartupMessages({
  library(caret)
  library(randomForest)
  library(pROC)
  library(PRROC)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(patchwork)
  library(scales)
})
cat("All packages loaded.\n\n")

SEED <- 2024
set.seed(SEED)

COLOR_DISCOVERY <- "#1F78B4"
COLOR_PLASMA_V  <- "#33A02C"
COLOR_BRCA1     <- "#E31A1C"
COLOR_YYF1      <- "#FF7F00"
COLOR_PALETTE   <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00",
                      "#6A3D9A", "#B15928", "#A6CEE3")

my_theme <- theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom",
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

# ==============================================================================
# 0.1 LOAD OUTPUTS FROM PHASES 0-2
# ==============================================================================
cat("Loading previous outputs...\n")
combat_df_all    <- readRDS("results/models/combat_df_all.rds")
train_df         <- readRDS("results/models/train_df.rds")
valid_yyfbatch2  <- readRDS("results/models/valid_yyfbatch2.rds")
top_feats        <- readRDS("results/models/final_features.rds")
model            <- readRDS("results/models/final_model.rds")
optimal_thr      <- readRDS("results/models/optimal_threshold.rds")

cat("  Signature:", paste(top_feats, collapse = ", "), "\n")
cat("  Optimal threshold:", round(optimal_thr, 4), "\n")
cat("  Training cohorts:", paste(unique(train_df$Batch), collapse = ", "), "\n")
cat("  Validation: yyfbatch2 (plasma)\n\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

evaluate_cohort <- function(df, model, features, threshold, cohort_name) {
  prob <- predict(model, df[, features, drop = FALSE], type = "prob")$Tumor
  y_true <- df$Group
  y_bin  <- ifelse(y_true == "Tumor", 1, 0)

  roc_obj <- roc(y_bin, prob, direction = "<", quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  auc_ci  <- as.numeric(ci.auc(roc_obj))

  pred_class <- factor(ifelse(prob >= threshold, "Tumor", "Normal"),
                       levels = c("Normal", "Tumor"))
  cm <- confusionMatrix(pred_class, y_true, positive = "Tumor")

  sens <- as.numeric(cm$byClass["Sensitivity"])
  spec <- as.numeric(cm$byClass["Specificity"])
  ppv  <- as.numeric(cm$byClass["Pos Pred Value"])
  npv  <- as.numeric(cm$byClass["Neg Pred Value"])
  acc  <- as.numeric(cm$overall["Accuracy"])
  f1   <- ifelse((ppv + sens) == 0 | is.na(ppv) | is.na(sens), NA,
                 2 * ppv * sens / (ppv + sens))

  pr_obj <- tryCatch(
    pr.curve(scores.class0 = prob[y_bin == 1],
             scores.class1 = prob[y_bin == 0], curve = TRUE),
    error = function(e) NULL
  )

  list(
    cohort = cohort_name, n = nrow(df),
    prob = prob, y_true = y_true, y_bin = y_bin,
    roc_obj = roc_obj, pr_obj = pr_obj,
    auc = auc_val, auc_lower = auc_ci[1], auc_upper = auc_ci[3],
    sens = sens, spec = spec, ppv = ppv, npv = npv, acc = acc, f1 = f1,
    cm = cm
  )
}

make_roc_df <- function(roc_obj, label) {
  data.frame(fpr = 1 - roc_obj$specificities,
             tpr = roc_obj$sensitivities,
             Set = label, stringsAsFactors = FALSE)
}

plot_roc <- function(roc_obj, title, label, color) {
  df <- data.frame(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)
  df <- df[order(df$fpr), ]
  ggplot(df, aes(fpr, tpr)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_path(color = color, linewidth = 1.2) +
    annotate("label", x = 0.95, y = 0.05, label = label,
             hjust = 1, vjust = 0, size = 4, fill = "white", alpha = 0.85,
             label.size = NA) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
    labs(title = title, x = "False Positive Rate (1-Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    my_theme
}

# ==============================================================================
# 1. DISCOVERY CV METRICS
# ==============================================================================
cat("========== 1. Discovery CV Metrics ==========\n")

pred_cv <- model$pred
bt <- model$bestTune
for (nm in names(bt)) pred_cv <- pred_cv[pred_cv[[nm]] == bt[[nm]], ]

y_true_cv  <- ifelse(pred_cv$obs == "Tumor", 1, 0)
y_score_cv <- pred_cv$Tumor
roc_cv     <- roc(y_true_cv, y_score_cv, direction = "<", quiet = TRUE)
auc_cv     <- as.numeric(auc(roc_cv))
auc_cv_ci  <- as.numeric(ci.auc(roc_cv))

cat(sprintf("  Discovery CV AUC: %.4f (95%% CI: %.4f–%.4f)\n",
            auc_cv, auc_cv_ci[1], auc_cv_ci[3]))

# ==============================================================================
# 2. INDEPENDENT VALIDATION — yyfbatch2 (PLASMA)
# ==============================================================================
cat("\n========== 2. Independent Validation: yyfbatch2 (Plasma) ==========\n")

eval_v2 <- evaluate_cohort(valid_yyfbatch2, model, top_feats, optimal_thr, "yyfbatch2")

cat(sprintf("  yyfbatch2 (plasma):  AUC=%.4f (%.4f–%.4f)\n",
            eval_v2$auc, eval_v2$auc_lower, eval_v2$auc_upper))
cat(sprintf("    Sensitivity=%.1f%%, Specificity=%.1f%%\n",
            eval_v2$sens * 100, eval_v2$spec * 100))
cat(sprintf("    PPV=%.1f%%, NPV=%.1f%%, F1=%.3f\n",
            eval_v2$ppv * 100, eval_v2$npv * 100, eval_v2$f1))

# Confusion matrix
tab <- as.data.frame(eval_v2$cm$table)
colnames(tab) <- c("Prediction", "Reference", "N")
tab <- tab %>% group_by(Reference) %>%
  mutate(Pct = sprintf("%.1f%%", 100 * N / sum(N))) %>% ungroup()
p_cm <- ggplot(tab, aes(x = Prediction, y = Reference, fill = N)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = paste0(N, "\n(", Pct, ")")), size = 5, fontface = "bold") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_y_discrete(limits = rev(c("Normal", "Tumor"))) +
  labs(title = "Confusion Matrix — yyfbatch2 (Plasma, Independent)",
       subtitle = sprintf("Threshold = %.3f", optimal_thr),
       x = "Predicted", y = "Actual") +
  my_theme + theme(legend.position = "none")
ggsave("results/validation/yyfbatch2_confusion.png", p_cm,
       width = 6, height = 5, dpi = 300, bg = "white")
cat("  Saved: results/validation/yyfbatch2_confusion.png\n")

# ==============================================================================
# 3. KEY FIGURE: 3-PANEL ROC
# ==============================================================================
cat("\n========== 3. Key Figure: ROC Panels ==========\n")

# Panel A: Discovery CV (BRCA1 + yyfbatch1 combined)
lab_a <- sprintf("AUC = %.3f\n(95%% CI: %.3f–%.3f)", auc_cv, auc_cv_ci[1], auc_cv_ci[3])
p_a <- plot_roc(roc_cv, "A: Discovery CV (BRCA1 + yyfbatch1)", lab_a, COLOR_DISCOVERY)

# Panel B: Independent Plasma Validation (yyfbatch2)
lab_b <- sprintf("AUC = %.3f\n(95%% CI: %.3f–%.3f)",
                 eval_v2$auc, eval_v2$auc_lower, eval_v2$auc_upper)
p_b <- plot_roc(eval_v2$roc_obj, "B: Independent Validation (yyfbatch2, Plasma)",
                lab_b, COLOR_PLASMA_V)

# Panel C: Overlay
roc_overlay <- rbind(
  make_roc_df(roc_cv,          sprintf("Discovery AUC=%.3f", auc_cv)),
  make_roc_df(eval_v2$roc_obj, sprintf("yyfbatch2 AUC=%.3f", eval_v2$auc))
)
roc_overlay$Set <- factor(roc_overlay$Set, levels = unique(roc_overlay$Set))

p_c <- ggplot(roc_overlay, aes(fpr, tpr, color = Set)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_path(linewidth = 1.1) +
  scale_color_manual(values = c(COLOR_DISCOVERY, COLOR_PLASMA_V)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(title = "C: Discovery vs Independent Validation",
       x = "False Positive Rate", y = "True Positive Rate", color = "") +
  my_theme +
  theme(legend.position = c(0.65, 0.2),
        legend.background = element_rect(fill = alpha("white", 0.9), color = "grey70"))

p_3panel <- (p_a | p_b) / (p_c + plot_spacer()) +
  plot_annotation(title = paste(length(top_feats), "piRNA Diagnostic Signature"),
                  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16)))

ggsave("results/validation/Figure_ROC_validation.png", p_3panel,
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave("results/validation/Figure_ROC_validation.pdf", p_3panel,
       width = 10, height = 8)
cat("  Saved: results/validation/Figure_ROC_validation.png\n")
cat("  Saved: results/validation/Figure_ROC_validation.pdf\n")

# ==============================================================================
# 4. PER-COHORT INDIVIDUAL PLOTS (yyfbatch2)
# ==============================================================================
cat("\n========== 4. Per-Cohort Individual Plots ==========\n")

# ROC
p_roc_v2 <- plot_roc(
  eval_v2$roc_obj, "ROC — yyfbatch2 (Plasma)",
  sprintf("AUC = %.3f (%.3f–%.3f)", eval_v2$auc, eval_v2$auc_lower, eval_v2$auc_upper),
  COLOR_PLASMA_V
)
ggsave("results/validation/yyfbatch2_roc.png", p_roc_v2,
       width = 7, height = 6, dpi = 300, bg = "white")

# PRC
if (!is.null(eval_v2$pr_obj)) {
  pr_df <- data.frame(recall = eval_v2$pr_obj$curve[, 1],
                      precision = eval_v2$pr_obj$curve[, 2])
  p_prc <- ggplot(pr_df, aes(recall, precision)) +
    geom_hline(yintercept = mean(eval_v2$y_bin), linetype = "dashed", color = "grey50") +
    geom_path(color = COLOR_PLASMA_V, linewidth = 1.2) +
    annotate("label", x = 0.95, y = 0.05,
             label = sprintf("AUPRC = %.3f", eval_v2$pr_obj$auc.integral),
             hjust = 1, vjust = 0, size = 4, fill = "white") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
    labs(title = "Precision-Recall — yyfbatch2 (Plasma)",
         x = "Recall", y = "Precision") +
    my_theme
  ggsave("results/validation/yyfbatch2_prc.png", p_prc,
         width = 7, height = 6, dpi = 300, bg = "white")
}

# Probability histogram
hist_df <- data.frame(prob = eval_v2$prob, Group = eval_v2$y_true)
hist_df$uncertain <- ifelse(hist_df$prob >= 0.4 & hist_df$prob <= 0.6,
                            "Uncertain", "Confident")
p_hist <- ggplot(hist_df, aes(x = prob, fill = Group)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = optimal_thr, linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(0.4, 0.6), linetype = "dotted", color = "grey50") +
  scale_fill_manual(values = c("Normal" = "#4393C3", "Tumor" = "#D6604D")) +
  labs(title = "Prediction Probability — yyfbatch2 (Plasma)",
       subtitle = sprintf("Uncertain zone [0.4–0.6]: %d samples (%.1f%%)",
                          sum(hist_df$uncertain == "Uncertain"),
                          100 * mean(hist_df$uncertain == "Uncertain")),
       x = "Predicted P(Tumor)", y = "Count") +
  my_theme
ggsave("results/validation/yyfbatch2_probability_hist.png", p_hist,
       width = 7, height = 6, dpi = 300, bg = "white")
cat("  Saved ROC, PRC, probability histogram for yyfbatch2\n")

# ==============================================================================
# 5. PER-BATCH ROC OVERLAY (all available cohorts)
# ==============================================================================
cat("\n========== 5. Per-Batch ROC Overlay ==========\n")

all_batches <- unique(combat_df_all$Batch)
roc_per_batch <- list()

for (b in all_batches) {
  sub <- combat_df_all[combat_df_all$Batch == b, ]
  if (nrow(sub) < 5) next
  prob <- predict(model, sub[, top_feats, drop = FALSE], type = "prob")$Tumor
  y01 <- ifelse(sub$Group == "Tumor", 1, 0)
  if (length(unique(y01)) < 2) next
  roc_obj <- roc(y01, prob, direction = "<", quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  specimen <- ifelse(b == "BRCA1", "Tissue", "Plasma")
  label <- sprintf("%s [%s] AUC=%.3f", b, specimen, auc_val)
  roc_per_batch[[b]] <- make_roc_df(roc_obj, label)
}
roc_per_df <- do.call(rbind, roc_per_batch)
roc_per_df$Set <- factor(roc_per_df$Set, levels = unique(roc_per_df$Set))

n_b <- length(unique(roc_per_df$Set))
p_all_batch <- ggplot(roc_per_df, aes(fpr, tpr, color = Set)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_path(linewidth = 0.9) +
  scale_color_manual(values = COLOR_PALETTE[seq_len(n_b)]) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(title = "ROC — All Cohorts (Individual Prediction)",
       subtitle = paste(length(top_feats), "piRNA signature"),
       x = "False Positive Rate", y = "True Positive Rate", color = "") +
  my_theme +
  theme(legend.position = "right", legend.text = element_text(size = 9))
ggsave("results/validation/roc_all_cohorts.png", p_all_batch,
       width = 9, height = 7, dpi = 300, bg = "white")
cat("  Saved: results/validation/roc_all_cohorts.png\n")

# ==============================================================================
# 6. ROBUSTNESS CHECKS (yyfbatch2 only)
# ==============================================================================
cat("\n========== 6. Robustness Checks — yyfbatch2 ==========\n")

# 6a. Bootstrap (1000 resamples)
cat("  Bootstrap AUC (1000 resamples)...\n")
set.seed(SEED)
boot_aucs <- numeric(1000)
for (i in seq_len(1000)) {
  idx <- sample(length(eval_v2$y_bin), replace = TRUE)
  yt <- eval_v2$y_bin[idx]; ys <- eval_v2$prob[idx]
  if (length(unique(yt)) < 2) { boot_aucs[i] <- NA; next }
  boot_aucs[i] <- as.numeric(auc(roc(yt, ys, direction = "<", quiet = TRUE)))
}
boot_aucs <- na.omit(boot_aucs)
boot_ci <- quantile(boot_aucs, c(0.025, 0.975))
cat(sprintf("    Bootstrap median AUC: %.4f, 95%% CI: [%.4f, %.4f]\n",
            median(boot_aucs), boot_ci[1], boot_ci[2]))

p_boot <- ggplot(data.frame(AUC = boot_aucs), aes(AUC)) +
  geom_histogram(bins = 50, fill = "#1F78B4", alpha = 0.7) +
  geom_vline(xintercept = eval_v2$auc, color = "#E31A1C", linewidth = 1) +
  geom_vline(xintercept = boot_ci, linetype = "dashed", color = "#E31A1C") +
  labs(title = "Bootstrap AUC — yyfbatch2 (Plasma)",
       subtitle = sprintf("Median=%.3f, 95%% CI=[%.3f, %.3f]",
                          median(boot_aucs), boot_ci[1], boot_ci[2]),
       x = "AUC", y = "Count") +
  my_theme
ggsave("results/validation/yyfbatch2_bootstrap.png", p_boot,
       width = 7, height = 5, dpi = 300, bg = "white")

# 6b. Permutation test (1000 label shuffles)
cat("  Permutation test (1000 shuffles)...\n")
set.seed(SEED)
null_aucs <- numeric(1000)
for (i in seq_len(1000)) {
  perm_y <- sample(eval_v2$y_bin)
  if (length(unique(perm_y)) < 2) { null_aucs[i] <- 0.5; next }
  null_aucs[i] <- as.numeric(auc(roc(perm_y, eval_v2$prob, direction = "<", quiet = TRUE)))
}
perm_p <- mean(null_aucs >= eval_v2$auc)
cat(sprintf("    Permutation p-value: %.4f\n", perm_p))

p_perm <- ggplot(data.frame(AUC = null_aucs), aes(AUC)) +
  geom_histogram(bins = 50, fill = "grey70", alpha = 0.7) +
  geom_vline(xintercept = eval_v2$auc, color = "#E31A1C", linewidth = 1.2) +
  labs(title = "Permutation Test — yyfbatch2 (Plasma)",
       subtitle = sprintf("Observed AUC=%.3f, p=%.4f", eval_v2$auc, perm_p),
       x = "Null AUC Distribution", y = "Count") +
  my_theme
ggsave("results/validation/yyfbatch2_permutation.png", p_perm,
       width = 7, height = 5, dpi = 300, bg = "white")
cat("  Saved bootstrap + permutation plots\n")

# ==============================================================================
# 7. PERFORMANCE DROP ASSESSMENT
# ==============================================================================
cat("\n========== 7. Performance Drop Assessment ==========\n")

delta_val <- auc_cv - eval_v2$auc

cat(sprintf("  Discovery CV AUC:     %.4f\n", auc_cv))
cat(sprintf("  yyfbatch2 AUC:        %.4f\n", eval_v2$auc))
cat(sprintf("  Delta(train–yyfbatch2): %.4f\n", delta_val))

if (delta_val > 0.10) cat("  SEVERE WARNING: AUC drop > 0.10!\n")
else if (delta_val > 0.05) cat("  WARNING: AUC drop > 0.05 — potential overfitting.\n")
else cat("  OK: AUC drop within acceptable range.\n")

# AUC bar chart
auc_bar_df <- data.frame(
  Cohort = c("Discovery CV\n(BRCA1+yyfbatch1)", "yyfbatch2\n(Plasma, Independent)"),
  AUC    = c(auc_cv, eval_v2$auc),
  Lower  = c(auc_cv_ci[1], eval_v2$auc_lower),
  Upper  = c(auc_cv_ci[3], eval_v2$auc_upper),
  Type   = c("Discovery", "Validation")
)
auc_bar_df$Cohort <- factor(auc_bar_df$Cohort, levels = auc_bar_df$Cohort)

p_bar <- ggplot(auc_bar_df, aes(x = Cohort, y = AUC, fill = Type)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.5, size = 4.5) +
  scale_fill_manual(values = c("Discovery" = COLOR_DISCOVERY,
                                "Validation" = COLOR_PLASMA_V)) +
  ylim(c(0, 1.1)) +
  labs(title = "AUC Comparison: Discovery vs Validation", x = "", y = "AUC") +
  my_theme + theme(legend.position = "none")
ggsave("results/validation/auc_comparison_bar.png", p_bar,
       width = 7, height = 6, dpi = 300, bg = "white")
cat("  Saved: results/validation/auc_comparison_bar.png\n")

# ==============================================================================
# 8. RESULTS TABLE
# ==============================================================================
cat("\n========== 8. Summary Results Table ==========\n")

results_table <- data.frame(
  Cohort  = c("Discovery CV", "yyfbatch2"),
  Type    = c("Tissue+Plasma", "Plasma (Independent)"),
  N       = c(nrow(pred_cv), eval_v2$n),
  AUC     = c(auc_cv, eval_v2$auc),
  AUC_CI  = c(sprintf("%.3f-%.3f", auc_cv_ci[1], auc_cv_ci[3]),
              sprintf("%.3f-%.3f", eval_v2$auc_lower, eval_v2$auc_upper)),
  Sens    = c(NA, eval_v2$sens),
  Spec    = c(NA, eval_v2$spec),
  PPV     = c(NA, eval_v2$ppv),
  NPV     = c(NA, eval_v2$npv),
  stringsAsFactors = FALSE
)

cat("\n")
print(results_table)

write.csv(results_table, "results/validation/validation_results.csv", row.names = FALSE)
cat("\n  Saved: results/validation/validation_results.csv\n")

# Save key objects for Phase 4
saveRDS(list(eval_v2 = eval_v2,
             auc_cv = auc_cv, auc_cv_ci = auc_cv_ci,
             roc_cv = roc_cv, optimal_thr = optimal_thr,
             boot_ci = boot_ci, perm_p = perm_p),
        "results/validation/validation_details.rds")

end_time <- Sys.time()
cat(sprintf("\n=== PHASE 3 complete. Runtime: %.1f minutes ===\n",
            as.numeric(difftime(end_time, start_time, units = "mins"))))
cat("=================================================================\n")
cat("  Next: run 05_clinical_interpretation_report.R\n")
cat("=================================================================\n")

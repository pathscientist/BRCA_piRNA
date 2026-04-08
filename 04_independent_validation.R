################################################################################
#                                                                              #
#   04 — Independent Validation: Tissue vs Plasma ROC (PHASE 3)               #
#                                                                              #
#   Evaluate final model on yyfbatch1 (tissue) and yyfbatch2 (plasma/serum)   #
#   SEPARATELY. Per-cohort ROC, combined 4-panel ROC figure, per-cohort       #
#   granular ROC, robustness checks (bootstrap + permutation),                #
#   performance drop assessment.                                               #
#                                                                              #
#   Requires: results from 01 + 02 + 03                                       #
#   Produces: results/validation/*                                             #
#                                                                              #
################################################################################

start_time <- Sys.time()
cat("=================================================================\n")
cat("  PHASE 3: Independent Validation — Tissue vs Plasma\n")
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
COLOR_TISSUE    <- "#E31A1C"
COLOR_PLASMA    <- "#33A02C"
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
valid_yyfbatch1  <- readRDS("results/models/valid_yyfbatch1.rds")
valid_yyfbatch2  <- readRDS("results/models/valid_yyfbatch2.rds")
top_feats        <- readRDS("results/models/final_features.rds")
model            <- readRDS("results/models/final_model.rds")
optimal_thr      <- readRDS("results/models/optimal_threshold.rds")

cat("  Signature:", paste(top_feats, collapse = ", "), "\n")
cat("  Optimal threshold:", round(optimal_thr, 4), "\n\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Compute full metrics for a validation set
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

  # PRC
  pr_obj <- tryCatch(
    pr.curve(scores.class0 = prob[y_bin == 1],
             scores.class1 = prob[y_bin == 0], curve = TRUE),
    error = function(e) NULL
  )

  list(
    cohort    = cohort_name,
    n         = nrow(df),
    prob      = prob,
    y_true    = y_true,
    y_bin     = y_bin,
    roc_obj   = roc_obj,
    pr_obj    = pr_obj,
    auc       = auc_val,
    auc_lower = auc_ci[1],
    auc_upper = auc_ci[3],
    sens = sens, spec = spec, ppv = ppv, npv = npv, acc = acc, f1 = f1,
    cm = cm
  )
}

# Build ROC data frame for plotting
make_roc_df <- function(roc_obj, label) {
  data.frame(fpr = 1 - roc_obj$specificities,
             tpr = roc_obj$sensitivities,
             Set = label, stringsAsFactors = FALSE)
}

# Plot single ROC with CI band
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
# 1. DISCOVERY CV METRICS (for comparison)
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
# 2. PER-COHORT EVALUATION
# ==============================================================================
cat("\n========== 2. Per-Cohort Evaluation ==========\n")

eval_v1 <- evaluate_cohort(valid_yyfbatch1, model, top_feats, optimal_thr, "yyfbatch1")
eval_v2 <- evaluate_cohort(valid_yyfbatch2, model, top_feats, optimal_thr, "yyfbatch2")

cat(sprintf("  yyfbatch1 (tissue):  AUC=%.4f (%.4f–%.4f), Sens=%.1f%%, Spec=%.1f%%\n",
            eval_v1$auc, eval_v1$auc_lower, eval_v1$auc_upper,
            eval_v1$sens * 100, eval_v1$spec * 100))
cat(sprintf("  yyfbatch2 (plasma):  AUC=%.4f (%.4f–%.4f), Sens=%.1f%%, Spec=%.1f%%\n",
            eval_v2$auc, eval_v2$auc_lower, eval_v2$auc_upper,
            eval_v2$sens * 100, eval_v2$spec * 100))

# Confusion matrices
for (ev in list(eval_v1, eval_v2)) {
  tab <- as.data.frame(ev$cm$table)
  colnames(tab) <- c("Prediction", "Reference", "N")
  tab <- tab %>% group_by(Reference) %>%
    mutate(Pct = sprintf("%.1f%%", 100 * N / sum(N))) %>% ungroup()
  p_cm <- ggplot(tab, aes(x = Prediction, y = Reference, fill = N)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = paste0(N, "\n(", Pct, ")")), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    scale_y_discrete(limits = rev(c("Normal", "Tumor"))) +
    labs(title = paste("Confusion Matrix —", ev$cohort),
         subtitle = sprintf("Threshold = %.3f", optimal_thr),
         x = "Predicted", y = "Actual") +
    my_theme + theme(legend.position = "none")
  ggsave(sprintf("results/validation/%s_confusion.png", ev$cohort),
         p_cm, width = 6, height = 5, dpi = 300, bg = "white")
  cat(sprintf("  Saved: results/validation/%s_confusion.png\n", ev$cohort))
}

# ==============================================================================
# 3. KEY FIGURE: 4-PANEL ROC (2x2 grid)
# ==============================================================================
cat("\n========== 3. Key Figure: 4-Panel ROC ==========\n")

# Panel A: Discovery CV
lab_a <- sprintf("AUC = %.3f\n(95%% CI: %.3f–%.3f)", auc_cv, auc_cv_ci[1], auc_cv_ci[3])
p_a <- plot_roc(roc_cv, "A: Discovery (5-Cohort Training, 10x5 CV)", lab_a, COLOR_DISCOVERY)

# Panel B: Tissue Validation (yyfbatch1)
lab_b <- sprintf("AUC = %.3f\n(95%% CI: %.3f–%.3f)",
                 eval_v1$auc, eval_v1$auc_lower, eval_v1$auc_upper)
p_b <- plot_roc(eval_v1$roc_obj, "B: Tissue Validation (yyfbatch1)", lab_b, COLOR_TISSUE)

# Panel C: Plasma Validation (yyfbatch2)
lab_c <- sprintf("AUC = %.3f\n(95%% CI: %.3f–%.3f)",
                 eval_v2$auc, eval_v2$auc_lower, eval_v2$auc_upper)
p_c <- plot_roc(eval_v2$roc_obj, "C: Plasma Validation (yyfbatch2)", lab_c, COLOR_PLASMA)

# Panel D: All Cohorts Overlay
roc_overlay <- rbind(
  make_roc_df(roc_cv,          sprintf("Discovery AUC=%.3f", auc_cv)),
  make_roc_df(eval_v1$roc_obj, sprintf("Tissue AUC=%.3f", eval_v1$auc)),
  make_roc_df(eval_v2$roc_obj, sprintf("Plasma AUC=%.3f", eval_v2$auc))
)
roc_overlay$Set <- factor(roc_overlay$Set, levels = unique(roc_overlay$Set))

p_d <- ggplot(roc_overlay, aes(fpr, tpr, color = Set)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_path(linewidth = 1.1) +
  scale_color_manual(values = c(COLOR_DISCOVERY, COLOR_TISSUE, COLOR_PLASMA)) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(title = "D: All Cohorts Overlay",
       x = "False Positive Rate", y = "True Positive Rate", color = "") +
  my_theme +
  theme(legend.position = c(0.65, 0.2),
        legend.background = element_rect(fill = alpha("white", 0.9), color = "grey70"))

# Combine 4 panels
p_4panel <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(title = paste(length(top_feats), "piRNA Diagnostic Signature"),
                  theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 16)))

ggsave("results/validation/Figure_ROC_tissue_vs_plasma.png", p_4panel,
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave("results/validation/Figure_ROC_tissue_vs_plasma.pdf", p_4panel,
       width = 10, height = 8)
cat("  Saved: results/validation/Figure_ROC_tissue_vs_plasma.png\n")
cat("  Saved: results/validation/Figure_ROC_tissue_vs_plasma.pdf\n")

# ==============================================================================
# 4. PER-COHORT INDIVIDUAL PLOTS
# ==============================================================================
cat("\n========== 4. Per-Cohort Individual Plots ==========\n")

for (ev in list(eval_v1, eval_v2)) {
  cohort <- ev$cohort

  # ROC with CI band
  p_roc_ind <- plot_roc(
    ev$roc_obj,
    paste("ROC —", cohort),
    sprintf("AUC = %.3f (%.3f–%.3f)", ev$auc, ev$auc_lower, ev$auc_upper),
    ifelse(cohort == "yyfbatch1", COLOR_TISSUE, COLOR_PLASMA)
  )
  ggsave(sprintf("results/validation/%s_roc.png", cohort),
         p_roc_ind, width = 7, height = 6, dpi = 300, bg = "white")

  # Precision-Recall curve
  if (!is.null(ev$pr_obj)) {
    pr_df <- data.frame(recall = ev$pr_obj$curve[, 1], precision = ev$pr_obj$curve[, 2])
    prevalence <- mean(ev$y_bin)
    p_prc <- ggplot(pr_df, aes(recall, precision)) +
      geom_hline(yintercept = prevalence, linetype = "dashed", color = "grey50") +
      geom_path(color = ifelse(cohort == "yyfbatch1", COLOR_TISSUE, COLOR_PLASMA),
                linewidth = 1.2) +
      annotate("label", x = 0.95, y = 0.05,
               label = sprintf("AUPRC = %.3f", ev$pr_obj$auc.integral),
               hjust = 1, vjust = 0, size = 4, fill = "white") +
      scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
      labs(title = paste("Precision-Recall —", cohort),
           x = "Recall", y = "Precision") +
      my_theme
    ggsave(sprintf("results/validation/%s_prc.png", cohort),
           p_prc, width = 7, height = 6, dpi = 300, bg = "white")
  }

  # Prediction probability histogram
  hist_df <- data.frame(prob = ev$prob, Group = ev$y_true)
  hist_df$uncertain <- ifelse(hist_df$prob >= 0.4 & hist_df$prob <= 0.6,
                              "Uncertain", "Confident")
  p_hist <- ggplot(hist_df, aes(x = prob, fill = Group)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
    geom_vline(xintercept = optimal_thr, linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(0.4, 0.6), linetype = "dotted", color = "grey50") +
    scale_fill_manual(values = c("Normal" = "#4393C3", "Tumor" = "#D6604D")) +
    labs(title = paste("Prediction Probability —", cohort),
         subtitle = sprintf("Uncertain zone [0.4–0.6]: %d samples (%.1f%%)",
                            sum(hist_df$uncertain == "Uncertain"),
                            100 * mean(hist_df$uncertain == "Uncertain")),
         x = "Predicted P(Tumor)", y = "Count") +
    my_theme
  ggsave(sprintf("results/validation/%s_probability_hist.png", cohort),
         p_hist, width = 7, height = 6, dpi = 300, bg = "white")

  cat(sprintf("  Saved plots for %s\n", cohort))
}

# ==============================================================================
# 5. PER-COHORT-PER-DATASET ROC (all 7 cohorts)
# ==============================================================================
cat("\n========== 5. All 7 Cohorts ROC Overlay ==========\n")

all_batches <- unique(combat_df_all$Batch)
roc_all7 <- list()

for (b in all_batches) {
  sub <- combat_df_all[combat_df_all$Batch == b, ]
  if (nrow(sub) < 5) next
  prob <- predict(model, sub[, top_feats, drop = FALSE], type = "prob")$Tumor
  y01 <- ifelse(sub$Group == "Tumor", 1, 0)
  if (length(unique(y01)) < 2) next
  roc_obj <- roc(y01, prob, direction = "<", quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  label <- sprintf("%s AUC=%.3f", b, auc_val)
  roc_all7[[b]] <- make_roc_df(roc_obj, label)
}
roc_all7_df <- do.call(rbind, roc_all7)
roc_all7_df$Set <- factor(roc_all7_df$Set, levels = unique(roc_all7_df$Set))

# Colors: tissue cohorts blue/red shades, plasma green
n_batches <- length(unique(roc_all7_df$Set))
cohort_colors <- COLOR_PALETTE[seq_len(n_batches)]

p_all7 <- ggplot(roc_all7_df, aes(fpr, tpr, color = Set)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_path(linewidth = 0.9) +
  scale_color_manual(values = cohort_colors) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
  labs(title = "ROC — All Cohorts (Individual Prediction)",
       subtitle = paste(length(top_feats), "piRNA signature"),
       x = "False Positive Rate", y = "True Positive Rate", color = "") +
  my_theme +
  theme(legend.position = "right", legend.text = element_text(size = 9))
ggsave("results/validation/roc_all_7_cohorts.png", p_all7,
       width = 9, height = 7, dpi = 300, bg = "white")
cat("  Saved: results/validation/roc_all_7_cohorts.png\n")

# ==============================================================================
# 6. ROBUSTNESS CHECKS
# ==============================================================================
cat("\n========== 6. Robustness Checks ==========\n")

for (ev in list(eval_v1, eval_v2)) {
  cohort <- ev$cohort
  cat(sprintf("\n--- %s ---\n", cohort))

  # 6a. Bootstrap (1000 resamples)
  cat("  Bootstrap AUC (1000 resamples)...\n")
  set.seed(SEED)
  boot_aucs <- numeric(1000)
  for (i in seq_len(1000)) {
    idx <- sample(length(ev$y_bin), replace = TRUE)
    yt <- ev$y_bin[idx]; ys <- ev$prob[idx]
    if (length(unique(yt)) < 2) { boot_aucs[i] <- NA; next }
    boot_aucs[i] <- as.numeric(auc(roc(yt, ys, direction = "<", quiet = TRUE)))
  }
  boot_aucs <- na.omit(boot_aucs)
  boot_ci <- quantile(boot_aucs, c(0.025, 0.975))
  cat(sprintf("    Bootstrap median AUC: %.4f, 95%% CI: [%.4f, %.4f]\n",
              median(boot_aucs), boot_ci[1], boot_ci[2]))

  p_boot <- ggplot(data.frame(AUC = boot_aucs), aes(AUC)) +
    geom_histogram(bins = 50, fill = "#1F78B4", alpha = 0.7) +
    geom_vline(xintercept = ev$auc, color = "#E31A1C", linewidth = 1) +
    geom_vline(xintercept = boot_ci, linetype = "dashed", color = "#E31A1C") +
    labs(title = paste("Bootstrap AUC —", cohort),
         subtitle = sprintf("Median=%.3f, 95%% CI=[%.3f, %.3f]",
                            median(boot_aucs), boot_ci[1], boot_ci[2]),
         x = "AUC", y = "Count") +
    my_theme
  ggsave(sprintf("results/validation/%s_bootstrap.png", cohort),
         p_boot, width = 7, height = 5, dpi = 300, bg = "white")

  # 6b. Permutation test (1000 label shuffles)
  cat("  Permutation test (1000 shuffles)...\n")
  set.seed(SEED)
  null_aucs <- numeric(1000)
  for (i in seq_len(1000)) {
    perm_y <- sample(ev$y_bin)
    if (length(unique(perm_y)) < 2) { null_aucs[i] <- 0.5; next }
    null_aucs[i] <- as.numeric(auc(roc(perm_y, ev$prob, direction = "<", quiet = TRUE)))
  }
  perm_p <- mean(null_aucs >= ev$auc)
  cat(sprintf("    Permutation p-value: %.4f\n", perm_p))

  p_perm <- ggplot(data.frame(AUC = null_aucs), aes(AUC)) +
    geom_histogram(bins = 50, fill = "grey70", alpha = 0.7) +
    geom_vline(xintercept = ev$auc, color = "#E31A1C", linewidth = 1.2) +
    labs(title = paste("Permutation Test —", cohort),
         subtitle = sprintf("Observed AUC=%.3f, p=%.4f", ev$auc, perm_p),
         x = "Null AUC Distribution", y = "Count") +
    my_theme
  ggsave(sprintf("results/validation/%s_permutation.png", cohort),
         p_perm, width = 7, height = 5, dpi = 300, bg = "white")

  cat(sprintf("  Saved bootstrap + permutation plots for %s\n", cohort))
}

# ==============================================================================
# 7. PERFORMANCE DROP ASSESSMENT
# ==============================================================================
cat("\n========== 7. Performance Drop Assessment ==========\n")

delta_tissue <- auc_cv - eval_v1$auc
delta_plasma <- auc_cv - eval_v2$auc

cat(sprintf("  Delta(train–tissue): %.4f\n", delta_tissue))
cat(sprintf("  Delta(train–plasma): %.4f\n", delta_plasma))

if (delta_tissue > 0.10) cat("  SEVERE WARNING: Tissue AUC drop > 0.10!\n")
else if (delta_tissue > 0.05) cat("  WARNING: Tissue AUC drop > 0.05 — potential overfitting.\n")

if (delta_plasma > 0.10) cat("  SEVERE WARNING: Plasma AUC drop > 0.10!\n")
else if (delta_plasma > 0.05) cat("  WARNING: Plasma AUC drop > 0.05 — potential overfitting.\n")

if (eval_v2$auc < eval_v1$auc) {
  cat("  NOTE: Plasma AUC < Tissue AUC (expected due to lower circulating piRNA abundance).\n")
  cat(sprintf("  Tissue–Plasma AUC gap: %.4f\n", eval_v1$auc - eval_v2$auc))
}

# AUC bar chart
auc_bar_df <- data.frame(
  Cohort = c("Discovery CV", "yyfbatch1\n(Tissue)", "yyfbatch2\n(Plasma)"),
  AUC    = c(auc_cv, eval_v1$auc, eval_v2$auc),
  Lower  = c(auc_cv_ci[1], eval_v1$auc_lower, eval_v2$auc_lower),
  Upper  = c(auc_cv_ci[3], eval_v1$auc_upper, eval_v2$auc_upper),
  Type   = c("Discovery", "Tissue", "Plasma")
)
auc_bar_df$Cohort <- factor(auc_bar_df$Cohort, levels = auc_bar_df$Cohort)

p_bar <- ggplot(auc_bar_df, aes(x = Cohort, y = AUC, fill = Type)) +
  geom_col(alpha = 0.8, width = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.5, size = 4.5) +
  scale_fill_manual(values = c("Discovery" = COLOR_DISCOVERY,
                                "Tissue" = COLOR_TISSUE,
                                "Plasma" = COLOR_PLASMA)) +
  ylim(c(0, 1.1)) +
  labs(title = "AUC Comparison Across Cohorts", x = "", y = "AUC") +
  my_theme + theme(legend.position = "none")
ggsave("results/validation/auc_comparison_bar.png", p_bar,
       width = 7, height = 6, dpi = 300, bg = "white")
cat("  Saved: results/validation/auc_comparison_bar.png\n")

# ==============================================================================
# 8. RESULTS TABLE
# ==============================================================================
cat("\n========== 8. Summary Results Table ==========\n")

results_table <- data.frame(
  Cohort  = c("Discovery CV", "yyfbatch1", "yyfbatch2"),
  Type    = c("Tissue", "Tissue", "Plasma/Serum"),
  N       = c(nrow(pred_cv), eval_v1$n, eval_v2$n),
  AUC     = c(auc_cv, eval_v1$auc, eval_v2$auc),
  AUC_CI  = c(sprintf("%.3f–%.3f", auc_cv_ci[1], auc_cv_ci[3]),
              sprintf("%.3f–%.3f", eval_v1$auc_lower, eval_v1$auc_upper),
              sprintf("%.3f–%.3f", eval_v2$auc_lower, eval_v2$auc_upper)),
  Sens    = c(NA, eval_v1$sens, eval_v2$sens),
  Spec    = c(NA, eval_v1$spec, eval_v2$spec),
  PPV     = c(NA, eval_v1$ppv, eval_v2$ppv),
  NPV     = c(NA, eval_v1$npv, eval_v2$npv),
  stringsAsFactors = FALSE
)

cat("\n")
print(results_table)

write.csv(results_table, "results/validation/validation_results.csv", row.names = FALSE)
cat("\n  Saved: results/validation/validation_results.csv\n")

# Save key objects for Phase 4
saveRDS(list(eval_v1 = eval_v1, eval_v2 = eval_v2,
             auc_cv = auc_cv, auc_cv_ci = auc_cv_ci,
             roc_cv = roc_cv, optimal_thr = optimal_thr),
        "results/validation/validation_details.rds")

end_time <- Sys.time()
cat(sprintf("\n=== PHASE 3 complete. Runtime: %.1f minutes ===\n",
            as.numeric(difftime(end_time, start_time, units = "mins"))))
cat("=================================================================\n")
cat("  Next: run 05_clinical_interpretation_report.R\n")
cat("=================================================================\n")

################################################################################
#                                                                              #
#   02 — Feature Selection → ~10 piRNA Consensus Set (PHASE 1)                #
#                                                                              #
#   8 feature selection methods on train_df ONLY.                             #
#   Consensus construction (frequency threshold + forward stepwise).          #
#   Forward/backward/swap refinement.                                         #
#   Final signature: 7–12 piRNAs.                                             #
#                                                                              #
#   Requires: results from 01_data_loading_batch_correction.R                 #
#   Produces: results/feature_selection/* , results/models/final_features.rds  #
#                                                                              #
################################################################################

start_time <- Sys.time()
cat("=================================================================\n")
cat("  PHASE 1: Feature Selection → ~10 piRNA Consensus Set\n")
cat("=================================================================\n\n")

# ==============================================================================
# 0. PACKAGES
# ==============================================================================
required_pkgs <- c(
  "limma", "caret", "randomForest", "glmnet", "Boruta", "praznik",
  "pROC", "ggplot2", "dplyr", "UpSetR", "pheatmap", "viridis",
  "xgboost"
)
bioc_pkgs <- c("limma")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}
for (pkg in setdiff(required_pkgs, bioc_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(limma)
  library(caret)
  library(randomForest)
  library(glmnet)
  library(Boruta)
  library(praznik)
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(UpSetR)
  library(pheatmap)
  library(viridis)
})
cat("All packages loaded.\n\n")

SEED <- 2024
set.seed(SEED)

my_theme <- theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        legend.position = "bottom")

# ==============================================================================
# 0.1 LOAD PHASE 0 OUTPUTS
# ==============================================================================
cat("Loading Phase 0 outputs...\n")
train_df   <- readRDS("results/models/train_df.rds")
gene_cols  <- readRDS("results/models/common_pirnas.rds")

x_train <- as.matrix(train_df[, gene_cols])
y_train <- train_df$Group

cat("  Training samples:", nrow(train_df), "\n")
cat("  piRNA features:", length(gene_cols), "\n")
cat("  Class distribution:", table(y_train)["Normal"], "Normal,",
    table(y_train)["Tumor"], "Tumor\n\n")

fs_results <- list()

# ==============================================================================
# 1. limma (empirical Bayes DE)
# ==============================================================================
cat("--- Method 1/8: limma ---\n")
tryCatch({
  design <- model.matrix(~ y_train)
  fit <- lmFit(t(x_train), design)
  fit <- eBayes(fit)
  limma_res <- topTable(fit, coef = 2, number = Inf, sort.by = "p")

  fs_limma <- rownames(limma_res)[limma_res$adj.P.Val < 0.01 & abs(limma_res$logFC) > 1.0]
  fs_results$limma <- fs_limma
  cat("  limma:", length(fs_limma), "piRNAs (adj.p<0.01, |log2FC|>1.0)\n")

  # Volcano plot
  limma_res$sig <- ifelse(limma_res$adj.P.Val < 0.01 & abs(limma_res$logFC) > 1.0,
                           ifelse(limma_res$logFC > 0, "Up", "Down"), "NS")
  p_vol <- ggplot(limma_res, aes(x = logFC, y = -log10(adj.P.Val), color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Down" = "#1F78B4", "NS" = "grey70", "Up" = "#E31A1C")) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", color = "grey40") +
    labs(title = "Volcano Plot — limma DE",
         x = "log2 Fold Change", y = "-log10(adj. P-value)", color = "") +
    my_theme
  ggsave("results/feature_selection/volcano_plot.png", p_vol,
         width = 8, height = 6, dpi = 300, bg = "white")
  cat("  Saved: results/feature_selection/volcano_plot.png\n")
}, error = function(e) {
  cat("  limma FAILED:", conditionMessage(e), "\n")
  fs_results$limma <<- character(0)
})

# ==============================================================================
# 2. Wilcoxon rank-sum
# ==============================================================================
cat("--- Method 2/8: Wilcoxon rank-sum ---\n")
tryCatch({
  pvals <- apply(x_train, 2, function(col) {
    suppressWarnings(wilcox.test(col ~ y_train)$p.value)
  })
  padj <- p.adjust(pvals, method = "BH")

  # Median fold change
  med_tumor  <- apply(x_train[y_train == "Tumor", ], 2, median)
  med_normal <- apply(x_train[y_train == "Normal", ], 2, median)
  med_fc     <- abs(med_tumor - med_normal)

  fs_wilcoxon <- names(padj)[padj < 0.01 & med_fc > 1.5]
  fs_results$wilcoxon <- fs_wilcoxon
  cat("  Wilcoxon:", length(fs_wilcoxon), "piRNAs (adj.p<0.01, |medFC|>1.5)\n")
}, error = function(e) {
  cat("  Wilcoxon FAILED:", conditionMessage(e), "\n")
  fs_results$wilcoxon <<- character(0)
})

# ==============================================================================
# 3. Random Forest importance
# ==============================================================================
cat("--- Method 3/8: RF importance ---\n")
tryCatch({
  set.seed(SEED)
  rf_fs <- randomForest(x = x_train, y = y_train, ntree = 1000, importance = TRUE)
  rf_imp <- importance(rf_fs, type = 1)  # MeanDecreaseAccuracy
  rf_gini <- importance(rf_fs, type = 2)  # MeanDecreaseGini
  rf_sorted <- sort(rf_imp[, 1], decreasing = TRUE)

  threshold <- mean(rf_imp[, 1]) + sd(rf_imp[, 1])
  fs_rf <- names(rf_sorted[rf_sorted > threshold])
  if (length(fs_rf) > 50) fs_rf <- names(rf_sorted[1:50])

  fs_results$rf_importance <- fs_rf
  cat("  RF importance:", length(fs_rf), "piRNAs (mean+1SD threshold, cap 50)\n")

  # Importance barplot (top 30)
  top30 <- head(data.frame(piRNA = names(rf_sorted), Importance = rf_sorted), 30)
  top30$piRNA <- factor(top30$piRNA, levels = rev(top30$piRNA))
  p_rf <- ggplot(top30, aes(x = piRNA, y = Importance)) +
    geom_bar(stat = "identity", fill = "#1F78B4") +
    coord_flip() +
    labs(title = "RF Feature Importance (Top 30)", x = "", y = "MeanDecreaseAccuracy") +
    my_theme
  ggsave("results/feature_selection/rf_importance.png", p_rf,
         width = 8, height = 8, dpi = 300, bg = "white")
  cat("  Saved: results/feature_selection/rf_importance.png\n")
}, error = function(e) {
  cat("  RF importance FAILED:", conditionMessage(e), "\n")
  fs_results$rf_importance <<- character(0)
})

# ==============================================================================
# 4. LASSO (glmnet, alpha=1)
# ==============================================================================
cat("--- Method 4/8: LASSO ---\n")
tryCatch({
  set.seed(SEED)
  y_bin <- as.numeric(y_train) - 1
  cv_lasso <- cv.glmnet(x_train, y_bin, family = "binomial", alpha = 1,
                         nfolds = 10, type.measure = "auc")
  lasso_coef <- coef(cv_lasso, s = "lambda.1se")
  fs_lasso <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
  fs_lasso <- setdiff(fs_lasso, "(Intercept)")
  fs_results$lasso <- fs_lasso
  cat("  LASSO:", length(fs_lasso), "piRNAs (non-zero coefs at lambda.1se)\n")
}, error = function(e) {
  cat("  LASSO FAILED:", conditionMessage(e), "\n")
  fs_results$lasso <<- character(0)
})

# ==============================================================================
# 5. Elastic Net (best alpha by CV)
# ==============================================================================
cat("--- Method 5/8: Elastic Net ---\n")
tryCatch({
  set.seed(SEED)
  if (!exists("y_bin")) y_bin <- as.numeric(y_train) - 1
  best_auc <- 0; best_alpha <- 0.5; best_en <- NULL
  for (a in seq(0.1, 0.9, by = 0.1)) {
    cv_en <- cv.glmnet(x_train, y_bin, family = "binomial", alpha = a,
                        nfolds = 10, type.measure = "auc")
    if (max(cv_en$cvm) > best_auc) {
      best_auc <- max(cv_en$cvm); best_alpha <- a; best_en <- cv_en
    }
  }
  en_coef <- coef(best_en, s = "lambda.1se")
  fs_elasticnet <- rownames(en_coef)[en_coef[, 1] != 0]
  fs_elasticnet <- setdiff(fs_elasticnet, "(Intercept)")
  fs_results$elasticnet <- fs_elasticnet
  cat("  Elastic Net (best alpha=", best_alpha, "):", length(fs_elasticnet), "piRNAs\n")
}, error = function(e) {
  cat("  Elastic Net FAILED:", conditionMessage(e), "\n")
  fs_results$elasticnet <<- character(0)
})

# ==============================================================================
# 6. RFE (caret::rfe, rfFuncs)
# ==============================================================================
cat("--- Method 6/8: RFE ---\n")
tryCatch({
  set.seed(SEED)
  rfe_ctrl <- rfeControl(functions = rfFuncs, method = "cv", number = 10, repeats = 3)
  rfe_sizes <- c(5, 10, 15, 20, 30, 50)
  rfe_sizes <- rfe_sizes[rfe_sizes <= ncol(x_train)]
  rfe_result <- rfe(x = data.frame(x_train), y = y_train,
                    sizes = rfe_sizes, rfeControl = rfe_ctrl)
  fs_rfe <- predictors(rfe_result)
  fs_results$rfe <- fs_rfe
  cat("  RFE:", length(fs_rfe), "piRNAs (optimal subset)\n")
}, error = function(e) {
  cat("  RFE FAILED:", conditionMessage(e), "\n")
  fs_results$rfe <<- character(0)
})

# ==============================================================================
# 7. Boruta
# ==============================================================================
cat("--- Method 7/8: Boruta ---\n")
tryCatch({
  set.seed(SEED)
  boruta_res <- Boruta(x = data.frame(x_train), y = y_train,
                       maxRuns = 300, doTrace = 0)
  fs_boruta <- getSelectedAttributes(boruta_res, withTentative = FALSE)
  fs_results$boruta <- fs_boruta
  cat("  Boruta:", length(fs_boruta), "Confirmed piRNAs\n")
}, error = function(e) {
  cat("  Boruta FAILED:", conditionMessage(e), "\n")
  fs_results$boruta <<- character(0)
})

# ==============================================================================
# 8. mRMR (praznik package)
# ==============================================================================
cat("--- Method 8/8: mRMR ---\n")
tryCatch({
  k_max <- min(50, ncol(x_train))
  mrmr_res <- MRMR(data.frame(x_train), y_train, k = k_max)
  fs_mrmr <- colnames(x_train)[mrmr_res$selection]
  fs_results$mrmr <- fs_mrmr
  cat("  mRMR:", length(fs_mrmr), "piRNAs (top 50)\n")
}, error = function(e) {
  cat("  mRMR FAILED:", conditionMessage(e), "\n")
  fs_results$mrmr <<- character(0)
})

# ==============================================================================
# 9. CONSENSUS CONSTRUCTION
# ==============================================================================
cat("\n========== STEP 9: Consensus Construction ==========\n")

active_methods <- names(fs_results)[sapply(fs_results, length) > 0]
n_methods <- length(active_methods)
cat("Active methods:", n_methods, "/8\n")

all_selected <- unique(unlist(fs_results[active_methods]))
freq_table <- data.frame(
  piRNA = all_selected,
  count = sapply(all_selected, function(m) {
    sum(sapply(fs_results[active_methods], function(fs) m %in% fs))
  }),
  stringsAsFactors = FALSE
)
freq_table <- freq_table[order(-freq_table$count), ]
rownames(freq_table) <- NULL
write.csv(freq_table, "results/feature_selection/fs_frequency_table.csv", row.names = FALSE)

cat("\nTop 20 piRNAs by selection frequency:\n")
print(head(freq_table, 20))

# --- UpSet plot ---
cat("\nGenerating UpSet plot...\n")
tryCatch({
  upset_list <- lapply(fs_results[active_methods], function(x) x)
  # Convert to binary matrix for UpSetR
  all_feats <- unique(unlist(upset_list))
  upset_mat <- data.frame(
    sapply(upset_list, function(s) as.integer(all_feats %in% s))
  )
  rownames(upset_mat) <- all_feats

  png("results/feature_selection/upset_plot.png", width = 10, height = 7,
      units = "in", res = 300)
  print(upset(upset_mat, sets = names(upset_list),
              order.by = "freq", nsets = n_methods,
              mainbar.y.label = "Intersection Size",
              sets.x.label = "Features per Method"))
  dev.off()
  cat("  Saved: results/feature_selection/upset_plot.png\n")
}, error = function(e) cat("  UpSet plot failed:", conditionMessage(e), "\n"))

# --- STRATEGY A: Frequency threshold ---
cat("\n--- Strategy A: Frequency Threshold ---\n")
for (thr in c(6, 5, 4, 3)) {
  strat_a <- freq_table$piRNA[freq_table$count >= thr]
  if (length(strat_a) >= 7 && length(strat_a) <= 12) {
    cat("  Threshold >=", thr, "/", n_methods, "=>", length(strat_a), "piRNAs [OK]\n")
    break
  } else {
    cat("  Threshold >=", thr, "/", n_methods, "=>", length(strat_a), "piRNAs")
    if (length(strat_a) > 12) {
      cat(" [too many, tightening]\n")
      next
    }
    if (length(strat_a) < 7) {
      cat(" [too few, relaxing]\n")
    }
  }
}
# Force into 7-12 range
if (length(strat_a) > 12) strat_a <- strat_a[1:12]
if (length(strat_a) < 7) strat_a <- freq_table$piRNA[1:min(10, nrow(freq_table))]
cat("  Strategy A set:", length(strat_a), "piRNAs\n")

# --- STRATEGY B: Forward stepwise on frequency-ranked list ---
cat("\n--- Strategy B: Forward Stepwise ---\n")
ranked_pirnas <- freq_table$piRNA

# Add RF importance as tiebreaker
if (exists("rf_imp")) {
  freq_table$rf_importance <- rf_imp[match(freq_table$piRNA, rownames(rf_imp)), 1]
  freq_table$rf_importance[is.na(freq_table$rf_importance)] <- 0
  freq_table <- freq_table[order(-freq_table$count, -freq_table$rf_importance), ]
  ranked_pirnas <- freq_table$piRNA
}

# Forward stepwise with 5-fold CV AUC
strat_b <- c()
best_auc_b <- 0
no_improve <- 0

set.seed(SEED)
folds <- createFolds(y_train, k = 5, returnTrain = TRUE)

quick_cv_auc <- function(features, x, y, folds, seed = SEED) {
  if (length(features) == 0) return(0.5)
  aucs <- numeric(length(folds))
  for (i in seq_along(folds)) {
    tr_idx <- folds[[i]]
    te_idx <- setdiff(seq_along(y), tr_idx)
    set.seed(seed)
    rf_tmp <- randomForest(x = x[tr_idx, features, drop = FALSE], y = y[tr_idx],
                           ntree = 500, importance = FALSE)
    prob <- predict(rf_tmp, x[te_idx, features, drop = FALSE], type = "prob")[, "Tumor"]
    y01 <- ifelse(y[te_idx] == "Tumor", 1, 0)
    if (length(unique(y01)) < 2) { aucs[i] <- 0.5; next }
    aucs[i] <- as.numeric(auc(roc(y01, prob, direction = "<", quiet = TRUE)))
  }
  mean(aucs)
}

for (i in seq_along(ranked_pirnas)) {
  if (length(strat_b) >= 12) break

  candidate <- c(strat_b, ranked_pirnas[i])
  cv_auc <- quick_cv_auc(candidate, x_train, y_train, folds)

  if (cv_auc > best_auc_b + 0.005 || length(strat_b) < 7) {
    if (cv_auc > best_auc_b) {
      best_auc_b <- cv_auc
      no_improve <- 0
    } else {
      no_improve <- no_improve + 1
    }
    strat_b <- candidate
    cat(sprintf("  + %s => %d feats, CV AUC=%.4f\n",
                ranked_pirnas[i], length(candidate), cv_auc))
  } else {
    no_improve <- no_improve + 1
  }

  # Stop if AUC improvement < 0.005 for 3 consecutive additions AND we have >= 7
  if (no_improve >= 3 && length(strat_b) >= 7) break
}
cat("  Strategy B set:", length(strat_b), "piRNAs, CV AUC:", round(best_auc_b, 4), "\n")

# --- Pick best strategy ---
cat("\n--- Picking best strategy ---\n")
auc_a <- quick_cv_auc(strat_a, x_train, y_train, folds)
auc_b <- best_auc_b
cat("  Strategy A:", length(strat_a), "piRNAs, CV AUC:", round(auc_a, 4), "\n")
cat("  Strategy B:", length(strat_b), "piRNAs, CV AUC:", round(auc_b, 4), "\n")

if (abs(auc_a - auc_b) < 0.01) {
  consensus_set <- if (length(strat_a) <= length(strat_b)) strat_a else strat_b
  cat("  Similar AUC — picking smaller set.\n")
} else if (auc_a > auc_b) {
  consensus_set <- strat_a
  cat("  Picking Strategy A (higher AUC).\n")
} else {
  consensus_set <- strat_b
  cat("  Picking Strategy B (higher AUC).\n")
}

# ==============================================================================
# 10. FORWARD/BACKWARD/SWAP REFINEMENT
# ==============================================================================
cat("\n========== STEP 10: Forward/Backward/Swap Refinement ==========\n")

# Backward: remove piRNAs that don't decrease AUC by >= 0.003
cat("Backward pruning...\n")
current_set <- consensus_set
current_auc <- quick_cv_auc(current_set, x_train, y_train, folds)
improved <- TRUE

while (improved && length(current_set) > 7) {
  improved <- FALSE
  for (j in seq_along(current_set)) {
    candidate <- current_set[-j]
    cand_auc <- quick_cv_auc(candidate, x_train, y_train, folds)
    if (current_auc - cand_auc < 0.003) {
      cat(sprintf("  - Removed %s => %d feats, AUC %.4f -> %.4f\n",
                  current_set[j], length(candidate), current_auc, cand_auc))
      current_set <- candidate
      current_auc <- cand_auc
      improved <- TRUE
      break
    }
  }
}

# Swap: try replacing each piRNA with the next-best not in the set
cat("Swap search...\n")
untested <- setdiff(ranked_pirnas[1:min(30, length(ranked_pirnas))], current_set)
swap_improved <- TRUE
while (swap_improved) {
  swap_improved <- FALSE
  for (j in seq_along(current_set)) {
    for (new_feat in untested) {
      candidate <- current_set
      candidate[j] <- new_feat
      cand_auc <- quick_cv_auc(candidate, x_train, y_train, folds)
      if (cand_auc > current_auc + 0.005) {
        cat(sprintf("  SWAP %s -> %s => AUC %.4f -> %.4f\n",
                    current_set[j], new_feat, current_auc, cand_auc))
        old_feat <- current_set[j]
        current_set[j] <- new_feat
        untested <- setdiff(untested, new_feat)
        untested <- c(untested, old_feat)
        current_auc <- cand_auc
        swap_improved <- TRUE
        break
      }
    }
    if (swap_improved) break
  }
}

# Force into 7–12 range
if (length(current_set) > 12) current_set <- current_set[1:12]
if (length(current_set) < 7) {
  needed <- 7 - length(current_set)
  extra <- setdiff(ranked_pirnas, current_set)[1:needed]
  current_set <- c(current_set, extra)
}

final_features <- current_set
final_auc <- quick_cv_auc(final_features, x_train, y_train, folds)

cat("\n*** FINAL piRNA SIGNATURE:", length(final_features), "piRNAs ***\n")
cat(paste(final_features, collapse = ", "), "\n")
cat("CV AUC:", round(final_auc, 4), "\n")

# ==============================================================================
# 11. JUSTIFICATION TABLE
# ==============================================================================
cat("\n========== STEP 11: Justification Table ==========\n")

justify <- data.frame(piRNA = final_features, stringsAsFactors = FALSE)

# N_methods
justify$N_methods <- sapply(justify$piRNA, function(p) {
  sum(sapply(fs_results[active_methods], function(fs) p %in% fs))
})

# which_methods
justify$which_methods <- sapply(justify$piRNA, function(p) {
  paste(names(fs_results[active_methods])[sapply(fs_results[active_methods], function(fs) p %in% fs)],
        collapse = ", ")
})

# limma stats
if (exists("limma_res")) {
  justify$limma_adjP   <- limma_res[justify$piRNA, "adj.P.Val"]
  justify$limma_log2FC <- limma_res[justify$piRNA, "logFC"]
} else {
  justify$limma_adjP   <- NA
  justify$limma_log2FC <- NA
}

# Individual AUC (univariate, 5-fold CV)
justify$individual_AUC <- sapply(justify$piRNA, function(p) {
  quick_cv_auc(p, x_train, y_train, folds)
})

# RF importance rank
if (exists("rf_imp")) {
  rf_rank <- rank(-rf_imp[, 1])
  justify$RF_importance_rank <- rf_rank[match(justify$piRNA, names(rf_rank))]
} else {
  justify$RF_importance_rank <- NA
}

# LASSO coefficient
if (exists("cv_lasso")) {
  lc <- coef(cv_lasso, s = "lambda.1se")
  lasso_vals <- as.numeric(lc[match(justify$piRNA, rownames(lc)), 1])
  lasso_vals[is.na(lasso_vals)] <- 0
  justify$LASSO_coef <- lasso_vals
} else {
  justify$LASSO_coef <- NA
}

# mRMR rank
if (exists("mrmr_res")) {
  mrmr_feats <- colnames(x_train)[mrmr_res$selection]
  justify$mRMR_rank <- match(justify$piRNA, mrmr_feats)
} else {
  justify$mRMR_rank <- NA
}

cat("\nFull Justification Table:\n")
print(justify)

# Check individual AUCs
low_auc <- justify$piRNA[justify$individual_AUC < 0.55]
if (length(low_auc) > 0) {
  cat("\nWARNING: piRNAs with individual AUC < 0.55:", paste(low_auc, collapse = ", "), "\n")
  for (p in low_auc) {
    n_sel <- justify$N_methods[justify$piRNA == p]
    if (n_sel >= 5) cat("  ", p, "kept (N_methods >= 5)\n")
    else cat("  ", p, "FLAGGED for review (N_methods =", n_sel, ")\n")
  }
}

write.csv(justify, "results/feature_selection/fs_summary_table.csv", row.names = FALSE)
cat("  Saved: results/feature_selection/fs_summary_table.csv\n")

# ==============================================================================
# 12. INCREMENTAL AUC CURVE
# ==============================================================================
cat("\n========== STEP 12: Incremental AUC Curve ==========\n")

incr_aucs <- numeric(length(final_features))
for (k in seq_along(final_features)) {
  incr_aucs[k] <- quick_cv_auc(final_features[1:k], x_train, y_train, folds)
}

incr_df <- data.frame(n_piRNAs = seq_along(final_features), CV_AUC = incr_aucs)
p_incr <- ggplot(incr_df, aes(x = n_piRNAs, y = CV_AUC)) +
  geom_line(linewidth = 1.2, color = "#E31A1C") +
  geom_point(size = 3, color = "#E31A1C") +
  geom_text(aes(label = sprintf("%.3f", CV_AUC)), vjust = -1, size = 3.5) +
  scale_x_continuous(breaks = seq_along(final_features)) +
  ylim(c(min(incr_aucs) - 0.05, 1)) +
  labs(title = "Incremental AUC as piRNAs Are Added",
       x = "Number of piRNAs", y = "5-Fold CV AUC") +
  my_theme
ggsave("results/feature_selection/incremental_auc.png", p_incr,
       width = 8, height = 6, dpi = 300, bg = "white")
cat("  Saved: results/feature_selection/incremental_auc.png\n")

# ==============================================================================
# 13. PAIRWISE CORRELATION HEATMAP
# ==============================================================================
cat("\n========== STEP 13: Correlation Heatmap ==========\n")

cor_mat <- cor(x_train[, final_features])
png("results/feature_selection/signature_correlation.png",
    width = 8, height = 7, units = "in", res = 300)
pheatmap(cor_mat,
         display_numbers = TRUE,
         number_format = "%.2f",
         color = viridis(50),
         main = "Pairwise Correlation — Final Signature",
         fontsize = 10)
dev.off()
cat("  Saved: results/feature_selection/signature_correlation.png\n")

# ==============================================================================
# 14. SAVE FINAL OUTPUTS
# ==============================================================================
cat("\n========== STEP 14: Saving Outputs ==========\n")

writeLines(final_features, "results/feature_selection/consensus_pirnas.txt")
saveRDS(final_features, "results/models/final_features.rds")
saveRDS(fs_results, "results/feature_selection/fs_results_all_methods.rds")

cat("  Saved: results/feature_selection/consensus_pirnas.txt\n")
cat("  Saved: results/models/final_features.rds\n")
cat("  Saved: results/feature_selection/fs_results_all_methods.rds\n")

cat(sprintf("\nFinal piRNA signature: %d piRNAs selected by at least %d/%d methods\n",
            length(final_features), min(justify$N_methods), n_methods))

end_time <- Sys.time()
cat(sprintf("\n=== PHASE 1 complete. Runtime: %.1f minutes ===\n",
            as.numeric(difftime(end_time, start_time, units = "mins"))))
cat("=================================================================\n")
cat("  Next: run 03_model_training.R\n")
cat("=================================================================\n")

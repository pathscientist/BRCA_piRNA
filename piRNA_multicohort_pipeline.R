################################################################################
#                                                                              #
#   Breast Cancer piRNA Multi-Cohort Diagnostic Pipeline                       #
#                                                                              #
#   Global ComBat → Feature Selection (<10) → RF Model →                       #
#   Independent Validation → Logistic Regression Forest Plot →                 #
#   Subgroup ROC & T-Score Boxplots                                            #
#                                                                              #
#   7 datasets: BRCA1, PRJNA294226, PRJNA482141, PRJNA808405,                  #
#               PRJNA934049, yyfbatch1, yyfbatch2                              #
#                                                                              #
################################################################################

start_time <- Sys.time()

# ==============================================================================
# 0. PACKAGES
# ==============================================================================
required_pkgs <- c(
  "sva", "caret", "randomForest", "glmnet", "pROC", "PRROC",
  "ggplot2", "dplyr", "tidyr", "gridExtra", "ggpubr",
  "forestplot", "survival", "doParallel", "foreach",
  "pheatmap", "Boruta", "praznik"
)

bioc_pkgs <- c("sva", "limma")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", quiet = TRUE)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

for (pkg in setdiff(required_pkgs, bioc_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE, quiet = TRUE)
}

suppressPackageStartupMessages({
  library(sva)
  library(limma)
  library(caret)
  library(randomForest)
  library(glmnet)
  library(pROC)
  library(PRROC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(ggpubr)
  library(survival)
  library(doParallel)
  library(pheatmap)
  library(Boruta)
  library(praznik)
})

cat("All packages loaded.\n")

# ==============================================================================
# 0.1 GLOBAL SETTINGS
# ==============================================================================
SEED <- 2024
set.seed(SEED)

COLOR_PALETTE <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00",
                    "#6A3D9A", "#B15928", "#A6CEE3", "#FB9A99")

my_theme <- theme_bw() +
  theme(
    panel.grid.major  = element_line(color = "grey90"),
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(size = 12, color = "black"),
    axis.title        = element_text(size = 14, face = "bold"),
    plot.title        = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle     = element_text(size = 11, hjust = 0.5, color = "grey40"),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1)
  )

dir.create("results/feature_selection", recursive = TRUE, showWarnings = FALSE)
dir.create("results/models",            recursive = TRUE, showWarnings = FALSE)
dir.create("results/validation",        recursive = TRUE, showWarnings = FALSE)
dir.create("results/forest_plot",       recursive = TRUE, showWarnings = FALSE)
dir.create("results/subgroup",          recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================
cat("\n========== PHASE 1: DATA LOADING ==========\n")

# >>> EDIT THIS: set your working directory and file paths <<<
# setwd("C:/HP/PhD/Projects/GEO BRCA/brca_dataset_1101tpm")

# Expression matrices (rows = piRNAs, cols = samples)
# BRCA1      <- read.table("BRCA1.txt",      header=TRUE, row.names=1, sep="\t")
# PRJNA294226<- read.table("PRJNA294226.txt", header=TRUE, row.names=1, sep="\t")
# PRJNA482141<- read.table("PRJNA482141.txt", header=TRUE, row.names=1, sep="\t")
# PRJNA808405<- read.table("PRJNA808405.txt", header=TRUE, row.names=1, sep="\t")
# PRJNA934049<- read.table("PRJNA934049.txt", header=TRUE, row.names=1, sep="\t")
# yyfbatch1  <- read.table("yyfbatch1.txt",   header=TRUE, row.names=1, sep="\t")
# yyfbatch2  <- read.table("yyfbatch2.txt",   header=TRUE, row.names=1, sep="\t")

# Phenotype files (rows = samples, must have "Group" column)
# BRCA1_pheno      <- read.table("BRCA1_pheno.txt",      header=TRUE, row.names=1, sep="\t")
# PRJNA294226_pheno<- read.table("PRJNA294226_pheno.txt", header=TRUE, row.names=1, sep="\t")
# ...

# --- Reformat helper: genes-as-rows → samples-as-rows, add Group column ---
reformat_dataset <- function(expr, pheno, group_col = "Group") {
  common <- intersect(colnames(expr), rownames(pheno))
  expr  <- expr[, common, drop = FALSE]
  pheno <- pheno[common, , drop = FALSE]
  df <- data.frame(Group = pheno[[group_col]], t(expr), check.names = FALSE)
  rownames(df) <- common
  df
}

recode_group <- function(df, group_col = "Group") {
  g <- as.character(df[[group_col]])
  g[g == "Benign"] <- "Normal"
  g[g == "Cancer"] <- "Tumor"
  df[[group_col]] <- factor(g, levels = c("Normal", "Tumor"))
  # Keep only Tumor and Normal
  df <- df[df[[group_col]] %in% c("Tumor", "Normal"), ]
  df[[group_col]] <- droplevels(df[[group_col]])
  df
}

# Build ready-list from your raw objects
# >>> UNCOMMENT and run after loading your data <<<
# expr_list <- list(BRCA1=BRCA1, PRJNA294226=PRJNA294226, PRJNA482141=PRJNA482141,
#                   PRJNA808405=PRJNA808405, PRJNA934049=PRJNA934049,
#                   yyfbatch1=yyfbatch1, yyfbatch2=yyfbatch2)
# pheno_list <- list(BRCA1=BRCA1_pheno, PRJNA294226=PRJNA294226_pheno,
#                    PRJNA482141=PRJNA482141_pheno, PRJNA808405=PRJNA808405_pheno,
#                    PRJNA934049=PRJNA934049_pheno, yyfbatch1=yyfbatch1_pheno,
#                    yyfbatch2=yyfbatch2_pheno)
#
# ready_list <- list()
# for (nm in names(expr_list)) {
#   ready_list[[nm]] <- recode_group(reformat_dataset(expr_list[[nm]], pheno_list[[nm]]))
# }

# --- OR load from processed CSV files ---
if (!exists("ready_list")) {
  if (dir.exists("processed_results")) {
    cat("Loading from processed_results/ directory...\n")
    files <- list.files("processed_results", pattern = "*.csv", full.names = TRUE)
    dataset_names <- gsub("processed_results/|_processed.csv", "", files)
    ready_list <- lapply(files, function(x) {
      read.csv(x, row.names = 1, stringsAsFactors = FALSE)
    })
    names(ready_list) <- dataset_names
    ready_list <- lapply(ready_list, recode_group)
  } else {
    stop("No data found. Please load your data into 'ready_list' first.\n",
         "See the commented-out code blocks above for instructions.")
  }
}

# Print dataset summary
cat("\nDataset summary:\n")
for (nm in names(ready_list)) {
  tb <- table(ready_list[[nm]]$Group)
  cat(sprintf("  %-15s: %d samples (Normal=%d, Tumor=%d)\n",
              nm, nrow(ready_list[[nm]]),
              tb["Normal"], tb["Tumor"]))
}


# ==============================================================================
# 2. BRCA1 BALANCING (Matched Pairs + 40% Remaining Tumors)
# ==============================================================================
cat("\n========== PHASE 2: BRCA1 BALANCING ==========\n")

balance_brca <- function(df, seed = 123) {
  set.seed(seed)
  idx_normal <- which(df$Group == "Normal")
  idx_tumor  <- which(df$Group == "Tumor")

  n_normal      <- length(idx_normal)
  n_tumor_total <- length(idx_tumor)

  # Step 1: matched pairs — N_normal tumors paired with N_normal normals
  n_matched <- min(n_normal, n_tumor_total)

  # Step 2: 40% of remaining tumors
  n_remaining <- n_tumor_total - n_matched
  n_extra     <- round(n_remaining * 0.40)

  # Total tumors to keep
  n_tumor_keep <- n_matched + n_extra
  idx_tumor_keep <- sample(idx_tumor, n_tumor_keep)

  df_balanced <- df[c(idx_normal, idx_tumor_keep), ]

  cat(sprintf("  BRCA1 balanced: Normal=%d, Tumor=%d (matched %d + 40%% of %d remaining = %d)\n",
              n_normal, n_tumor_keep, n_matched, n_remaining, n_extra))
  cat(sprintf("  Ratio Normal:Tumor = 1:%.1f\n", n_tumor_keep / n_normal))

  df_balanced
}

if ("BRCA1" %in% names(ready_list)) {
  ready_list[["BRCA1"]] <- balance_brca(ready_list[["BRCA1"]], seed = SEED)
} else {
  cat("  No BRCA1 dataset found; skipping balancing.\n")
}

# ==============================================================================
# 3. COMMON GENES + LOG2 + GLOBAL ComBat + Z-SCORE
# ==============================================================================
cat("\n========== PHASE 3: BATCH CORRECTION (Global ComBat) ==========\n")

# Find common piRNAs across all datasets
all_gene_sets <- lapply(ready_list, function(df) setdiff(colnames(df), "Group"))
common_genes  <- Reduce(intersect, all_gene_sets)
cat("Common piRNAs across all datasets:", length(common_genes), "\n")

# Subset to common genes + log2 if needed
clean_list <- lapply(ready_list, function(df) {
  df_sub <- df[, c("Group", common_genes)]
  mat <- df_sub[, -1]
  # Auto-detect if log2 is needed (TPM values typically > 50)
  if (max(mat, na.rm = TRUE) > 50) {
    df_sub[, -1] <- log2(mat + 1)
  }
  df_sub
})

# Build combined expression matrix for ComBat (genes x samples)
expr_matrices <- lapply(clean_list, function(df) t(as.matrix(df[, -1])))
combined_expr <- do.call(cbind, expr_matrices)

batch_vec <- unlist(lapply(names(clean_list), function(n) {
  rep(n, nrow(clean_list[[n]]))
}))
group_vec <- unlist(lapply(clean_list, function(df) as.character(df$Group)))

# Design matrix preserving biological signal (Group)
mod <- model.matrix(~ as.factor(group_vec))

cat("Running ComBat on", ncol(combined_expr), "samples across",
    length(unique(batch_vec)), "batches...\n")

combat_expr <- ComBat(dat = combined_expr, batch = batch_vec,
                      mod = mod, par.prior = TRUE)

# Reconstruct data frame (samples x genes)
combat_df_all <- data.frame(
  Group = factor(group_vec, levels = c("Normal", "Tumor")),
  t(combat_expr),
  check.names = FALSE
)
combat_df_all$Batch <- batch_vec
rownames(combat_df_all) <- colnames(combat_expr)

# Global Z-score normalization
gene_cols <- setdiff(colnames(combat_df_all), c("Group", "Batch"))
combat_df_all[, gene_cols] <- scale(combat_df_all[, gene_cols])
combat_df_all[is.na(combat_df_all)] <- 0

cat("ComBat + Z-score complete. Total:", nrow(combat_df_all), "samples x",
    length(gene_cols), "piRNAs\n")

# Post-ComBat class distribution
cat("\nPost-ComBat class distribution by batch:\n")
print(table(combat_df_all$Batch, combat_df_all$Group))


# ==============================================================================
# 4. DATA SPLITTING: DISCOVERY / HOLD-OUT / INDEPENDENT
# ==============================================================================
cat("\n========== PHASE 4: DATA SPLITTING ==========\n")

# >>> EDIT: choose your independent validation set <<<
independent_set <- "yyfbatch1"  # or "yyfbatch2"

cat("Independent validation set:", independent_set, "\n")

# Separate independent set
data_independent <- combat_df_all[combat_df_all$Batch == independent_set, ]
data_pool        <- combat_df_all[combat_df_all$Batch != independent_set, ]

# 70/30 stratified split of the pool → Discovery + Hold-out
set.seed(SEED)
train_idx <- createDataPartition(data_pool$Group, p = 0.7, list = FALSE)
data_discovery <- data_pool[train_idx, ]
data_holdout   <- data_pool[-train_idx, ]

cat("  Discovery  :", nrow(data_discovery), "samples\n")
cat("  Hold-out   :", nrow(data_holdout),   "samples\n")
cat("  Independent:", nrow(data_independent), "samples (", independent_set, ")\n")

cat("\nDiscovery class distribution:\n")
print(table(data_discovery$Group))


# ==============================================================================
# 5. FEATURE SELECTION (<10 piRNAs)
# ==============================================================================
cat("\n========== PHASE 5: FEATURE SELECTION ==========\n")

x_disc <- as.matrix(data_discovery[, gene_cols])
y_disc <- data_discovery$Group

fs_results <- list()

# --- 5.1 Differential Expression (limma) ---
cat("--- 5.1 limma ---\n")
tryCatch({
  design <- model.matrix(~ y_disc)
  fit <- lmFit(t(x_disc), design)
  fit <- eBayes(fit)
  limma_res <- topTable(fit, coef = 2, number = Inf, sort.by = "p")

  fs_limma <- rownames(limma_res)[limma_res$adj.P.Val < 0.01 & abs(limma_res$logFC) > 0.5]
  fs_results$limma <- fs_limma
  cat("  limma:", length(fs_limma), "piRNAs\n")
}, error = function(e) {
  cat("  limma failed:", conditionMessage(e), "\n")
  fs_results$limma <<- character(0)
})

# --- 5.2 Wilcoxon Rank-Sum ---
cat("--- 5.2 Wilcoxon ---\n")
tryCatch({
  pvals <- apply(x_disc, 2, function(col) wilcox.test(col ~ y_disc)$p.value)
  padj  <- p.adjust(pvals, method = "BH")
  fs_wilcoxon <- names(padj)[padj < 0.01]
  fs_results$wilcoxon <- fs_wilcoxon
  cat("  Wilcoxon:", length(fs_wilcoxon), "piRNAs\n")
}, error = function(e) {
  cat("  Wilcoxon failed:", conditionMessage(e), "\n")
  fs_results$wilcoxon <<- character(0)
})

# --- 5.3 Random Forest Importance ---
cat("--- 5.3 RF Importance ---\n")
tryCatch({
  set.seed(SEED)
  rf_fs <- randomForest(x = x_disc, y = y_disc, ntree = 1000, importance = TRUE)
  rf_imp <- importance(rf_fs, type = 1)
  rf_sorted <- sort(rf_imp[, 1], decreasing = TRUE)

  # Top piRNAs above mean + 1 SD, capped at 50
  threshold <- mean(rf_imp[, 1]) + sd(rf_imp[, 1])
  fs_rf <- names(rf_sorted[rf_sorted > threshold])
  if (length(fs_rf) > 50) fs_rf <- names(rf_sorted[1:50])

  fs_results$rf_importance <- fs_rf
  cat("  RF importance:", length(fs_rf), "piRNAs\n")
}, error = function(e) {
  cat("  RF importance failed:", conditionMessage(e), "\n")
  fs_results$rf_importance <<- character(0)
})

# --- 5.4 LASSO ---
cat("--- 5.4 LASSO ---\n")
tryCatch({
  set.seed(SEED)
  y_bin <- as.numeric(y_disc) - 1
  cv_lasso <- cv.glmnet(x_disc, y_bin, family = "binomial", alpha = 1,
                         nfolds = 10, type.measure = "auc")
  lasso_coef <- coef(cv_lasso, s = "lambda.1se")
  fs_lasso <- rownames(lasso_coef)[lasso_coef[, 1] != 0]
  fs_lasso <- setdiff(fs_lasso, "(Intercept)")

  fs_results$lasso <- fs_lasso
  cat("  LASSO:", length(fs_lasso), "piRNAs\n")
}, error = function(e) {
  cat("  LASSO failed:", conditionMessage(e), "\n")
  fs_results$lasso <<- character(0)
})

# --- 5.5 Elastic Net ---
cat("--- 5.5 Elastic Net ---\n")
tryCatch({
  set.seed(SEED)
  best_auc <- 0; best_alpha <- 0.5; best_model <- NULL
  for (a in seq(0.1, 0.9, by = 0.1)) {
    cv_en <- cv.glmnet(x_disc, y_bin, family = "binomial", alpha = a,
                        nfolds = 10, type.measure = "auc")
    if (max(cv_en$cvm) > best_auc) {
      best_auc <- max(cv_en$cvm); best_alpha <- a; best_model <- cv_en
    }
  }
  en_coef <- coef(best_model, s = "lambda.1se")
  fs_elasticnet <- rownames(en_coef)[en_coef[, 1] != 0]
  fs_elasticnet <- setdiff(fs_elasticnet, "(Intercept)")

  fs_results$elasticnet <- fs_elasticnet
  cat("  Elastic Net (alpha=", best_alpha, "):", length(fs_elasticnet), "piRNAs\n")
}, error = function(e) {
  cat("  Elastic Net failed:", conditionMessage(e), "\n")
  fs_results$elasticnet <<- character(0)
})

# --- 5.6 Boruta ---
cat("--- 5.6 Boruta ---\n")
tryCatch({
  set.seed(SEED)
  boruta_res <- Boruta(x = data.frame(x_disc), y = y_disc,
                       maxRuns = 300, doTrace = 0)
  fs_boruta <- getSelectedAttributes(boruta_res, withTentative = FALSE)
  fs_results$boruta <- fs_boruta
  cat("  Boruta:", length(fs_boruta), "piRNAs\n")
}, error = function(e) {
  cat("  Boruta failed:", conditionMessage(e), "\n")
  fs_results$boruta <<- character(0)
})

# --- 5.7 mRMR ---
cat("--- 5.7 mRMR ---\n")
tryCatch({
  k_max <- min(50, ncol(x_disc))
  mrmr_res <- MRMR(data.frame(x_disc), y_disc, k = k_max)
  fs_mrmr <- colnames(x_disc)[mrmr_res$selection]
  fs_results$mrmr <- fs_mrmr
  cat("  mRMR:", length(fs_mrmr), "piRNAs\n")
}, error = function(e) {
  cat("  mRMR failed:", conditionMessage(e), "\n")
  fs_results$mrmr <<- character(0)
})

# --- 5.8 XGBoost Gain ---
cat("--- 5.8 XGBoost Gain ---\n")
tryCatch({
  library(xgboost)
  set.seed(SEED)
  dtrain <- xgb.DMatrix(data = x_disc, label = as.numeric(y_disc) - 1)
  xgb_mod <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc"),
                        data = dtrain, nrounds = 100, verbose = 0)
  xgb_imp <- xgb.importance(model = xgb_mod)
  fs_xgb <- xgb_imp$Feature[1:min(50, nrow(xgb_imp))]
  fs_results$xgboost <- fs_xgb
  cat("  XGBoost:", length(fs_xgb), "piRNAs\n")
}, error = function(e) {
  cat("  XGBoost failed:", conditionMessage(e), "\n")
  fs_results$xgboost <<- character(0)
})


# --- 5.9 Consensus: frequency-based selection, target < 10 ---
cat("\n--- 5.9 Consensus Feature Set ---\n")

active_methods <- names(fs_results)[sapply(fs_results, length) > 0]
cat("Active methods:", length(active_methods), "\n")

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

cat("piRNAs by selection frequency (top 20):\n")
print(head(freq_table, 20))

# Forward stepwise to find optimal set < 10 features
cat("\nForward stepwise selection (target < 10 features)...\n")
ranked_pirnas <- freq_table$piRNA
best_auc_fw <- 0
no_improve <- 0
selected_fw <- c()

set.seed(SEED)
for (i in seq_along(ranked_pirnas)) {
  if (length(selected_fw) >= 9) break  # Hard cap at 9

  candidate <- c(selected_fw, ranked_pirnas[i])
  tmp_df <- data.frame(data_discovery[, candidate, drop = FALSE], Group = y_disc)

  ctrl_fw <- trainControl(method = "cv", number = 5, classProbs = TRUE,
                          summaryFunction = twoClassSummary, verboseIter = FALSE)
  mtry_val <- max(1, min(floor(sqrt(length(candidate))), length(candidate)))

  set.seed(SEED)
  tmp_model <- tryCatch({
    train(Group ~ ., data = tmp_df, method = "rf", trControl = ctrl_fw,
          metric = "ROC", tuneGrid = data.frame(mtry = mtry_val), ntree = 500)
  }, error = function(e) NULL)

  cur_auc <- if (!is.null(tmp_model)) max(tmp_model$results$ROC) else best_auc_fw

  if (cur_auc > best_auc_fw + 0.003) {
    best_auc_fw <- cur_auc
    selected_fw <- candidate
    no_improve <- 0
  } else {
    no_improve <- no_improve + 1
  }

  if (no_improve >= 3) break
}

final_features <- selected_fw
cat("\n*** FINAL piRNA SIGNATURE:", length(final_features), "features ***\n")
cat(paste(final_features, collapse = ", "), "\n")
cat("Forward stepwise CV AUC:", round(best_auc_fw, 4), "\n")

# Save feature list
write.csv(freq_table, "results/feature_selection/fs_frequency_table.csv", row.names = FALSE)
writeLines(final_features, "results/feature_selection/final_features.txt")


# ==============================================================================
# 6. MODEL TRAINING (Discovery Set)
# ==============================================================================
cat("\n========== PHASE 6: MODEL TRAINING ==========\n")

top_feats <- final_features

# Train control: 5-fold CV, save out-of-fold predictions for later ROC
fitControl <- trainControl(
  method          = "cv",
  number          = 5,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final",
  sampling        = "down"  # down-sampling inside CV for class balance
)

# Train Random Forest
set.seed(SEED)
model <- train(
  Group ~ .,
  data       = data_discovery[, c("Group", top_feats)],
  method     = "rf",
  metric     = "ROC",
  trControl  = fitControl,
  ntree      = 500,
  importance = TRUE
)

cat("Model trained.\n")
cat("Best mtry:", model$bestTune$mtry, "\n")
cat("CV AUC:", round(max(model$results$ROC), 4), "\n")

# Save model
saveRDS(model, "results/models/final_model.rds")
saveRDS(top_feats, "results/models/final_features.rds")


# ==============================================================================
# 7. EVALUATION: DISCOVERY (CV), HOLD-OUT, INDEPENDENT
# ==============================================================================
cat("\n========== PHASE 7: EVALUATION ==========\n")

# Helper: compute metrics + bootstrap CI
calc_metrics <- function(y_true01, y_score, n_boot = 1000, seed = 42) {
  roc_obj  <- roc(y_true01, y_score, direction = "<", quiet = TRUE)
  auc_val  <- as.numeric(auc(roc_obj))

  pr_obj   <- pr.curve(scores.class0 = y_score[y_true01 == 1],
                        scores.class1 = y_score[y_true01 == 0], curve = TRUE)
  auprc_val <- pr_obj$auc.integral

  set.seed(seed)
  b_auc <- b_auprc <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample(seq_along(y_true01), length(y_true01), replace = TRUE)
    yt <- y_true01[idx]; ys <- y_score[idx]
    if (length(unique(yt)) < 2) { b_auc[i] <- NA; b_auprc[i] <- NA; next }
    b_auc[i] <- as.numeric(auc(roc(yt, ys, direction = "<", quiet = TRUE)))
    pr_tmp <- pr.curve(scores.class0 = ys[yt == 1],
                        scores.class1 = ys[yt == 0], curve = FALSE)
    b_auprc[i] <- pr_tmp$auc.integral
  }
  b_auc <- na.omit(b_auc); b_auprc <- na.omit(b_auprc)

  list(roc_obj  = roc_obj,   auc   = auc_val,
       auc_ci   = quantile(b_auc, c(0.025, 0.975)),
       pr_obj   = pr_obj,    auprc = auprc_val,
       auprc_ci = quantile(b_auprc, c(0.025, 0.975)))
}

# Helper: draw ROC
plot_roc <- function(roc_obj, label, title, color = "#E41A1C") {
  df <- data.frame(fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities)
  df <- df[order(df$fpr), ]
  ggplot(df, aes(fpr, tpr)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_path(color = color, linewidth = 1.2) +
    annotate("label", x = 0.95, y = 0.05, label = label,
             hjust = 1, vjust = 0, size = 4.5, fill = "white", alpha = 0.85,
             label.size = NA) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
    labs(title = title, x = "False Positive Rate (1-Specificity)",
         y = "True Positive Rate (Sensitivity)") +
    my_theme
}

# Helper: draw PRC
plot_prc <- function(pr_obj, prevalence, label, title, color = "#377EB8") {
  df <- data.frame(recall = pr_obj$curve[, 1], precision = pr_obj$curve[, 2])
  ggplot(df, aes(recall, precision)) +
    geom_hline(yintercept = prevalence, linetype = "dashed", color = "grey50") +
    geom_path(color = color, linewidth = 1.2) +
    annotate("label", x = 0.95, y = 0.05, label = label,
             hjust = 1, vjust = 0, size = 4.5, fill = "white", alpha = 0.85,
             label.size = NA) +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1.02)) +
    labs(title = title, x = "Recall", y = "Precision") +
    my_theme
}

# --- 7.1 Training (CV out-of-fold predictions) ---
cat("--- 7.1 Training (CV) ---\n")
pred_cv <- model$pred
bt <- model$bestTune
for (nm in names(bt)) pred_cv <- pred_cv[pred_cv[[nm]] == bt[[nm]], ]

y_true_tr  <- ifelse(pred_cv$obs == "Tumor", 1, 0)
y_score_tr <- pred_cv$Tumor
mt_tr <- calc_metrics(y_true_tr, y_score_tr, n_boot = 1000, seed = 101)

lab_tr_roc <- sprintf("AUC = %.3f\n(95%% CI: %.3f-%.3f)",
                      mt_tr$auc, mt_tr$auc_ci[1], mt_tr$auc_ci[2])
p_roc_tr <- plot_roc(mt_tr$roc_obj, lab_tr_roc, "ROC - Discovery (Training, CV)")
ggsave("results/validation/ROC_discovery_CV.png", p_roc_tr, width = 7, height = 6, dpi = 300)

lab_tr_prc <- sprintf("AUPRC = %.3f\n(95%% CI: %.3f-%.3f)",
                      mt_tr$auprc, mt_tr$auprc_ci[1], mt_tr$auprc_ci[2])
p_prc_tr <- plot_prc(mt_tr$pr_obj, mean(y_true_tr), lab_tr_prc,
                     "PRC - Discovery (Training, CV)")
ggsave("results/validation/PRC_discovery_CV.png", p_prc_tr, width = 7, height = 6, dpi = 300)
print(p_roc_tr); print(p_prc_tr)

cat(sprintf("  Training CV AUC: %.3f (%.3f-%.3f)\n",
            mt_tr$auc, mt_tr$auc_ci[1], mt_tr$auc_ci[2]))

# --- 7.2 Hold-out ---
cat("--- 7.2 Hold-out ---\n")
prob_ho     <- predict(model, data_holdout[, top_feats], type = "prob")$Tumor
y_true_ho   <- ifelse(data_holdout$Group == "Tumor", 1, 0)
mt_ho <- calc_metrics(y_true_ho, prob_ho, n_boot = 1000, seed = 202)

lab_ho_roc <- sprintf("AUC = %.3f\n(95%% CI: %.3f-%.3f)",
                      mt_ho$auc, mt_ho$auc_ci[1], mt_ho$auc_ci[2])
p_roc_ho <- plot_roc(mt_ho$roc_obj, lab_ho_roc, "ROC - Hold-out (Testing)")
ggsave("results/validation/ROC_holdout.png", p_roc_ho, width = 7, height = 6, dpi = 300)

lab_ho_prc <- sprintf("AUPRC = %.3f\n(95%% CI: %.3f-%.3f)",
                      mt_ho$auprc, mt_ho$auprc_ci[1], mt_ho$auprc_ci[2])
p_prc_ho <- plot_prc(mt_ho$pr_obj, mean(y_true_ho), lab_ho_prc,
                     "PRC - Hold-out (Testing)")
ggsave("results/validation/PRC_holdout.png", p_prc_ho, width = 7, height = 6, dpi = 300)
print(p_roc_ho); print(p_prc_ho)

cat(sprintf("  Hold-out AUC: %.3f (%.3f-%.3f)\n",
            mt_ho$auc, mt_ho$auc_ci[1], mt_ho$auc_ci[2]))

# --- 7.3 Independent Validation ---
cat("--- 7.3 Independent (", independent_set, ") ---\n")
prob_iv     <- predict(model, data_independent[, top_feats], type = "prob")$Tumor
y_true_iv   <- ifelse(data_independent$Group == "Tumor", 1, 0)
mt_iv <- calc_metrics(y_true_iv, prob_iv, n_boot = 1000, seed = 303)

lab_iv_roc <- sprintf("AUC = %.3f\n(95%% CI: %.3f-%.3f)",
                      mt_iv$auc, mt_iv$auc_ci[1], mt_iv$auc_ci[2])
p_roc_iv <- plot_roc(mt_iv$roc_obj, lab_iv_roc,
                     paste0("ROC - Independent (", independent_set, ")"))
ggsave("results/validation/ROC_independent.png", p_roc_iv, width = 7, height = 6, dpi = 300)

lab_iv_prc <- sprintf("AUPRC = %.3f\n(95%% CI: %.3f-%.3f)",
                      mt_iv$auprc, mt_iv$auprc_ci[1], mt_iv$auprc_ci[2])
p_prc_iv <- plot_prc(mt_iv$pr_obj, mean(y_true_iv), lab_iv_prc,
                     paste0("PRC - Independent (", independent_set, ")"))
ggsave("results/validation/PRC_independent.png", p_prc_iv, width = 7, height = 6, dpi = 300)
print(p_roc_iv); print(p_prc_iv)

cat(sprintf("  Independent AUC: %.3f (%.3f-%.3f)\n",
            mt_iv$auc, mt_iv$auc_ci[1], mt_iv$auc_ci[2]))

# --- 7.4 Confusion Matrices ---
cat("\n--- 7.4 Confusion Matrices ---\n")

plot_confusion <- function(truth_factor, prob_tumor, title, thr = 0.5) {
  pred <- factor(ifelse(prob_tumor >= thr, "Tumor", "Normal"),
                 levels = c("Normal", "Tumor"))
  cm <- confusionMatrix(pred, truth_factor, positive = "Tumor")
  tab <- as.data.frame(cm$table)
  colnames(tab) <- c("Prediction", "Reference", "N")
  tab <- tab %>% group_by(Reference) %>%
    mutate(Pct = sprintf("%.1f%%", 100 * N / sum(N))) %>% ungroup()

  p <- ggplot(tab, aes(x = Prediction, y = Reference, fill = N)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = paste0(N, "\n(", Pct, ")")), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "steelblue") +
    scale_y_discrete(limits = rev(c("Normal", "Tumor"))) +
    labs(title = title, x = "Predicted", y = "Actual") +
    my_theme + theme(legend.position = "none")
  list(cm = cm, plot = p)
}

cm_ho <- plot_confusion(data_holdout$Group, prob_ho, "CM - Hold-out")
cm_iv <- plot_confusion(data_independent$Group, prob_iv,
                        paste0("CM - Independent (", independent_set, ")"))
ggsave("results/validation/CM_holdout.png", cm_ho$plot, width = 6, height = 5, dpi = 300)
ggsave("results/validation/CM_independent.png", cm_iv$plot, width = 6, height = 5, dpi = 300)
print(cm_ho$plot); print(cm_iv$plot)
print(cm_ho$cm); print(cm_iv$cm)


# ==============================================================================
# 8. LOGISTIC REGRESSION FOREST PLOT (Binary Variables)
# ==============================================================================
cat("\n========== PHASE 8: LOGISTIC REGRESSION + FOREST PLOT ==========\n")

# We merge BRCA1 + yyfbatch1 + yyfbatch2 data (subset of combat_df_all)
# and use the model's T-score as a binary predictor

# Compute T-score for ALL samples in combat_df_all
all_prob <- predict(model, combat_df_all[, top_feats], type = "prob")$Tumor
combat_df_all$T_Score <- all_prob

# Binarize T-score: High (>=median) vs Low (<median)
tscore_cutoff <- median(combat_df_all$T_Score)
combat_df_all$T_Score_binary <- factor(
  ifelse(combat_df_all$T_Score >= tscore_cutoff, "High", "Low"),
  levels = c("Low", "High")
)

cat("T-Score cutoff (median):", round(tscore_cutoff, 4), "\n")

# --- Load or create clinical data ---
# >>> EDIT THIS: load your real clinical data <<<
# clinical <- read.csv("clinical_data.csv", stringsAsFactors = FALSE)
# Merge: clinical must have a "SampleID" column matching rownames(combat_df_all)

# For now, create placeholder columns if clinical data isn't available
# (You should replace these with your real clinical annotations)
if (!"Age" %in% colnames(combat_df_all)) {
  cat("WARNING: No clinical data found. Using simulated clinical data for demonstration.\n")
  cat("         Replace this block with your real clinical data.\n")
  set.seed(42)
  n <- nrow(combat_df_all)
  combat_df_all$Age <- sample(30:80, n, replace = TRUE)
  combat_df_all$Stage <- NA
  is_tumor <- combat_df_all$Group == "Tumor"
  combat_df_all$Stage[is_tumor] <- sample(
    c("Stage I", "Stage II", "Stage III", "Stage IV"),
    sum(is_tumor), replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.2))
  combat_df_all$Subtype <- NA
  combat_df_all$Subtype[is_tumor] <- sample(
    c("Luminal A", "Luminal B", "HER2+", "Basal-like"),
    sum(is_tumor), replace = TRUE)
}

# Binarize clinical variables
combat_df_all$Age_binary <- factor(
  ifelse(combat_df_all$Age >= 60, ">=60", "<60"),
  levels = c("<60", ">=60")
)
combat_df_all$Stage_binary <- factor(
  ifelse(combat_df_all$Stage %in% c("Stage III", "Stage IV"), "Late", "Early"),
  levels = c("Early", "Late")
)

# Binary outcome for logistic regression
combat_df_all$Outcome01 <- ifelse(combat_df_all$Group == "Tumor", 1, 0)

# --- Univariate logistic regression ---
cat("\nUnivariate logistic regression:\n")

run_univariate_lr <- function(data, var_name, ref_level = NULL) {
  # Remove NA for this variable
  df <- data[!is.na(data[[var_name]]), ]
  if (nrow(df) < 10) return(NULL)

  fml <- as.formula(paste0("Outcome01 ~ ", var_name))
  fit <- glm(fml, data = df, family = binomial)
  s   <- summary(fit)
  ci  <- confint.default(fit)

  # Extract the non-intercept term(s)
  coef_rows <- rownames(s$coefficients)[-1]  # skip Intercept
  results <- lapply(coef_rows, function(term) {
    or    <- exp(s$coefficients[term, "Estimate"])
    lower <- exp(ci[term, 1])
    upper <- exp(ci[term, 2])
    pval  <- s$coefficients[term, "Pr(>|z|)"]
    data.frame(Variable = term, OR = or, Lower = lower, Upper = upper,
               P = pval, stringsAsFactors = FALSE)
  })
  do.call(rbind, results)
}

# Variables to test
uni_vars <- c("T_Score_binary", "Age_binary", "Stage_binary")
uni_results <- do.call(rbind, lapply(uni_vars, function(v) {
  run_univariate_lr(combat_df_all, v)
}))

if (!is.null(uni_results)) {
  uni_results$Sig <- ifelse(uni_results$P < 0.001, "***",
                     ifelse(uni_results$P < 0.01, "**",
                     ifelse(uni_results$P < 0.05, "*", "ns")))
  cat("\nUnivariate results:\n")
  print(uni_results)
  write.csv(uni_results, "results/forest_plot/univariate_results.csv", row.names = FALSE)
}

# --- Multivariate logistic regression ---
cat("\nMultivariate logistic regression:\n")

# Use only samples with complete clinical data
df_multi <- combat_df_all[complete.cases(combat_df_all[, c("Outcome01", uni_vars)]), ]
fml_multi <- as.formula(paste0("Outcome01 ~ ", paste(uni_vars, collapse = " + ")))
fit_multi <- glm(fml_multi, data = df_multi, family = binomial)
s_multi   <- summary(fit_multi)
ci_multi  <- confint.default(fit_multi)

multi_results <- data.frame()
coef_rows <- rownames(s_multi$coefficients)[-1]
for (term in coef_rows) {
  multi_results <- rbind(multi_results, data.frame(
    Variable = term,
    OR    = exp(s_multi$coefficients[term, "Estimate"]),
    Lower = exp(ci_multi[term, 1]),
    Upper = exp(ci_multi[term, 2]),
    P     = s_multi$coefficients[term, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  ))
}
multi_results$Sig <- ifelse(multi_results$P < 0.001, "***",
                    ifelse(multi_results$P < 0.01, "**",
                    ifelse(multi_results$P < 0.05, "*", "ns")))

cat("Multivariate results:\n")
print(multi_results)
write.csv(multi_results, "results/forest_plot/multivariate_results.csv", row.names = FALSE)

# --- Forest Plot (ggplot2 version) ---
cat("\nGenerating forest plot...\n")

# Combine univariate and multivariate
uni_results$Analysis    <- "Univariate"
multi_results$Analysis  <- "Multivariate"
forest_data <- rbind(uni_results, multi_results)

# Clean variable names for display
forest_data$Label <- gsub("T_Score_binary", "T-Score (High vs Low)", forest_data$Variable)
forest_data$Label <- gsub("Age_binary", "Age (>=60 vs <60)", forest_data$Label)
forest_data$Label <- gsub("Stage_binary", "Stage (Late vs Early)", forest_data$Label)
# Remove factor level suffix from variable names
forest_data$Label <- gsub("High$|>=60$|Late$", "", forest_data$Label)
forest_data$Label <- trimws(forest_data$Label)

# Deduplicate labels
forest_data$Label[duplicated(forest_data$Label) & forest_data$Analysis == "Multivariate"] <-
  paste0(forest_data$Label[duplicated(forest_data$Label) & forest_data$Analysis == "Multivariate"], " ")

forest_data$Analysis <- factor(forest_data$Analysis,
                               levels = c("Univariate", "Multivariate"))

# Build annotation text
forest_data$OR_text <- sprintf("%.2f (%.2f-%.2f)", forest_data$OR,
                               forest_data$Lower, forest_data$Upper)
forest_data$P_text <- ifelse(forest_data$P < 0.001,
                             sprintf("%.1e", forest_data$P),
                             sprintf("%.3f", forest_data$P))

# Y-axis ordering
forest_data$y_pos <- seq(nrow(forest_data), 1)

p_forest <- ggplot(forest_data, aes(x = OR, y = y_pos, color = Analysis)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_errorbarh(aes(xmin = Lower, xmax = Upper),
                 height = 0.2, linewidth = 0.8,
                 position = position_dodge(width = 0.5)) +
  scale_y_continuous(breaks = forest_data$y_pos,
                     labels = paste0(forest_data$Label, " [",
                                     forest_data$Analysis, "]")) +
  scale_color_manual(values = c("Univariate" = "#E41A1C", "Multivariate" = "#377EB8")) +
  scale_x_log10() +
  # Add OR (CI) and p-value as text annotations
  geom_text(aes(x = max(forest_data$Upper, na.rm = TRUE) * 1.5,
                label = paste0(OR_text, "  p=", P_text)),
            hjust = 0, size = 3.2, show.legend = FALSE,
            position = position_dodge(width = 0.5)) +
  labs(title = "Forest Plot: Logistic Regression",
       subtitle = "Odds Ratio (95% CI) for Tumor vs Normal",
       x = "Odds Ratio (log scale)", y = "", color = "Analysis") +
  my_theme +
  theme(legend.position = "top",
        axis.text.y = element_text(size = 10, hjust = 1),
        plot.margin = margin(10, 120, 10, 10))  # extra right margin for text

ggsave("results/forest_plot/forest_plot.png", p_forest,
       width = 12, height = max(4, nrow(forest_data) * 0.6 + 2), dpi = 300)
print(p_forest)

cat("Forest plot saved.\n")


# ==============================================================================
# 9. SUBGROUP ANALYSIS: ROC CURVES
#    (Merge BRCA1 + yyfbatch1 + yyfbatch2 → stratify by Age, Stage, Subtype)
# ==============================================================================
cat("\n========== PHASE 9: SUBGROUP ROC ANALYSIS ==========\n")

# Subset to BRCA1 + yyfbatch1 + yyfbatch2
target_batches <- c("BRCA1", "yyfbatch1", "yyfbatch2")
sub_df <- combat_df_all[combat_df_all$Batch %in% target_batches, ]

cat("Merged subset (BRCA1 + yyfbatch1 + yyfbatch2):", nrow(sub_df), "samples\n")
print(table(sub_df$Batch, sub_df$Group))

# Ensure T-score is computed for all
if (!"T_Score" %in% colnames(sub_df)) {
  sub_df$T_Score <- predict(model, sub_df[, top_feats], type = "prob")$Tumor
}

# Subgroup ROC function:
# For each subgroup value, compare: Tumor(in that subgroup) vs ALL Normals
plot_subgroup_roc <- function(df, group_col, title_text) {
  df_normal <- df[df$Group == "Normal", ]
  df_tumor  <- df[df$Group == "Tumor", ]

  subgroups <- sort(unique(na.omit(df_tumor[[group_col]])))

  roc_data_all <- data.frame()
  legend_labels <- c()

  for (grp in subgroups) {
    sub_tumor <- df_tumor[df_tumor[[group_col]] == grp & !is.na(df_tumor[[group_col]]), ]
    if (nrow(sub_tumor) < 5) next

    combined <- rbind(sub_tumor, df_normal)
    y01 <- ifelse(combined$Group == "Tumor", 1, 0)
    roc_obj <- roc(y01, combined$T_Score, direction = "<", quiet = TRUE)
    auc_val <- round(as.numeric(auc(roc_obj)), 3)

    roc_df <- data.frame(
      fpr      = 1 - roc_obj$specificities,
      tpr      = roc_obj$sensitivities,
      Subgroup = paste0(grp, " (AUC=", auc_val, ", n=", nrow(sub_tumor), ")")
    )
    roc_data_all <- rbind(roc_data_all, roc_df)
  }

  if (nrow(roc_data_all) == 0) {
    cat("  No valid subgroups for", group_col, "\n")
    return(NULL)
  }

  ggplot(roc_data_all, aes(x = fpr, y = tpr, color = Subgroup)) +
    geom_path(linewidth = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    labs(title = title_text,
         x = "1 - Specificity", y = "Sensitivity", color = "") +
    my_theme +
    theme(legend.position = c(0.65, 0.25),
          legend.background = element_rect(fill = "white", color = "grey80"),
          legend.text = element_text(size = 10))
}

# --- 9.1 By Age Group ---
sub_df$Age_Group <- ifelse(sub_df$Age >= 60, "Age >= 60", "Age < 60")
p_roc_age <- plot_subgroup_roc(sub_df, "Age_Group", "ROC Stratified by Age Group")
if (!is.null(p_roc_age)) {
  ggsave("results/subgroup/ROC_by_Age.png", p_roc_age, width = 8, height = 7, dpi = 300)
  print(p_roc_age)
}

# --- 9.2 By Stage ---
p_roc_stage <- plot_subgroup_roc(sub_df, "Stage", "ROC Stratified by Cancer Stage")
if (!is.null(p_roc_stage)) {
  ggsave("results/subgroup/ROC_by_Stage.png", p_roc_stage, width = 8, height = 7, dpi = 300)
  print(p_roc_stage)
}

# --- 9.3 By Subtype ---
p_roc_subtype <- plot_subgroup_roc(sub_df, "Subtype", "ROC Stratified by Molecular Subtype")
if (!is.null(p_roc_subtype)) {
  ggsave("results/subgroup/ROC_by_Subtype.png", p_roc_subtype, width = 8, height = 7, dpi = 300)
  print(p_roc_subtype)
}


# ==============================================================================
# 10. T-SCORE BOXPLOTS BY SUBGROUP
# ==============================================================================
cat("\n========== PHASE 10: T-SCORE BOXPLOTS ==========\n")

# Subgroup boxplot function:
# For each subgroup level, show Tumor vs Normal T-scores side by side
plot_subgroup_boxplot <- function(df, subgroup_col, title_text) {
  df_normal <- df[df$Group == "Normal", ]
  df_tumor  <- df[df$Group == "Tumor", ]

  valid_subgroups <- sort(unique(na.omit(df_tumor[[subgroup_col]])))

  plot_data <- data.frame()
  for (grp in valid_subgroups) {
    sub_tumor <- df_tumor[df_tumor[[subgroup_col]] == grp, ]
    if (nrow(sub_tumor) < 3) next

    # Normal samples use all normals, tagged with this subgroup label
    tmp_normal <- df_normal
    tmp_normal[[subgroup_col]] <- grp
    plot_data <- rbind(plot_data, rbind(sub_tumor, tmp_normal))
  }

  if (nrow(plot_data) == 0) return(NULL)

  plot_data[[subgroup_col]] <- factor(plot_data[[subgroup_col]], levels = valid_subgroups)
  my_colors <- c("Normal" = "#868686FF", "Tumor" = "#CD534CFF")

  ggplot(plot_data, aes_string(x = subgroup_col, y = "T_Score", fill = "Group")) +
    stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(0.8)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8,
                 position = position_dodge(0.8)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
                size = 0.8, alpha = 0.3, color = "black") +
    stat_compare_means(aes(group = Group), method = "wilcox.test",
                       label = "p.signif", label.y = 1.05, size = 5) +
    scale_fill_manual(values = my_colors) +
    scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2)) +
    labs(title = title_text, x = "", y = "Predicted T-Score (P(Tumor))") +
    my_theme +
    theme(axis.text.x = element_text(size = 11, face = "bold", angle = 15, hjust = 1),
          legend.position = "top", legend.title = element_blank(),
          panel.grid = element_blank())
}

# --- 10.1 T-Score by Age ---
p_box_age <- plot_subgroup_boxplot(sub_df, "Age_Group", "T-Score by Age Group")
if (!is.null(p_box_age)) {
  ggsave("results/subgroup/Boxplot_TScore_Age.png", p_box_age,
         width = 7, height = 6, dpi = 300)
  print(p_box_age)
}

# --- 10.2 T-Score by Stage ---
# Ensure stage ordering
if ("Stage" %in% colnames(sub_df)) {
  valid_stages <- c("Stage I", "Stage II", "Stage III", "Stage IV")
  sub_df$Stage <- factor(sub_df$Stage, levels = intersect(valid_stages, unique(sub_df$Stage)))
}
p_box_stage <- plot_subgroup_boxplot(sub_df, "Stage", "T-Score by Cancer Stage")
if (!is.null(p_box_stage)) {
  ggsave("results/subgroup/Boxplot_TScore_Stage.png", p_box_stage,
         width = 8, height = 6, dpi = 300)
  print(p_box_stage)
}

# --- 10.3 T-Score by Subtype ---
p_box_subtype <- plot_subgroup_boxplot(sub_df, "Subtype", "T-Score by Molecular Subtype")
if (!is.null(p_box_subtype)) {
  ggsave("results/subgroup/Boxplot_TScore_Subtype.png", p_box_subtype,
         width = 8, height = 6, dpi = 300)
  print(p_box_subtype)
}

# --- 10.4 T-Score by Phase (Discovery / Hold-out / Independent) ---
cat("\nGenerating T-Score boxplot by study phase...\n")

score_phase <- function(df, phase_name) {
  prob <- predict(model, df[, top_feats], type = "prob")$Tumor
  data.frame(
    Group   = df$Group,
    T_Score = prob,
    Phase   = phase_name,
    stringsAsFactors = FALSE
  )
}

phase_scores <- rbind(
  score_phase(data_discovery,   "Discovery"),
  score_phase(data_holdout,     "Hold-out"),
  score_phase(data_independent, "Independent")
)
phase_scores$Phase <- factor(phase_scores$Phase,
                             levels = c("Discovery", "Hold-out", "Independent"))

my_colors_2 <- c("Normal" = "#868686FF", "Tumor" = "#CD534CFF")

p_box_phase <- ggplot(phase_scores, aes(x = Phase, y = T_Score, fill = Group)) +
  stat_boxplot(geom = "errorbar", width = 0.2) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2),
              size = 1, alpha = 0.3, color = "black") +
  stat_compare_means(aes(group = Group), method = "wilcox.test",
                     label = "p.signif", label.y = 1.05, size = 6) +
  scale_fill_manual(values = my_colors_2) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.2)) +
  labs(title = "T-Score: Tumor vs Normal across Study Phases",
       subtitle = paste("Independent set:", independent_set),
       x = "", y = "Predicted T-Score (P(Tumor))") +
  my_theme +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "top", legend.title = element_blank(),
        panel.grid = element_blank())

ggsave("results/subgroup/Boxplot_TScore_Phase.png", p_box_phase,
       width = 8, height = 6, dpi = 300)
print(p_box_phase)


# ==============================================================================
# 11. FINAL SUMMARY
# ==============================================================================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  FINAL SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("piRNA Signature (", length(top_feats), "features ):\n")
cat("  ", paste(top_feats, collapse = ", "), "\n\n")

cat("Model: Random Forest (down-sampled CV)\n")
cat(sprintf("  Discovery CV AUC:  %.3f (95%% CI: %.3f-%.3f)\n",
            mt_tr$auc, mt_tr$auc_ci[1], mt_tr$auc_ci[2]))
cat(sprintf("  Hold-out AUC:      %.3f (95%% CI: %.3f-%.3f)\n",
            mt_ho$auc, mt_ho$auc_ci[1], mt_ho$auc_ci[2]))
cat(sprintf("  Independent AUC:   %.3f (95%% CI: %.3f-%.3f) [%s]\n",
            mt_iv$auc, mt_iv$auc_ci[1], mt_iv$auc_ci[2], independent_set))

cat("\nBatch correction: Global ComBat (all 7 datasets)\n")
cat("BRCA1 balancing: Matched pairs + 40% remaining tumors\n")
cat("T-Score cutoff (median):", round(tscore_cutoff, 4), "\n")

end_time <- Sys.time()
cat("\nTotal runtime:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")
cat("\nAll results saved to results/ directory.\n")

# Write final piRNA list
write.csv(data.frame(piRNA = top_feats), "Final_piRNA_Signature.csv", row.names = FALSE)

cat("\n*** PIPELINE COMPLETE ***\n")

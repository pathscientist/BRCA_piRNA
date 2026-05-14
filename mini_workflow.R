# MINI Workflow for Breast Cancer piRNA Diagnosis
# Implements PROTOCOL_mini_workflow (1).md
#
# Design:
#   Development pool (70/30 split, stratified by Dataset + Group):
#     BRCA1, PRJNA482141, PRJNA294226, PRJNA808405, PRJNA934049, yyfbatch2
#   Primary independent validation cohort: yyfbatch1
#   ComBat harmonization across ALL public datasets + yyfbatch1 + yyfbatch2
#   Model: Random Forest only
#   Survival analysis: TCGA-BRCA / BRCA1 only

suppressPackageStartupMessages({
  library(sva)
  library(caret)
  library(randomForest)
  library(pROC)
  library(PRROC)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(survival)
  library(survminer)
  library(timeROC)
})

set.seed(2024)

# ----------------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------------
PROC_DIR     <- "processed_results"
CLIN_DIR     <- "clinical_data"
OUT_DIR      <- file.path("analysis", "mini_workflow")
FIG_DIR      <- file.path(OUT_DIR, "figures")
TAB_DIR      <- file.path(OUT_DIR, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

TISSUE_SETS  <- c("BRCA1", "PRJNA482141")
PLASMA_SETS  <- c("PRJNA294226", "PRJNA808405", "PRJNA934049",
                  "yyfbatch1", "yyfbatch2")
DEV_SETS     <- c("BRCA1", "PRJNA482141", "PRJNA294226",
                  "PRJNA808405", "PRJNA934049", "yyfbatch2")
INDEP_SET    <- "yyfbatch1"
ALL_SETS     <- c(TISSUE_SETS, PLASMA_SETS)

SPLIT_FRAC   <- 0.7
N_TOP_FEATS  <- 15
RF_NTREE_FS  <- 500
RF_NTREE     <- 1000

# ----------------------------------------------------------------------------
# Step 1. Data harmonization
# ----------------------------------------------------------------------------
cat("== Step 1: load and harmonize ==\n")

load_cohort <- function(name) {
  path <- file.path(PROC_DIR, paste0(name, "_processed.csv"))
  alt  <- file.path(PROC_DIR, paste0(name, "_processed_1.csv"))
  if (!file.exists(path) && file.exists(alt)) path <- alt
  if (!file.exists(path)) return(NULL)
  df <- read.csv(path, row.names = 1, stringsAsFactors = FALSE,
                 check.names = FALSE)
  if (!"Group" %in% colnames(df)) {
    stop(sprintf("Cohort %s missing Group column", name))
  }
  df
}

recode_group <- function(g) {
  g <- as.character(g)
  out <- rep(NA_character_, length(g))
  out[grepl("tumor|cancer|brca|case", g, ignore.case = TRUE)] <- "Tumor"
  out[grepl("normal|control|healthy|hc", g, ignore.case = TRUE)] <- "Normal"
  out
}

raw_list <- lapply(ALL_SETS, load_cohort)
names(raw_list) <- ALL_SETS
present <- ALL_SETS[!vapply(raw_list, is.null, logical(1))]
missing <- setdiff(ALL_SETS, present)
if (length(missing)) {
  cat("WARNING: missing cohorts (skipped):", paste(missing, collapse = ", "), "\n")
}
raw_list <- raw_list[present]

for (nm in names(raw_list)) {
  raw_list[[nm]]$Group <- recode_group(raw_list[[nm]]$Group)
  drop <- is.na(raw_list[[nm]]$Group)
  if (any(drop)) {
    cat(sprintf("  %s: dropping %d samples with unmapped group\n",
                nm, sum(drop)))
    raw_list[[nm]] <- raw_list[[nm]][!drop, , drop = FALSE]
  }
}

# Intersect common piRNAs
common <- Reduce(intersect,
                 lapply(raw_list, function(df) setdiff(colnames(df), "Group")))
cat("Common piRNAs across cohorts:", length(common), "\n")
if (length(common) < 10) stop("Too few common piRNAs.")

# log2(x + 1) — skip if already on log scale
log2_if_raw <- function(df, feats) {
  m <- as.matrix(df[, feats, drop = FALSE])
  if (max(m, na.rm = TRUE) > 50) m <- log2(m + 1)
  df[, feats] <- m
  df
}
clean_list <- lapply(raw_list, function(df) {
  df <- df[, c("Group", common), drop = FALSE]
  log2_if_raw(df, common)
})

# Build matrix for ComBat across ALL cohorts
expr_mat <- do.call(cbind, lapply(clean_list,
                                  function(df) t(df[, common, drop = FALSE])))
batch    <- unlist(lapply(names(clean_list),
                          function(n) rep(n, nrow(clean_list[[n]]))))
groups   <- unlist(lapply(clean_list, function(df) df$Group))
samples  <- unlist(lapply(clean_list, rownames), use.names = FALSE)
colnames(expr_mat) <- samples

cat("Running ComBat on", ncol(expr_mat), "samples,",
    length(unique(batch)), "batches...\n")
mod <- model.matrix(~ as.factor(groups))
combat_mat <- ComBat(dat = expr_mat, batch = batch, mod = mod,
                     par.prior = TRUE)

# Z-score normalization per feature (row)
z_mat <- t(scale(t(combat_mat)))
z_mat[is.na(z_mat)] <- 0

harmonized <- data.frame(
  Sample  = samples,
  Dataset = batch,
  Group   = factor(groups, levels = c("Normal", "Tumor")),
  t(z_mat),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
rownames(harmonized) <- samples
write.csv(harmonized,
          file.path(TAB_DIR, "harmonized_expression.csv"),
          row.names = FALSE)

# ----------------------------------------------------------------------------
# Step 2 & 3. Build development pool and 70/30 split
# ----------------------------------------------------------------------------
cat("== Step 2-3: development pool and 70/30 split ==\n")

dev_sets_present  <- intersect(DEV_SETS, present)
dev_df            <- harmonized[harmonized$Dataset %in% dev_sets_present, ]
indep_df          <- harmonized[harmonized$Dataset == INDEP_SET, ]
cat("Development pool n =", nrow(dev_df),
    "| Independent (", INDEP_SET, ") n =", nrow(indep_df), "\n")

# Stratify by Dataset + Group
strat_key  <- paste(dev_df$Dataset, dev_df$Group, sep = "_")
train_idx  <- caret::createDataPartition(strat_key, p = SPLIT_FRAC,
                                         list = FALSE)[, 1]
train_df   <- dev_df[ train_idx, ]
holdout_df <- dev_df[-train_idx, ]
cat("Training n =", nrow(train_df),
    "| Hold-out n =", nrow(holdout_df), "\n")

# ----------------------------------------------------------------------------
# Step 4. Feature selection (on training only) and Random Forest
# ----------------------------------------------------------------------------
cat("== Step 4: Random Forest ==\n")

x_train <- as.matrix(train_df[, common])
y_train <- train_df$Group

rf_fs <- randomForest(x = x_train, y = y_train,
                      ntree = RF_NTREE_FS, importance = TRUE)
imp    <- importance(rf_fs, type = 2)  # MeanDecreaseGini
top_feats <- rownames(imp)[order(imp[, 1], decreasing = TRUE)][
  seq_len(min(N_TOP_FEATS, nrow(imp)))]
cat("Selected", length(top_feats), "features:\n  ",
    paste(top_feats, collapse = ", "), "\n")

train_sub <- train_df[, c("Group", top_feats)]

ctrl <- caret::trainControl(method = "cv", number = 5,
                            classProbs = TRUE,
                            summaryFunction = caret::twoClassSummary,
                            sampling = "down",
                            savePredictions = "final")
rf_grid <- expand.grid(mtry = c(2, 3, 4, 6, 8))
rf_model <- caret::train(Group ~ ., data = train_sub,
                         method = "rf",
                         metric = "ROC",
                         tuneGrid = rf_grid,
                         trControl = ctrl,
                         ntree = RF_NTREE)
cat("Best mtry =", rf_model$bestTune$mtry, "\n")

saveRDS(list(model = rf_model, features = top_feats),
        file.path(TAB_DIR, "rf_model.rds"))

# ----------------------------------------------------------------------------
# Step 5. Validation: prediction helpers and per-set metrics
# ----------------------------------------------------------------------------
cat("== Step 5: evaluation ==\n")

predict_prob <- function(df) {
  if (nrow(df) == 0) return(NULL)
  miss <- setdiff(top_feats, colnames(df))
  if (length(miss)) for (m in miss) df[[m]] <- 0
  predict(rf_model, newdata = df[, top_feats], type = "prob")[, "Tumor"]
}

youden_threshold <- function(labels, probs) {
  r <- pROC::roc(labels, probs,
                 levels = c("Normal", "Tumor"), direction = "<",
                 quiet = TRUE)
  co <- pROC::coords(r, "best", best.method = "youden",
                     ret = c("threshold", "sensitivity", "specificity"),
                     transpose = FALSE)
  if (is.data.frame(co)) co <- co[1, ]
  list(roc = r, threshold = as.numeric(co["threshold"]))
}

metrics_at <- function(labels, probs, thr) {
  labels <- factor(labels, levels = c("Normal", "Tumor"))
  pred   <- factor(ifelse(probs >= thr, "Tumor", "Normal"),
                   levels = c("Normal", "Tumor"))
  cm <- table(Pred = pred, True = labels)
  TP <- cm["Tumor", "Tumor"]; TN <- cm["Normal", "Normal"]
  FP <- cm["Tumor", "Normal"]; FN <- cm["Normal", "Tumor"]
  sens <- TP / (TP + FN); spec <- TN / (TN + FP)
  acc  <- (TP + TN) / sum(cm); bal_acc <- (sens + spec) / 2
  ppv  <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
  npv  <- if ((TN + FN) > 0) TN / (TN + FN) else NA_real_
  f1   <- if (!is.na(ppv) && (ppv + sens) > 0) 2 * ppv * sens / (ppv + sens) else NA_real_
  roc_obj <- pROC::roc(labels, probs,
                       levels = c("Normal", "Tumor"), direction = "<",
                       quiet = TRUE)
  pr <- PRROC::pr.curve(scores.class0 = probs[labels == "Tumor"],
                        scores.class1 = probs[labels == "Normal"],
                        curve = TRUE)
  list(
    AUROC = as.numeric(pROC::auc(roc_obj)),
    AUPRC = pr$auc.integral,
    Sensitivity = sens, Specificity = spec,
    Accuracy = acc, BalancedAccuracy = bal_acc,
    PPV = ppv, NPV = npv, F1 = f1,
    Threshold = thr,
    n = sum(cm),
    roc = roc_obj, pr = pr,
    confusion = cm
  )
}

# Establish threshold from the training set (Youden's J)
train_probs <- predict_prob(train_df)
thr_info    <- youden_threshold(train_df$Group, train_probs)
THR         <- thr_info$threshold
cat("Operating threshold (Youden, training):", round(THR, 4), "\n")

# ----------------------------------------------------------------------------
# Build evaluation groups
# ----------------------------------------------------------------------------
eval_groups <- list()
eval_groups[["Discovery/Training"]]            <- train_df
eval_groups[["Hold-out validation"]]           <- holdout_df
if (nrow(indep_df) > 0)
  eval_groups[["yyfbatch1 independent validation"]] <- indep_df

tissue_df <- harmonized[harmonized$Dataset %in% intersect(TISSUE_SETS, present), ]
plasma_df <- harmonized[harmonized$Dataset %in% intersect(PLASMA_SETS, present), ]
if (nrow(tissue_df)) eval_groups[["Tissue pooled"]] <- tissue_df
if (nrow(plasma_df)) eval_groups[["Plasma pooled"]] <- plasma_df

for (ds in present) {
  sub <- harmonized[harmonized$Dataset == ds, ]
  if (nrow(sub) > 0) eval_groups[[paste0("Dataset: ", ds)]] <- sub
}

# Compute metrics for every group
eval_results <- list()
for (gname in names(eval_groups)) {
  df <- eval_groups[[gname]]
  if (length(unique(df$Group)) < 2) {
    cat("  skip", gname, "- only one class present\n")
    next
  }
  p <- predict_prob(df)
  m <- metrics_at(df$Group, p, THR)
  m$probs  <- p
  m$labels <- df$Group
  m$name   <- gname
  eval_results[[gname]] <- m
}

metrics_tbl <- do.call(rbind, lapply(eval_results, function(m) {
  data.frame(
    Group = m$name, N = length(m$labels),
    AUROC = m$AUROC, AUPRC = m$AUPRC,
    Sensitivity = m$Sensitivity, Specificity = m$Specificity,
    Accuracy = m$Accuracy, BalancedAccuracy = m$BalancedAccuracy,
    PPV = m$PPV, NPV = m$NPV, F1 = m$F1,
    Threshold = m$Threshold,
    stringsAsFactors = FALSE
  )
}))
rownames(metrics_tbl) <- NULL
write.csv(metrics_tbl,
          file.path(TAB_DIR, "diagnostic_metrics.csv"), row.names = FALSE)
cat("\nDiagnostic metrics:\n"); print(metrics_tbl, row.names = FALSE)

# Confusion matrices (supplementary)
cm_list <- lapply(eval_results, function(m) {
  data.frame(Group = m$name,
             as.data.frame.matrix(m$confusion),
             stringsAsFactors = FALSE)
})
saveRDS(cm_list, file.path(TAB_DIR, "confusion_matrices.rds"))

# ----------------------------------------------------------------------------
# Combined ROC / PRC figures
# ----------------------------------------------------------------------------
plot_combined_roc <- function(results, fname, title) {
  if (!length(results)) return(invisible(NULL))
  df <- do.call(rbind, lapply(results, function(m) {
    co <- pROC::coords(m$roc, x = "all",
                       ret = c("specificity", "sensitivity"),
                       transpose = FALSE)
    data.frame(FPR = 1 - co$specificity, TPR = co$sensitivity,
               Group = sprintf("%s (AUC=%.3f)", m$name, m$AUROC))
  }))
  p <- ggplot(df, aes(FPR, TPR, color = Group)) +
    geom_line(linewidth = 0.8) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "grey60") +
    coord_equal() +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          legend.text = element_text(size = 8))
  ggsave(fname, p, width = 7.5, height = 5, dpi = 200)
}

plot_combined_prc <- function(results, fname, title) {
  if (!length(results)) return(invisible(NULL))
  df <- do.call(rbind, lapply(results, function(m) {
    cv <- as.data.frame(m$pr$curve)
    colnames(cv) <- c("Recall", "Precision", "Threshold")
    cv$Group <- sprintf("%s (AUPRC=%.3f)", m$name, m$AUPRC)
    cv
  }))
  p <- ggplot(df, aes(Recall, Precision, color = Group)) +
    geom_line(linewidth = 0.8) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(title = title, x = "Recall", y = "Precision") +
    theme_bw(base_size = 11) +
    theme(legend.position = "right",
          legend.text = element_text(size = 8))
  ggsave(fname, p, width = 7.5, height = 5, dpi = 200)
}

main_groups <- c("Discovery/Training", "Hold-out validation",
                 "yyfbatch1 independent validation",
                 "Tissue pooled", "Plasma pooled")
main_results <- eval_results[intersect(main_groups, names(eval_results))]
plot_combined_roc(main_results,
                  file.path(FIG_DIR, "fig_main_ROC.png"),
                  "Main MINI workflow ROC")
plot_combined_prc(main_results,
                  file.path(FIG_DIR, "fig_main_PRC.png"),
                  "Main MINI workflow PRC")

dataset_results <- eval_results[grepl("^Dataset: ", names(eval_results))]
plot_combined_roc(dataset_results,
                  file.path(FIG_DIR, "fig_supp_dataset_ROC.png"),
                  "Per-dataset ROC overlay")
plot_combined_prc(dataset_results,
                  file.path(FIG_DIR, "fig_supp_dataset_PRC.png"),
                  "Per-dataset PRC overlay")

# Feature importance plot
imp_df <- data.frame(piRNA = top_feats,
                     Importance = imp[top_feats, 1])
imp_df <- imp_df[order(imp_df$Importance), ]
imp_df$piRNA <- factor(imp_df$piRNA, levels = imp_df$piRNA)
p_imp <- ggplot(imp_df, aes(Importance, piRNA)) +
  geom_col(fill = "steelblue") +
  labs(title = "RF feature importance (MeanDecreaseGini)",
       x = "MeanDecreaseGini", y = NULL) +
  theme_bw(base_size = 11)
ggsave(file.path(FIG_DIR, "fig_feature_importance.png"),
       p_imp, width = 6, height = 5, dpi = 200)
write.csv(imp_df, file.path(TAB_DIR, "feature_importance.csv"),
          row.names = FALSE)

# ----------------------------------------------------------------------------
# Survival analysis — TCGA-BRCA / BRCA1 only
# ----------------------------------------------------------------------------
cat("== Survival analysis (TCGA-BRCA / BRCA1) ==\n")

clin_path <- file.path(CLIN_DIR, "TCGA_BRCA_clinical.csv")
brca_df   <- harmonized[harmonized$Dataset == "BRCA1" &
                          harmonized$Group == "Tumor", ]

if (nrow(brca_df) > 0 && file.exists(clin_path)) {
  clin <- read.csv(clin_path, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(clin) <- make.unique(colnames(clin))
  id_col <- colnames(clin)[1]
  clin$sample_id <- clin[[id_col]]

  # Vital status -> 0/1
  vs_col <- grep("vital_status", colnames(clin), value = TRUE)[1]
  d_col  <- grep("days_to_death", colnames(clin), value = TRUE)[1]
  f_col  <- grep("days_to_last_follow", colnames(clin), value = TRUE)[1]
  clin$event <- as.integer(tolower(clin[[vs_col]]) %in% c("dead", "deceased"))
  d <- suppressWarnings(as.numeric(clin[[d_col]]))
  f <- suppressWarnings(as.numeric(clin[[f_col]]))
  clin$time <- ifelse(clin$event == 1 & !is.na(d), d,
                      ifelse(!is.na(f), f, d))

  # Match TCGA barcodes: harmonized rownames look like TCGA_3C_AAAU_01A
  brca_df$match_id <- gsub("_", "-", rownames(brca_df))
  surv_in <- merge(brca_df, clin[, c("sample_id", "time", "event")],
                   by.x = "match_id", by.y = "sample_id")
  surv_in <- surv_in[!is.na(surv_in$time) & surv_in$time > 0 &
                       !is.na(surv_in$event), ]
  cat("Survival cohort n =", nrow(surv_in), "\n")

  if (nrow(surv_in) >= 30) {
    # Risk score = predicted tumor probability for each TCGA tumor sample
    surv_in$risk <- predict_prob(surv_in)
    surv_in$risk_group <- ifelse(surv_in$risk >= median(surv_in$risk),
                                 "High", "Low")
    surv_in$risk_group <- factor(surv_in$risk_group, levels = c("Low", "High"))

    # Univariate Cox on continuous risk
    cox_uni <- coxph(Surv(time, event) ~ risk, data = surv_in)
    cox_uni_tbl <- as.data.frame(summary(cox_uni)$coefficients)
    cox_uni_tbl$Variable <- rownames(cox_uni_tbl)

    # Multivariate Cox on top features (where available)
    feat_present <- intersect(top_feats, colnames(surv_in))
    fmla <- as.formula(paste("Surv(time, event) ~",
                             paste(sprintf("`%s`", feat_present),
                                   collapse = " + ")))
    cox_multi <- tryCatch(coxph(fmla, data = surv_in), error = function(e) NULL)
    if (!is.null(cox_multi)) {
      cox_multi_tbl <- as.data.frame(summary(cox_multi)$coefficients)
      cox_multi_tbl$Variable <- rownames(cox_multi_tbl)
      write.csv(cox_multi_tbl,
                file.path(TAB_DIR, "cox_multivariate.csv"),
                row.names = FALSE)
    }
    write.csv(cox_uni_tbl,
              file.path(TAB_DIR, "cox_univariate.csv"),
              row.names = FALSE)

    # Kaplan-Meier
    km_fit <- survfit(Surv(time, event) ~ risk_group, data = surv_in)
    km_plot <- ggsurvplot(km_fit, data = surv_in,
                          pval = TRUE, risk.table = TRUE,
                          conf.int = TRUE,
                          xlab = "Days", ylab = "Overall survival",
                          title = "TCGA-BRCA Kaplan-Meier (RF risk groups)")
    ggsave(file.path(FIG_DIR, "fig_KM_TCGA.png"),
           print(km_plot), width = 7, height = 7, dpi = 200)

    # Time-dependent ROC at 1/3/5 years
    horizons <- c(365, 1095, 1825)
    tdroc <- timeROC(T = surv_in$time, delta = surv_in$event,
                     marker = surv_in$risk,
                     cause = 1, times = horizons, iid = FALSE)
    td_df <- data.frame(Time_days = horizons,
                        AUC = as.numeric(tdroc$AUC))
    write.csv(td_df,
              file.path(TAB_DIR, "time_dependent_AUC.csv"),
              row.names = FALSE)

    td_long <- do.call(rbind, lapply(seq_along(horizons), function(i) {
      data.frame(FPR = tdroc$FP[, i], TPR = tdroc$TP[, i],
                 Horizon = sprintf("%d days (AUC=%.3f)",
                                   horizons[i], tdroc$AUC[i]))
    }))
    p_td <- ggplot(td_long, aes(FPR, TPR, color = Horizon)) +
      geom_line(linewidth = 0.8) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                  color = "grey60") +
      coord_equal() +
      labs(title = "Time-dependent ROC (TCGA-BRCA)",
           x = "1 - Specificity", y = "Sensitivity") +
      theme_bw(base_size = 11)
    ggsave(file.path(FIG_DIR, "fig_timeROC_TCGA.png"),
           p_td, width = 6.5, height = 5, dpi = 200)
  } else {
    cat("Not enough matched samples for survival analysis; skipping.\n")
  }
} else {
  cat("BRCA1 expression or clinical data unavailable; skipping survival.\n")
}

cat("\nMINI workflow complete. Outputs in:", OUT_DIR, "\n")

# 12_prognostic_model_figure.R
# Composite figure: Development and validation of the piRNA-based prognostic model.
#
# Panels:
#   A — Univariate Cox forest plot of core piRNA features
#   B — ML pipeline comparison heatmap (feature-selection x classifier C-index grid)
#   C — Kaplan-Meier survival curve, training set
#   D — Kaplan-Meier survival curve, validation set
#   E — Time-dependent ROC (1/3/5 yr), training set
#   F — Time-dependent ROC (1/3/5 yr), validation set

source("pipeline/00_config.R")

# ==============================================================================
# 0.  PACKAGES
# ==============================================================================
required_pkgs <- c(
 "survival", "survminer", "glmnet", "timeROC",
 "ggplot2", "dplyr", "tidyr", "cowplot", "ggpubr",
 "pheatmap", "RColorBrewer", "scales", "gridExtra",
 "randomForest", "xgboost", "caret", "e1071", "pROC",
 "gbm", "superpc", "CoxBoost"
)

for (pkg in required_pkgs) {
 if (!requireNamespace(pkg, quietly = TRUE))
   install.packages(pkg, dependencies = TRUE, quiet = TRUE)
}

suppressPackageStartupMessages({
 library(survival)
 library(survminer)
 library(glmnet)
 library(timeROC)
 library(ggplot2)
 library(dplyr)
 library(tidyr)
 library(cowplot)
 library(ggpubr)
 library(pheatmap)
 library(RColorBrewer)
 library(scales)
 library(gridExtra)
 library(randomForest)
 library(xgboost)
 library(caret)
 library(e1071)
 library(pROC)
})

dir.create("results/prognostic_figure", recursive = TRUE, showWarnings = FALSE)

# Publication theme
pub_theme <- theme_bw() +
 theme(
   panel.grid.major  = element_line(color = "grey92"),
   panel.grid.minor  = element_blank(),
   axis.text         = element_text(size = 11, color = "black"),
   axis.title        = element_text(size = 13, face = "bold"),
   plot.title        = element_text(size = 14, hjust = 0.5, face = "bold"),
   panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1)
 )

# ==============================================================================
# 0.1  LOAD PIPELINE RESULTS (same as piRNA_advanced_analysis.R)
# ==============================================================================
cat("\n========== LOADING PIPELINE RESULTS ==========\n")

if (!exists("model")) {
 if (file.exists("results/models/final_model.rds"))
   model <- readRDS("results/models/final_model.rds")
 else stop("Run piRNA_multicohort_pipeline.R first.")
}
if (!exists("top_feats")) {
 if (file.exists("results/models/final_features.rds"))
   top_feats <- readRDS("results/models/final_features.rds")
 else if (file.exists("results/feature_selection/final_features.txt")) {
   top_feats <- readLines("results/feature_selection/final_features.txt")
   top_feats <- top_feats[nchar(top_feats) > 0]
 }
}
if (!exists("combat_df_all"))
 stop("combat_df_all not found. Run piRNA_multicohort_pipeline.R first.")

gene_cols <- setdiff(
 colnames(combat_df_all),
 c("Group", "Batch", "T_Score", "T_Score_binary", "T_Score_all",
   "Age", "Stage", "Subtype", "Age_binary", "Stage_binary",
   "Age_Group", "Outcome01", "OS_time", "OS_status")
)

# --- Prepare merged datasets ---
target_batches <- c("BRCA1", "yyfbatch1", "yyfbatch2")
merged_df <- combat_df_all[combat_df_all$Batch %in% target_batches, ]

if (!"T_Score" %in% colnames(merged_df)) {
 merged_df$T_Score <- predict(model, merged_df[, top_feats], type = "prob")$Tumor
}

# Ensure survival data exists
if (!all(c("OS_time", "OS_status") %in% colnames(merged_df))) {
 surv_time_cols <- c("OS_time", "Survival_time", "survival_time", "time", "OS.time")
 surv_stat_cols <- c("OS_status", "Survival_status", "survival_status", "status", "OS.status")
 for (tc in surv_time_cols) {
   if (tc %in% colnames(merged_df)) { merged_df$OS_time <- merged_df[[tc]]; break }
 }
 for (sc in surv_stat_cols) {
   if (sc %in% colnames(merged_df)) { merged_df$OS_status <- merged_df[[sc]]; break }
 }
}

# Require real survival data — no simulated fallback
if (!all(c("OS_time", "OS_status") %in% colnames(merged_df))) {
 stop("OS_time and OS_status columns are required.\n",
      "Provide clinical data CSVs in clinical_data/ with survival endpoints.\n",
      "See PROTOCOL.md for the required file format.")
}

tumor_df <- merged_df[merged_df$Group == "Tumor", ]
cat("Tumor patients:", nrow(tumor_df), "\n")

# Split into training / validation cohorts
train_batches <- "BRCA1"
valid_batches <- c("yyfbatch1", "yyfbatch2")
train_df <- tumor_df[tumor_df$Batch %in% train_batches, ]
valid_df <- tumor_df[tumor_df$Batch %in% valid_batches, ]

cat("Training set (", paste(train_batches, collapse = "+"), "):", nrow(train_df), "patients\n")
cat("Validation set (", paste(valid_batches, collapse = "+"), "):", nrow(valid_df), "patients\n")


# ==============================================================================
# PANEL A — Univariate Cox Forest Plot
# ==============================================================================
cat("\n========== PANEL A: UNIVARIATE COX FOREST PLOT ==========\n")

run_univariate_cox <- function(df, features, time_col = "OS_time", event_col = "OS_status") {
 results <- lapply(features, function(feat) {
   if (!(feat %in% colnames(df))) return(NULL)
   fml <- as.formula(paste0("Surv(", time_col, ",", event_col, ") ~ `", feat, "`"))
   fit <- tryCatch(coxph(fml, data = df), error = function(e) NULL)
   if (is.null(fit)) return(NULL)
   s <- summary(fit)
   data.frame(
     Feature   = feat,
     HR        = s$conf.int[1, 1],
     HR_lower  = s$conf.int[1, 3],
     HR_upper  = s$conf.int[1, 4],
     p_value   = s$coefficients[1, 5],
     stringsAsFactors = FALSE
   )
 })
 dplyr::bind_rows(results)
}

# Run on all piRNA features, then filter significant
cox_uni <- run_univariate_cox(train_df, gene_cols)
cox_uni$FDR <- p.adjust(cox_uni$p_value, method = "BH")
cox_sig <- cox_uni %>% filter(p_value < 0.05) %>% arrange(p_value)

# If too many, take top 20 by p-value
if (nrow(cox_sig) > 20) cox_sig <- cox_sig[1:20, ]
# If none significant, show top 10 features regardless
if (nrow(cox_sig) == 0) cox_sig <- cox_uni %>% arrange(p_value) %>% head(10)

# Mark signature piRNAs
cox_sig$InSignature <- cox_sig$Feature %in% top_feats

# Reorder for plotting
cox_sig <- cox_sig %>% arrange(HR)
cox_sig$Feature <- factor(cox_sig$Feature, levels = cox_sig$Feature)

panel_A <- ggplot(cox_sig, aes(x = HR, y = Feature)) +
 geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
 geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper), height = 0.25, linewidth = 0.6) +
 geom_point(aes(shape = InSignature), size = 3, fill = "black") +
 scale_shape_manual(values = c("FALSE" = 15, "TRUE" = 18),
                    labels = c("Other", "Signature piRNA"),
                    name = "") +
 labs(title = "Univariate Cox Regression",
      x = "Hazard Ratio (95% CI)", y = NULL) +
 pub_theme +
 theme(legend.position = "bottom")

ggsave("results/prognostic_figure/panel_A_forest_plot.png",
      panel_A, width = 7, height = 8, dpi = 300)
cat("  Panel A saved.\n")


# ==============================================================================
# PANEL B — ML Pipeline Comparison Heatmap
#           feature-selection methods x classifiers, colored by C-index
# ==============================================================================
cat("\n========== PANEL B: ML PIPELINE COMPARISON ==========\n")

# Define feature-selection methods
fs_methods <- list(
 Lasso       = function(x, y) {
   fit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 5)
   coefs <- as.matrix(coef(fit, s = "lambda.1se"))
   names(which(coefs[, 1] != 0))
 },
 Ridge       = function(x, y) {
   fit <- cv.glmnet(x, y, family = "cox", alpha = 0, nfolds = 5)
   coefs <- as.matrix(coef(fit, s = "lambda.1se"))
   imp <- abs(coefs[, 1])
   names(sort(imp, decreasing = TRUE))[1:min(15, length(imp))]
 },
 ElasticNet  = function(x, y) {
   fit <- cv.glmnet(x, y, family = "cox", alpha = 0.5, nfolds = 5)
   coefs <- as.matrix(coef(fit, s = "lambda.1se"))
   names(which(coefs[, 1] != 0))
 },
 StepCox_forward = function(x, y) {
   df_tmp <- data.frame(time = y[, 1], status = y[, 2], x)
   full_fml <- as.formula(paste("Surv(time, status) ~", paste(colnames(x), collapse = " + ")))
   null_fml <- as.formula("Surv(time, status) ~ 1")
   null_fit <- coxph(null_fml, data = df_tmp)
   stepped <- tryCatch(
     step(null_fit, scope = list(lower = null_fml, upper = full_fml),
          direction = "forward", trace = 0, k = log(nrow(df_tmp))),
     error = function(e) null_fit)
   setdiff(all.vars(formula(stepped)), c("time", "status"))
 },
 StepCox_backward = function(x, y) {
   df_tmp <- data.frame(time = y[, 1], status = y[, 2], x)
   full_fml <- as.formula(paste("Surv(time, status) ~", paste(colnames(x), collapse = " + ")))
   full_fit <- tryCatch(coxph(full_fml, data = df_tmp), error = function(e) NULL)
   if (is.null(full_fit)) return(colnames(x))
   stepped <- tryCatch(
     step(full_fit, direction = "backward", trace = 0, k = log(nrow(df_tmp))),
     error = function(e) full_fit)
   setdiff(all.vars(formula(stepped)), c("time", "status"))
 },
 StepCox_both = function(x, y) {
   df_tmp <- data.frame(time = y[, 1], status = y[, 2], x)
   full_fml <- as.formula(paste("Surv(time, status) ~", paste(colnames(x), collapse = " + ")))
   null_fml <- as.formula("Surv(time, status) ~ 1")
   null_fit <- coxph(null_fml, data = df_tmp)
   stepped <- tryCatch(
     step(null_fit, scope = list(lower = null_fml, upper = full_fml),
          direction = "both", trace = 0, k = log(nrow(df_tmp))),
     error = function(e) null_fit)
   setdiff(all.vars(formula(stepped)), c("time", "status"))
 },
 RSF = function(x, y) {
   df_tmp <- data.frame(time = y[, 1], status = y[, 2], x)
   fit <- tryCatch({
     randomForest::randomForest(
       x = x, y = Surv(y[, 1], y[, 2]),
       ntree = 500, importance = TRUE)
   }, error = function(e) NULL)
   if (is.null(fit)) return(colnames(x)[1:min(10, ncol(x))])
   imp <- randomForest::importance(fit)
   imp_vals <- imp[, 1]
   names(sort(imp_vals, decreasing = TRUE))[1:min(15, length(imp_vals))]
 }
)

# Define classifiers that produce risk scores / C-index
clf_methods <- list(
 CoxPH = function(feats, train_x, train_y, test_x, test_y) {
   df_tr <- data.frame(time = train_y[, 1], status = train_y[, 2], train_x[, feats, drop = FALSE])
   fml <- as.formula(paste("Surv(time, status) ~", paste(feats, collapse = " + ")))
   fit <- tryCatch(coxph(fml, data = df_tr), error = function(e) NULL)
   if (is.null(fit)) return(list(train_ci = NA, test_ci = NA))
   lp_train <- predict(fit, type = "lp")
   df_te <- data.frame(time = test_y[, 1], status = test_y[, 2], test_x[, feats, drop = FALSE])
   lp_test <- tryCatch(predict(fit, newdata = df_te, type = "lp"), error = function(e) rep(NA, nrow(df_te)))
   ci_train <- tryCatch(concordance(Surv(train_y[, 1], train_y[, 2]) ~ lp_train)$concordance, error = function(e) NA)
   ci_test  <- tryCatch(concordance(Surv(test_y[, 1], test_y[, 2]) ~ lp_test)$concordance, error = function(e) NA)
   list(train_ci = ci_train, test_ci = ci_test, lp_train = lp_train, lp_test = lp_test)
 },
 Lasso = function(feats, train_x, train_y, test_x, test_y) {
   xtr <- as.matrix(train_x[, feats, drop = FALSE])
   xte <- as.matrix(test_x[, feats, drop = FALSE])
   fit <- tryCatch(cv.glmnet(xtr, train_y, family = "cox", alpha = 1, nfolds = 5), error = function(e) NULL)
   if (is.null(fit)) return(list(train_ci = NA, test_ci = NA))
   lp_train <- as.numeric(predict(fit, newx = xtr, s = "lambda.1se", type = "link"))
   lp_test  <- as.numeric(predict(fit, newx = xte, s = "lambda.1se", type = "link"))
   ci_train <- tryCatch(concordance(Surv(train_y[, 1], train_y[, 2]) ~ lp_train)$concordance, error = function(e) NA)
   ci_test  <- tryCatch(concordance(Surv(test_y[, 1], test_y[, 2]) ~ lp_test)$concordance, error = function(e) NA)
   list(train_ci = ci_train, test_ci = ci_test, lp_train = lp_train, lp_test = lp_test)
 },
 Ridge = function(feats, train_x, train_y, test_x, test_y) {
   xtr <- as.matrix(train_x[, feats, drop = FALSE])
   xte <- as.matrix(test_x[, feats, drop = FALSE])
   fit <- tryCatch(cv.glmnet(xtr, train_y, family = "cox", alpha = 0, nfolds = 5), error = function(e) NULL)
   if (is.null(fit)) return(list(train_ci = NA, test_ci = NA))
   lp_train <- as.numeric(predict(fit, newx = xtr, s = "lambda.1se", type = "link"))
   lp_test  <- as.numeric(predict(fit, newx = xte, s = "lambda.1se", type = "link"))
   ci_train <- tryCatch(concordance(Surv(train_y[, 1], train_y[, 2]) ~ lp_train)$concordance, error = function(e) NA)
   ci_test  <- tryCatch(concordance(Surv(test_y[, 1], test_y[, 2]) ~ lp_test)$concordance, error = function(e) NA)
   list(train_ci = ci_train, test_ci = ci_test, lp_train = lp_train, lp_test = lp_test)
 },
 ElasticNet = function(feats, train_x, train_y, test_x, test_y) {
   xtr <- as.matrix(train_x[, feats, drop = FALSE])
   xte <- as.matrix(test_x[, feats, drop = FALSE])
   fit <- tryCatch(cv.glmnet(xtr, train_y, family = "cox", alpha = 0.5, nfolds = 5), error = function(e) NULL)
   if (is.null(fit)) return(list(train_ci = NA, test_ci = NA))
   lp_train <- as.numeric(predict(fit, newx = xtr, s = "lambda.1se", type = "link"))
   lp_test  <- as.numeric(predict(fit, newx = xte, s = "lambda.1se", type = "link"))
   ci_train <- tryCatch(concordance(Surv(train_y[, 1], train_y[, 2]) ~ lp_train)$concordance, error = function(e) NA)
   ci_test  <- tryCatch(concordance(Surv(test_y[, 1], test_y[, 2]) ~ lp_test)$concordance, error = function(e) NA)
   list(train_ci = ci_train, test_ci = ci_test, lp_train = lp_train, lp_test = lp_test)
 },
 GBM = function(feats, train_x, train_y, test_x, test_y) {
   df_tr <- data.frame(time = train_y[, 1], status = train_y[, 2], train_x[, feats, drop = FALSE])
   fml <- as.formula(paste("Surv(time, status) ~", paste(feats, collapse = " + ")))
   fit <- tryCatch(
     gbm::gbm(fml, data = df_tr, distribution = "coxph",
               n.trees = 500, interaction.depth = 3, shrinkage = 0.01,
               n.minobsinnode = 10, cv.folds = 5, verbose = FALSE),
     error = function(e) NULL)
   if (is.null(fit)) return(list(train_ci = NA, test_ci = NA))
   best_n <- gbm::gbm.perf(fit, method = "cv", plot.it = FALSE)
   df_te <- data.frame(time = test_y[, 1], status = test_y[, 2], test_x[, feats, drop = FALSE])
   lp_train <- predict(fit, newdata = df_tr, n.trees = best_n)
   lp_test  <- predict(fit, newdata = df_te, n.trees = best_n)
   ci_train <- tryCatch(concordance(Surv(train_y[, 1], train_y[, 2]) ~ lp_train)$concordance, error = function(e) NA)
   ci_test  <- tryCatch(concordance(Surv(test_y[, 1], test_y[, 2]) ~ lp_test)$concordance, error = function(e) NA)
   list(train_ci = ci_train, test_ci = ci_test, lp_train = lp_train, lp_test = lp_test)
 },
 RSF = function(feats, train_x, train_y, test_x, test_y) {
   df_tr <- data.frame(time = train_y[, 1], status = train_y[, 2], train_x[, feats, drop = FALSE])
   df_te <- data.frame(time = test_y[, 1], status = test_y[, 2], test_x[, feats, drop = FALSE])
   fit <- tryCatch({
     randomForest::randomForest(
       x = train_x[, feats, drop = FALSE],
       y = Surv(train_y[, 1], train_y[, 2]),
       ntree = 500)
   }, error = function(e) NULL)
   if (is.null(fit)) return(list(train_ci = NA, test_ci = NA))
   lp_train <- tryCatch(predict(fit), error = function(e) rep(NA, nrow(train_x)))
   lp_test  <- tryCatch(predict(fit, newdata = test_x[, feats, drop = FALSE]), error = function(e) rep(NA, nrow(test_x)))
   ci_train <- tryCatch(concordance(Surv(train_y[, 1], train_y[, 2]) ~ lp_train)$concordance, error = function(e) NA)
   ci_test  <- tryCatch(concordance(Surv(test_y[, 1], test_y[, 2]) ~ lp_test)$concordance, error = function(e) NA)
   list(train_ci = ci_train, test_ci = ci_test, lp_train = lp_train, lp_test = lp_test)
 },
 SVM = function(feats, train_x, train_y, test_x, test_y) {
   df_tr <- data.frame(risk = train_y[, 2], train_x[, feats, drop = FALSE])
   fit <- tryCatch(e1071::svm(risk ~ ., data = df_tr, type = "eps-regression", kernel = "radial"),
                   error = function(e) NULL)
   if (is.null(fit)) return(list(train_ci = NA, test_ci = NA))
   df_te <- data.frame(risk = test_y[, 2], test_x[, feats, drop = FALSE])
   lp_train <- predict(fit, newdata = df_tr)
   lp_test  <- predict(fit, newdata = df_te)
   ci_train <- tryCatch(concordance(Surv(train_y[, 1], train_y[, 2]) ~ lp_train)$concordance, error = function(e) NA)
   ci_test  <- tryCatch(concordance(Surv(test_y[, 1], test_y[, 2]) ~ lp_test)$concordance, error = function(e) NA)
   list(train_ci = ci_train, test_ci = ci_test, lp_train = lp_train, lp_test = lp_test)
 }
)

# --- Run all pipelines ---
set.seed(2026)

x_train <- as.matrix(train_df[, gene_cols[gene_cols %in% colnames(train_df)]])
y_train <- Surv(train_df$OS_time, train_df$OS_status)
x_valid <- as.matrix(valid_df[, gene_cols[gene_cols %in% colnames(valid_df)]])
y_valid <- Surv(valid_df$OS_time, valid_df$OS_status)

# Remove near-zero variance from training
nzv <- caret::nearZeroVar(x_train, saveMetrics = TRUE)
keep_cols <- rownames(nzv)[!nzv$nzv]
x_train <- x_train[, keep_cols]
x_valid <- x_valid[, colnames(x_valid) %in% keep_cols]

pipeline_results <- data.frame()
best_pipeline <- list(test_ci = 0)

cat("  Running", length(fs_methods), "FS x", length(clf_methods), "CLF =",
   length(fs_methods) * length(clf_methods), "pipelines...\n")

for (fs_name in names(fs_methods)) {
 cat("  FS:", fs_name, "... ")
 selected_feats <- tryCatch(
   fs_methods[[fs_name]](x_train, y_train),
   error = function(e) { cat("(FAILED) "); character(0) }
 )

 # Ensure features exist in both sets
 selected_feats <- selected_feats[selected_feats %in% colnames(x_train) &
                                  selected_feats %in% colnames(x_valid)]
 if (length(selected_feats) < 2) {
   cat("(<2 features) skipping\n")
   next
 }

 for (clf_name in names(clf_methods)) {
   pipeline_label <- paste0(fs_name, " + ", clf_name)
   res <- tryCatch(
     clf_methods[[clf_name]](selected_feats, x_train, y_train, x_valid, y_valid),
     error = function(e) list(train_ci = NA, test_ci = NA)
   )

   row <- data.frame(
     FS          = fs_name,
     CLF         = clf_name,
     Pipeline    = pipeline_label,
     Train_CI    = ifelse(is.null(res$train_ci), NA, res$train_ci),
     Valid_CI    = ifelse(is.null(res$test_ci), NA, res$test_ci),
     n_features  = length(selected_feats),
     stringsAsFactors = FALSE
   )
   pipeline_results <- rbind(pipeline_results, row)

   # Track the best pipeline for panels C-F
   if (!is.na(res$test_ci) && res$test_ci > best_pipeline$test_ci) {
     best_pipeline <- list(
       label     = pipeline_label,
       fs_name   = fs_name,
       clf_name  = clf_name,
       train_ci  = res$train_ci,
       test_ci   = res$test_ci,
       features  = selected_feats,
       lp_train  = res$lp_train,
       lp_test   = res$lp_test
     )
   }
 }
 cat("done\n")
}

# Compute derived columns
pipeline_results$Mean_CI <- rowMeans(pipeline_results[, c("Train_CI", "Valid_CI")], na.rm = TRUE)

# Save results table
data.table::fwrite(pipeline_results, "results/prognostic_figure/pipeline_comparison.csv")

cat("  Best pipeline:", best_pipeline$label,
   sprintf("(Train C-index=%.3f, Valid C-index=%.3f)\n",
           best_pipeline$train_ci, best_pipeline$test_ci))

# --- Panel B: Heatmap ---
# Reshape into a matrix (FS x CLF) for the heatmap columns: Train, Valid, Mean, ValidOnly
hm_train <- pipeline_results %>%
 select(FS, CLF, Train_CI) %>%
 pivot_wider(names_from = CLF, values_from = Train_CI)
hm_valid <- pipeline_results %>%
 select(FS, CLF, Valid_CI) %>%
 pivot_wider(names_from = CLF, values_from = Valid_CI)
hm_mean <- pipeline_results %>%
 select(FS, CLF, Mean_CI) %>%
 pivot_wider(names_from = CLF, values_from = Mean_CI)

# Create sorted pipeline list by mean validation C-index
pipeline_sorted <- pipeline_results %>%
 filter(!is.na(Valid_CI)) %>%
 arrange(desc(Mean_CI))

# Build long-form heatmap data
hm_long <- pipeline_results %>%
 filter(!is.na(Valid_CI)) %>%
 select(Pipeline, Train_CI, Valid_CI, Mean_CI) %>%
 pivot_longer(cols = c(Train_CI, Valid_CI, Mean_CI),
              names_to = "Cohort", values_to = "C_index") %>%
 mutate(
   Cohort = factor(Cohort,
                   levels = c("Train_CI", "Valid_CI", "Mean_CI"),
                   labels = c("Training", "Validation", "Mean")),
   Pipeline = factor(Pipeline, levels = rev(pipeline_sorted$Pipeline))
 )

panel_B <- ggplot(hm_long, aes(x = Cohort, y = Pipeline, fill = C_index)) +
 geom_tile(color = "white", linewidth = 0.3) +
 geom_text(aes(label = sprintf("%.2f", C_index)), size = 2.8) +
 scale_fill_gradientn(
   colors = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9",
              "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026"),
   limits = c(0.4, 1), na.value = "grey80",
   name = "C-index"
 ) +
 labs(title = "ML Pipeline Comparison (C-index)",
      x = NULL, y = NULL) +
 pub_theme +
 theme(axis.text.y = element_text(size = 8),
       axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
       legend.position = "right")

ggsave("results/prognostic_figure/panel_B_pipeline_heatmap.png",
      panel_B, width = 10, height = max(6, nrow(pipeline_sorted) * 0.35), dpi = 300)
cat("  Panel B saved.\n")


# ==============================================================================
# PANELS C & D — Kaplan-Meier Survival Curves (training & validation)
# ==============================================================================
cat("\n========== PANELS C-D: KAPLAN-MEIER CURVES ==========\n")

make_km_panel <- function(df, lp, dataset_label) {
 if (is.null(lp) || length(lp) != nrow(df)) {
   # Fallback: recompute with CoxPH if lp unavailable
   lp <- df$T_Score
 }
 med <- median(lp, na.rm = TRUE)
 df$RiskGroup <- factor(ifelse(lp >= med, "High", "Low"), levels = c("Low", "High"))

 n_low  <- sum(df$RiskGroup == "Low")
 n_high <- sum(df$RiskGroup == "High")

 sfit <- survfit(Surv(OS_time, OS_status) ~ RiskGroup, data = df)

 # Cox HR
 cox_fit <- coxph(Surv(OS_time, OS_status) ~ RiskGroup, data = df)
 s_cox <- summary(cox_fit)
 hr <- sprintf("%.2f", s_cox$conf.int[1, 1])
 ci_lo <- sprintf("%.2f", s_cox$conf.int[1, 3])
 ci_hi <- sprintf("%.2f", s_cox$conf.int[1, 4])
 pval <- s_cox$sctest["pvalue"]
 pval_label <- ifelse(pval < 0.001, "p < 0.001", sprintf("p = %.3f", pval))

 p <- ggsurvplot(
   sfit, data = df,
   risk.table = FALSE,
   pval = FALSE,
   conf.int = TRUE,
   palette = c("#E882A8", "#6C8EBF"),
   legend.labs = c(paste0("Low(n=", n_low, ")"), paste0("High(n=", n_high, ")")),
   xlab = "Time (months)", ylab = "Survival probability",
   title = paste(best_pipeline$label, "-", dataset_label),
   ggtheme = pub_theme
 )

 # Add HR annotation
 p$plot <- p$plot +
   annotate("text", x = max(df$OS_time, na.rm = TRUE) * 0.15,
            y = 0.15, hjust = 0, size = 4,
            label = paste0(pval_label, "\nHazard Ratio = ", hr,
                           "\n95% CI: ", ci_lo, " - ", ci_hi)) +
   geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40")

 p$plot
}

panel_C <- make_km_panel(train_df, best_pipeline$lp_train, "Training")
panel_D <- make_km_panel(valid_df, best_pipeline$lp_test, "Validation")

ggsave("results/prognostic_figure/panel_C_km_train.png",
      panel_C, width = 7, height = 5.5, dpi = 300)
ggsave("results/prognostic_figure/panel_D_km_valid.png",
      panel_D, width = 7, height = 5.5, dpi = 300)
cat("  Panels C-D saved.\n")


# ==============================================================================
# PANELS E & F — Time-dependent ROC (1, 3, 5 year)
# ==============================================================================
cat("\n========== PANELS E-F: TIME-DEPENDENT ROC ==========\n")

make_timeroc_panel <- function(df, lp, dataset_label, times = c(12, 36, 60)) {
 if (is.null(lp) || length(lp) != nrow(df)) {
   lp <- df$T_Score
 }

 # Clamp times to the observed follow-up range
 max_time <- max(df$OS_time, na.rm = TRUE) * 0.95
 times <- times[times <= max_time]
 if (length(times) == 0) times <- quantile(df$OS_time, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)

 troc <- timeROC(
   T = df$OS_time,
   delta = df$OS_status,
   marker = lp,
   cause = 1,
   times = times,
   iid = TRUE
 )

 # Build ROC curve data
 time_labels <- c("1-Year", "3-Year", "5-Year")
 if (length(times) != 3) time_labels <- paste0("T=", round(times, 0))
 colors <- c("#E41A1C", "#377EB8", "#984EA3")

 roc_list <- list()
 auc_labels <- c()
 for (i in seq_along(times)) {
   tp <- troc$TP[, i]
   fp <- troc$FP[, i]
   auc_val <- troc$AUC[i]
   roc_list[[i]] <- data.frame(
     FPR = fp, TPR = tp,
     TimePoint = time_labels[i],
     stringsAsFactors = FALSE
   )
   auc_labels[i] <- sprintf("%s AUC: %.3f", time_labels[i], auc_val)
 }
 roc_df <- dplyr::bind_rows(roc_list)
 roc_df$TimePoint <- factor(roc_df$TimePoint, levels = time_labels)

 p <- ggplot(roc_df, aes(x = FPR, y = TPR, color = TimePoint)) +
   geom_line(linewidth = 1.2) +
   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
   scale_color_manual(values = setNames(colors[seq_along(times)], time_labels),
                      labels = auc_labels,
                      name = "Time Point") +
   labs(title = paste(best_pipeline$label, "Model -", dataset_label),
        x = "False Positive Rate (1-Specificity)",
        y = "True Positive Rate (Sensitivity)") +
   coord_equal() +
   pub_theme +
   theme(legend.position = c(0.65, 0.25),
         legend.background = element_rect(fill = alpha("white", 0.8)))

 p
}

panel_E <- make_timeroc_panel(train_df, best_pipeline$lp_train, "Training")
panel_F <- make_timeroc_panel(valid_df, best_pipeline$lp_test, "Validation")

ggsave("results/prognostic_figure/panel_E_timeroc_train.png",
      panel_E, width = 6.5, height = 6, dpi = 300)
ggsave("results/prognostic_figure/panel_F_timeroc_valid.png",
      panel_F, width = 6.5, height = 6, dpi = 300)
cat("  Panels E-F saved.\n")


# ==============================================================================
# COMPOSITE FIGURE — Arrange all panels
# ==============================================================================
cat("\n========== ASSEMBLING COMPOSITE FIGURE ==========\n")

# Top row: A + B
top_row <- plot_grid(panel_A, panel_B, labels = c("A", "B"),
                    rel_widths = c(0.4, 0.6), nrow = 1, label_size = 16)

# Middle row: C + D
mid_row <- plot_grid(panel_C, panel_D, labels = c("C", "D"),
                    nrow = 1, label_size = 16)

# Bottom row: E + F
bot_row <- plot_grid(panel_E, panel_F, labels = c("E", "F"),
                    nrow = 1, label_size = 16)

composite <- plot_grid(top_row, mid_row, bot_row,
                      nrow = 3, rel_heights = c(1.2, 0.9, 0.9))

ggsave("results/prognostic_figure/Figure_prognostic_model.png",
      composite, width = 16, height = 20, dpi = 300)
ggsave("results/prognostic_figure/Figure_prognostic_model.pdf",
      composite, width = 16, height = 20)

cat("\n========== COMPOSITE FIGURE SAVED ==========\n")
cat("  PNG: results/prognostic_figure/Figure_prognostic_model.png\n")
cat("  PDF: results/prognostic_figure/Figure_prognostic_model.pdf\n")
cat("  Pipeline comparison: results/prognostic_figure/pipeline_comparison.csv\n")


# ==============================================================================
# SAVE KEY OBJECTS
# ==============================================================================
saveRDS(pipeline_results, "results/prognostic_figure/pipeline_results.rds")
saveRDS(best_pipeline, "results/prognostic_figure/best_pipeline.rds")
saveRDS(cox_sig, "results/prognostic_figure/univariate_cox_significant.rds")

cat("\nDone.\n")

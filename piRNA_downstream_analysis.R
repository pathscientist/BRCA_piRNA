################################################################################
#                                                                              #
#   Breast Cancer piRNA Downstream Analysis                                    #
#                                                                              #
#   This script runs AFTER piRNA_multicohort_pipeline.R and adds:              #
#                                                                              #
#   1. Cox Regression (univariate + multivariate) with binary variables        #
#   2. Forest Plot from Cox regression                                         #
#   3. Kaplan-Meier survival curves (individual piRNAs + composite risk score) #
#   4. Subgroup ROC curves (Age, Stage, Cancer Subtype) with 95% CI           #
#   5. T-Score Boxplots by subgroup (Mann-Whitney U, exact p-values)          #
#                                                                              #
#   Requires objects from main pipeline:                                       #
#     - combat_df_all, model, top_feats, independent_sets                      #
#   OR loads from saved results if run independently                           #
#                                                                              #
################################################################################

start_time_ds <- Sys.time()

# ==============================================================================
# 0. PACKAGES
# ==============================================================================
required_pkgs <- c(
  "survival", "survminer", "pROC", "ggplot2", "dplyr", "tidyr",
  "gridExtra", "ggpubr", "caret", "randomForest", "forestmodel",
  "grid", "scales", "cowplot"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE, quiet = TRUE)
}

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(pROC)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(ggpubr)
  library(caret)
  library(randomForest)
  library(grid)
  library(scales)
  library(cowplot)
})

cat("All downstream analysis packages loaded.\n")

# ==============================================================================
# 0.1 SETTINGS & THEME
# ==============================================================================
SEED <- 2024
set.seed(SEED)

# Publication-quality theme
pub_theme <- theme_bw() +
  theme(
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(size = 11, color = "black"),
    axis.title        = element_text(size = 13, face = "bold"),
    plot.title        = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle     = element_text(size = 10, hjust = 0.5, color = "grey40"),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.text       = element_text(size = 10),
    strip.text        = element_text(size = 11, face = "bold")
  )

# Color palettes
SUBGROUP_COLORS <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
  "#FF7F00", "#A65628", "#F781BF", "#999999"
)
DIAG_COLORS <- c("Normal" = "#4393C3", "Tumor" = "#D6604D")

dir.create("results/downstream",    recursive = TRUE, showWarnings = FALSE)
dir.create("results/cox_analysis",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/km_curves",     recursive = TRUE, showWarnings = FALSE)
dir.create("results/subgroup_roc",  recursive = TRUE, showWarnings = FALSE)
dir.create("results/subgroup_box",  recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 0.2 LOAD PIPELINE RESULTS (if not already in memory)
# ==============================================================================
cat("\n========== LOADING PIPELINE RESULTS ==========\n")

# Check if main pipeline objects exist; if not, try loading from saved files
if (!exists("model") || !exists("combat_df_all")) {
  cat("Pipeline objects not in memory. Attempting to load from saved files...\n")
  if (file.exists("results/models/final_model.rds")) {
    model <- readRDS("results/models/final_model.rds")
    cat("  Loaded model from results/models/final_model.rds\n")
  } else {
    stop("Cannot find saved model. Please run piRNA_multicohort_pipeline.R first.")
  }
  if (file.exists("results/models/final_features.rds")) {
    top_feats <- readRDS("results/models/final_features.rds")
    cat("  Loaded features from results/models/final_features.rds\n")
  }
}

if (!exists("top_feats")) {
  if (file.exists("results/feature_selection/final_features.txt")) {
    top_feats <- readLines("results/feature_selection/final_features.txt")
    top_feats <- top_feats[nchar(top_feats) > 0]
    cat("  Loaded", length(top_feats), "features.\n")
  }
}

if (!exists("independent_sets")) {
  independent_sets <- c("yyfbatch1", "yyfbatch2")
  cat("  Using default independent_sets:", paste(independent_sets, collapse = ", "), "\n")
}

cat("piRNA signature:", paste(top_feats, collapse = ", "), "\n")


# ==============================================================================
# 1. PREPARE MERGED ANALYSIS DATASET
#    (BRCA1 + yyfbatch1 + yyfbatch2, with T-Score and clinical variables)
# ==============================================================================
cat("\n========== PHASE 1: PREPARE ANALYSIS DATASET ==========\n")

# Compute T-Score for all samples if not already done
if (!"T_Score" %in% colnames(combat_df_all)) {
  combat_df_all$T_Score <- predict(model, combat_df_all[, top_feats],
                                   type = "prob")$Tumor
}

# Subset to the three key datasets for downstream analysis
target_batches <- c("BRCA1", "yyfbatch1", "yyfbatch2")
merged_df <- combat_df_all[combat_df_all$Batch %in% target_batches, ]

cat("Merged dataset (BRCA1 + yyfbatch1 + yyfbatch2):", nrow(merged_df), "samples\n")
print(table(merged_df$Batch, merged_df$Group))

# --- Binarize T-Score ---
tscore_cutoff <- median(merged_df$T_Score)
merged_df$T_Score_binary <- factor(
  ifelse(merged_df$T_Score >= tscore_cutoff, "High", "Low"),
  levels = c("Low", "High")
)
cat("T-Score binary cutoff (median):", round(tscore_cutoff, 4), "\n")

# --- Clinical variables ---
# Check and create clinical variables if not present
# >>> REPLACE THIS BLOCK with your real clinical data loading <<<
if (!"Age" %in% colnames(merged_df)) {
  cat("WARNING: No clinical data found. Using simulated data for demonstration.\n")
  cat("         Replace this block with your real clinical data.\n")
  set.seed(42)
  n <- nrow(merged_df)
  merged_df$Age <- sample(30:80, n, replace = TRUE)
  merged_df$Stage <- NA
  is_tumor <- merged_df$Group == "Tumor"
  merged_df$Stage[is_tumor] <- sample(
    c("Stage I", "Stage II", "Stage III", "Stage IV"),
    sum(is_tumor), replace = TRUE, prob = c(0.35, 0.30, 0.20, 0.15))
  merged_df$Subtype <- NA
  merged_df$Subtype[is_tumor] <- sample(
    c("Luminal A", "Luminal B", "HER2+", "Triple-negative"),
    sum(is_tumor), replace = TRUE, prob = c(0.40, 0.20, 0.15, 0.25))
}

# Check for survival data
has_survival <- all(c("OS_time", "OS_status") %in% colnames(merged_df))
if (!has_survival) {
  # Try alternative column names
  surv_time_cols <- c("OS_time", "Survival_time", "survival_time", "time",
                      "OS.time", "overall_survival_time", "DFS_time")
  surv_stat_cols <- c("OS_status", "Survival_status", "survival_status", "status",
                      "OS.status", "vital_status", "event", "DFS_status")

  for (tc in surv_time_cols) {
    if (tc %in% colnames(merged_df)) {
      merged_df$OS_time <- merged_df[[tc]]
      break
    }
  }
  for (sc in surv_stat_cols) {
    if (sc %in% colnames(merged_df)) {
      merged_df$OS_status <- merged_df[[sc]]
      break
    }
  }
  has_survival <- all(c("OS_time", "OS_status") %in% colnames(merged_df))
}

if (!has_survival) {
  cat("WARNING: No survival data found in phenotype files.\n")
  cat("  Simulating survival data for demonstration.\n")
  cat("  >>> Replace with real OS_time and OS_status columns <<<\n")
  set.seed(42)
  n <- nrow(merged_df)
  is_tumor <- merged_df$Group == "Tumor"
  # Simulate survival: tumors have shorter survival; high T-score = worse prognosis
  base_time <- ifelse(is_tumor,
                      rexp(n, rate = 1/48),    # tumor: median ~48 months
                      rexp(n, rate = 1/120))   # normal: median ~120 months
  # High T-score modifies hazard for tumor patients
  hazard_modifier <- ifelse(merged_df$T_Score > tscore_cutoff & is_tumor,
                            0.7, 1.0)
  merged_df$OS_time <- pmin(round(base_time * hazard_modifier, 1), 72)
  # Censoring at 72 months
  merged_df$OS_status <- ifelse(merged_df$OS_time >= 72, 0, 1)
  merged_df$OS_time[merged_df$OS_time >= 72] <- 72
  # Additional censoring (~20%)
  censor_idx <- sample(which(merged_df$OS_status == 1),
                       round(0.2 * sum(merged_df$OS_status == 1)))
  merged_df$OS_status[censor_idx] <- 0
}

# Binarize clinical variables
merged_df$Age_binary <- factor(
  ifelse(merged_df$Age >= 60, ">=60", "<60"),
  levels = c("<60", ">=60")
)
merged_df$Age_Group <- factor(
  ifelse(merged_df$Age >= 60, "Age >= 60", "Age < 60"),
  levels = c("Age < 60", "Age >= 60")
)
merged_df$Stage_binary <- factor(
  ifelse(merged_df$Stage %in% c("Stage III", "Stage IV"), "Late (III-IV)", "Early (I-II)"),
  levels = c("Early (I-II)", "Late (III-IV)")
)

# Ensure Stage is ordered
valid_stages <- c("Stage I", "Stage II", "Stage III", "Stage IV")
if ("Stage" %in% colnames(merged_df)) {
  present_stages <- intersect(valid_stages, unique(merged_df$Stage))
  merged_df$Stage <- factor(merged_df$Stage, levels = present_stages)
}

# Binary outcome
merged_df$Outcome01 <- ifelse(merged_df$Group == "Tumor", 1, 0)

cat("\nAnalysis dataset summary:\n")
cat("  Total samples:", nrow(merged_df), "\n")
cat("  Tumor:", sum(merged_df$Group == "Tumor"),
    "Normal:", sum(merged_df$Group == "Normal"), "\n")
cat("  Median OS time:", round(median(merged_df$OS_time, na.rm = TRUE), 1), "months\n")
cat("  Events:", sum(merged_df$OS_status == 1), "\n")


# ==============================================================================
# 2. COX REGRESSION ANALYSIS (Binary Variables)
# ==============================================================================
cat("\n========== PHASE 2: COX REGRESSION ANALYSIS ==========\n")

# Work with tumor patients only for survival analysis
tumor_df <- merged_df[merged_df$Group == "Tumor", ]
cat("Tumor patients for Cox analysis:", nrow(tumor_df), "\n")

# --- Binarize piRNA expression for each feature ---
cat("\nBinarizing piRNA expression (median split within tumors)...\n")
for (feat in top_feats) {
  if (feat %in% colnames(tumor_df)) {
    med_val <- median(tumor_df[[feat]], na.rm = TRUE)
    tumor_df[[paste0(feat, "_binary")]] <- factor(
      ifelse(tumor_df[[feat]] >= med_val, "High", "Low"),
      levels = c("Low", "High")
    )
    cat(sprintf("  %s cutoff (median): %.4f  High=%d, Low=%d\n",
                feat, med_val,
                sum(tumor_df[[paste0(feat, "_binary")]] == "High"),
                sum(tumor_df[[paste0(feat, "_binary")]] == "Low")))
  }
}

# --- Re-binarize T-Score within tumor patients ---
tscore_tumor_cutoff <- median(tumor_df$T_Score)
tumor_df$T_Score_binary <- factor(
  ifelse(tumor_df$T_Score >= tscore_tumor_cutoff, "High", "Low"),
  levels = c("Low", "High")
)
cat("\nT-Score cutoff within tumors (median):", round(tscore_tumor_cutoff, 4), "\n")

# --- 2.1 Univariate Cox Regression ---
cat("\n--- 2.1 Univariate Cox Regression ---\n")

run_univariate_cox <- function(data, var_name) {
  df <- data[!is.na(data[[var_name]]), ]
  if (nrow(df) < 20 || length(unique(df[[var_name]])) < 2) return(NULL)

  fml <- as.formula(paste0("Surv(OS_time, OS_status) ~ `", var_name, "`"))
  fit <- tryCatch(coxph(fml, data = df), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  s <- summary(fit)
  coef_names <- rownames(s$coefficients)

  results <- lapply(coef_names, function(term) {
    hr    <- s$conf.int[term, "exp(coef)"]
    lower <- s$conf.int[term, "lower .95"]
    upper <- s$conf.int[term, "upper .95"]
    pval  <- s$coefficients[term, "Pr(>|z|)"]
    data.frame(
      Variable    = var_name,
      Level       = gsub(var_name, "", term),
      HR          = hr,
      Lower_95CI  = lower,
      Upper_95CI  = upper,
      P_value     = pval,
      N           = nrow(df),
      Events      = sum(df$OS_status == 1),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}

# Variables for Cox regression (all binary)
cox_vars <- c("T_Score_binary", "Age_binary", "Stage_binary")
# Add binary piRNA features
pirna_binary_vars <- paste0(top_feats, "_binary")
pirna_binary_vars <- pirna_binary_vars[pirna_binary_vars %in% colnames(tumor_df)]
cox_vars <- c(cox_vars, pirna_binary_vars)

uni_cox_results <- do.call(rbind, lapply(cox_vars, function(v) {
  run_univariate_cox(tumor_df, v)
}))

if (!is.null(uni_cox_results) && nrow(uni_cox_results) > 0) {
  uni_cox_results$Sig <- ifelse(uni_cox_results$P_value < 0.001, "***",
                         ifelse(uni_cox_results$P_value < 0.01, "**",
                         ifelse(uni_cox_results$P_value < 0.05, "*", "ns")))
  cat("\nUnivariate Cox Regression Results:\n")
  print(uni_cox_results, row.names = FALSE)
  write.csv(uni_cox_results, "results/cox_analysis/univariate_cox_results.csv",
            row.names = FALSE)
}

# --- 2.2 Multivariate Cox Regression ---
cat("\n--- 2.2 Multivariate Cox Regression ---\n")

# Select significant univariate variables (p < 0.1) for multivariate model
if (!is.null(uni_cox_results)) {
  sig_vars <- unique(uni_cox_results$Variable[uni_cox_results$P_value < 0.1])
} else {
  sig_vars <- c("T_Score_binary", "Age_binary", "Stage_binary")
}

# Always include T_Score_binary even if not significant
if (!"T_Score_binary" %in% sig_vars) {
  sig_vars <- c("T_Score_binary", sig_vars)
}

# Limit to clinical + T-score variables for multivariate (avoid collinearity)
multi_vars <- intersect(sig_vars, c("T_Score_binary", "Age_binary", "Stage_binary"))
if (length(multi_vars) < 2) {
  multi_vars <- c("T_Score_binary", "Age_binary", "Stage_binary")
}

cat("Multivariate variables:", paste(multi_vars, collapse = ", "), "\n")

df_multi_cox <- tumor_df[complete.cases(tumor_df[, c("OS_time", "OS_status", multi_vars)]), ]
fml_multi_cox <- as.formula(
  paste0("Surv(OS_time, OS_status) ~ ", paste0("`", multi_vars, "`", collapse = " + "))
)
fit_multi_cox <- coxph(fml_multi_cox, data = df_multi_cox)
s_multi_cox <- summary(fit_multi_cox)

multi_cox_results <- data.frame()
coef_names <- rownames(s_multi_cox$coefficients)
for (term in coef_names) {
  multi_cox_results <- rbind(multi_cox_results, data.frame(
    Variable    = term,
    HR          = s_multi_cox$conf.int[term, "exp(coef)"],
    Lower_95CI  = s_multi_cox$conf.int[term, "lower .95"],
    Upper_95CI  = s_multi_cox$conf.int[term, "upper .95"],
    P_value     = s_multi_cox$coefficients[term, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  ))
}
multi_cox_results$Sig <- ifelse(multi_cox_results$P_value < 0.001, "***",
                         ifelse(multi_cox_results$P_value < 0.01, "**",
                         ifelse(multi_cox_results$P_value < 0.05, "*", "ns")))

cat("\nMultivariate Cox Regression Results:\n")
print(multi_cox_results, row.names = FALSE)
write.csv(multi_cox_results, "results/cox_analysis/multivariate_cox_results.csv",
          row.names = FALSE)

# Concordance and model fit
cat(sprintf("\nModel concordance: %.3f (se=%.3f)\n",
            s_multi_cox$concordance[1], s_multi_cox$concordance[2]))
cat(sprintf("Likelihood ratio test: p=%.2e\n", s_multi_cox$logtest["pvalue"]))
cat(sprintf("Wald test: p=%.2e\n", s_multi_cox$waldtest["pvalue"]))


# ==============================================================================
# 2.3 COX REGRESSION FOREST PLOT (Publication-style table + forest)
# ==============================================================================
cat("\n--- 2.3 Cox Regression Forest Plot ---\n")

# Install forestplot if needed
if (!requireNamespace("forestplot", quietly = TRUE))
  install.packages("forestplot", dependencies = TRUE, quiet = TRUE)
library(forestplot)

# --- Prepare data for both panels ---
label_map <- c(
  "T_Score_binary"  = "risk_score",
  "Age_binary"      = "Age",
  "Stage_binary"    = "Stage"
)

# Filter univariate to clinical + T-score for forest plot
uni_for_forest <- uni_cox_results[uni_cox_results$Variable %in%
                    c("T_Score_binary", "Age_binary", "Stage_binary"), ]

# --- Helper: build one forest panel (returns a grob) ---
draw_forest_panel <- function(df, panel_title) {
  # Clean labels
  df$Names <- sapply(df$Variable, function(v) {
    base_var <- gsub("High$|>=60$|Late \\(III-IV\\)$", "", v)
    if (base_var %in% names(label_map)) label_map[base_var] else v
  })

  # Format text columns
  df$p_text  <- ifelse(df$P_value < 0.001, "<0.001",
                       sprintf("%.3f", df$P_value))
  df$hr_text <- sprintf("%.3f(%.3f,%.3f)", df$HR, df$Lower_95CI, df$Upper_95CI)

  n_rows <- nrow(df)
  tabletext <- cbind(
    c("Names",   df$Names),
    c("p.value", df$p_text),
    c("Hazard Ratio(95% CI)", df$hr_text)
  )

  mean_vals  <- c(NA, df$HR)
  lower_vals <- c(NA, df$Lower_95CI)
  upper_vals <- c(NA, df$Upper_95CI)

  # Determine sensible x-axis clip range
  all_vals <- c(df$HR, df$Lower_95CI, df$Upper_95CI)
  all_vals <- all_vals[is.finite(all_vals) & all_vals > 0]
  clip_lower <- max(min(all_vals) * 0.5, 0.01)
  clip_upper <- min(max(all_vals) * 2, 100)

  fp <- forestplot(
    labeltext   = tabletext,
    mean        = mean_vals,
    lower       = lower_vals,
    upper       = upper_vals,
    zero        = 1,
    xlog        = TRUE,
    clip        = c(clip_lower, clip_upper),
    col         = fpColors(box = "red", line = "black", zero = "gray60"),
    boxsize     = 0.25,
    lwd.zero    = 1,
    lwd.ci      = 1.5,
    ci.vertices = TRUE,
    ci.vertices.height = 0.12,
    txt_gp      = fpTxtGp(
      label  = gpar(fontfamily = "sans", cex = 0.95),
      ticks  = gpar(fontfamily = "sans", cex = 0.8),
      xlab   = gpar(fontfamily = "sans", cex = 0.9),
      title  = gpar(fontfamily = "sans", cex = 1.1, fontface = "bold")
    ),
    xlab        = "HR",
    graph.pos   = 3,
    graphwidth  = unit(5, "cm"),
    title       = panel_title,
    is.summary  = c(TRUE, rep(FALSE, n_rows)),
    hrzl_lines  = list("2" = gpar(lty = 1, lwd = 1, col = "black")),
    new_page    = FALSE
  )
  fp
}

# --- Draw both panels to a single PNG ---
png("results/cox_analysis/cox_forest_plot.png",
    width = 2400, height = 1800, res = 300)

# Split the page: top half for univariate, bottom half for multivariate
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(1, 1), "null"))))

# Top panel — Univariate
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw_forest_panel(uni_for_forest, "Univariable Cox regression")
upViewport()

# Bottom panel — Multivariate
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw_forest_panel(multi_cox_results, "Multivariable Cox regression")
upViewport()

dev.off()
cat("Cox forest plot saved: results/cox_analysis/cox_forest_plot.png\n")


# ==============================================================================
# 3. KAPLAN-MEIER SURVIVAL CURVES
# ==============================================================================
cat("\n========== PHASE 3: KAPLAN-MEIER SURVIVAL CURVES ==========\n")

# --- 3.1 KM by individual piRNA expression (binary: High vs Low) ---
cat("--- 3.1 KM curves for individual piRNAs ---\n")

km_plots <- list()
for (feat in top_feats) {
  binary_col <- paste0(feat, "_binary")
  if (!binary_col %in% colnames(tumor_df)) next

  fit_km <- survfit(Surv(OS_time, OS_status) ~ tumor_df[[binary_col]], data = tumor_df)
  # Log-rank test
  diff_test <- survdiff(Surv(OS_time, OS_status) ~ tumor_df[[binary_col]], data = tumor_df)
  pval_km <- 1 - pchisq(diff_test$chisq, df = 1)

  # Compute HR from Cox
  cox_single <- coxph(
    as.formula(paste0("Surv(OS_time, OS_status) ~ `", binary_col, "`")),
    data = tumor_df
  )
  hr_val <- exp(coef(cox_single))
  ci_cox <- exp(confint(cox_single))

  p_km <- ggsurvplot(
    fit_km,
    data = tumor_df,
    palette = c("#2166AC", "#B2182B"),  # Low=blue, High=red
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.height = 0.25,
    ggtheme = pub_theme,
    xlab = "Time (months)",
    ylab = "Overall Survival Probability",
    title = feat,
    legend.title = paste0(feat, " Expression"),
    legend.labs = c("Low", "High"),
    surv.median.line = "hv",
    pval = FALSE,  # We'll add custom annotation
    conf.int = FALSE,
    break.time.by = 12,
    xlim = c(0, 72)
  )

  # Add HR and p-value annotation
  hr_text <- sprintf("HR = %.2f (95%% CI: %.2f-%.2f)\nLog-rank p = %.2e",
                     hr_val, ci_cox[1], ci_cox[2], pval_km)
  p_km$plot <- p_km$plot +
    annotate("text", x = 50, y = 0.95, label = hr_text,
             hjust = 0, vjust = 1, size = 3.5, fontface = "italic")

  # Save individual plot
  fname <- gsub("[^a-zA-Z0-9_.-]", "_", feat)
  ggsave(paste0("results/km_curves/KM_", fname, ".png"),
         print(p_km), width = 8, height = 7, dpi = 300)

  km_plots[[feat]] <- p_km
  cat(sprintf("  %s: HR=%.2f (%.2f-%.2f), p=%.2e\n",
              feat, hr_val, ci_cox[1], ci_cox[2], pval_km))
}

# --- 3.2 KM by composite T-Score (binary: High vs Low risk) ---
cat("\n--- 3.2 KM by composite T-Score (risk score) ---\n")

fit_km_tscore <- survfit(Surv(OS_time, OS_status) ~ T_Score_binary, data = tumor_df)
diff_tscore <- survdiff(Surv(OS_time, OS_status) ~ T_Score_binary, data = tumor_df)
pval_tscore <- 1 - pchisq(diff_tscore$chisq, df = 1)

cox_tscore <- coxph(Surv(OS_time, OS_status) ~ T_Score_binary, data = tumor_df)
hr_tscore <- exp(coef(cox_tscore))
ci_tscore <- exp(confint(cox_tscore))

p_km_risk <- ggsurvplot(
  fit_km_tscore,
  data = tumor_df,
  palette = c("#2166AC", "#F4A582"),  # Low=blue, High=gold
  risk.table = TRUE,
  risk.table.col = "strata",
  risk.table.height = 0.25,
  ggtheme = pub_theme,
  xlab = "Time (months)",
  ylab = "Overall Survival Probability",
  title = "Overall Survival by Composite Risk Score (T-Score)",
  legend.title = "Risk Group",
  legend.labs = c("Low Risk (T-Score Low)", "High Risk (T-Score High)"),
  surv.median.line = "hv",
  pval = FALSE,
  conf.int = FALSE,
  break.time.by = 12,
  xlim = c(0, 72)
)

hr_risk_text <- sprintf("HR = %.2f (95%% CI: %.2f-%.2f)\nLog-rank p = %.2e",
                        hr_tscore, ci_tscore[1], ci_tscore[2], pval_tscore)
p_km_risk$plot <- p_km_risk$plot +
  annotate("text", x = 48, y = 0.95, label = hr_risk_text,
           hjust = 0, vjust = 1, size = 3.8, fontface = "italic")

ggsave("results/km_curves/KM_TScore_risk.png",
       print(p_km_risk), width = 8, height = 7, dpi = 300)

cat(sprintf("  T-Score risk: HR=%.2f (%.2f-%.2f), Log-rank p=%.2e\n",
            hr_tscore, ci_tscore[1], ci_tscore[2], pval_tscore))


# ==============================================================================
# 4. SUBGROUP ROC ANALYSIS
#    (Merged BRCA1 + yyfbatch1 + yyfbatch2, stratified by Age, Stage, Subtype)
# ==============================================================================
cat("\n========== PHASE 4: SUBGROUP ROC ANALYSIS ==========\n")

# --- Helper function: compute ROC with bootstrap 95% CI ---
roc_with_ci <- function(y_true01, y_score, n_boot = 2000, seed = 42) {
  roc_obj <- roc(y_true01, y_score, direction = "<", quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))

  set.seed(seed)
  boot_aucs <- numeric(n_boot)
  for (i in seq_len(n_boot)) {
    idx <- sample(seq_along(y_true01), replace = TRUE)
    yt <- y_true01[idx]; ys <- y_score[idx]
    if (length(unique(yt)) < 2) { boot_aucs[i] <- NA; next }
    boot_aucs[i] <- as.numeric(auc(roc(yt, ys, direction = "<", quiet = TRUE)))
  }
  boot_aucs <- na.omit(boot_aucs)
  ci <- quantile(boot_aucs, c(0.025, 0.975))

  list(roc_obj = roc_obj, auc = auc_val,
       ci_lower = ci[1], ci_upper = ci[2])
}

# --- Subgroup ROC plotting function with 95% CI ---
plot_subgroup_roc_ci <- function(df, group_col, title_text, filename) {
  df_normal <- df[df$Group == "Normal", ]
  df_tumor  <- df[df$Group == "Tumor", ]

  subgroups <- sort(unique(na.omit(df_tumor[[group_col]])))

  roc_curves <- list()
  legend_entries <- c()
  colors_used <- c()

  for (i in seq_along(subgroups)) {
    grp <- subgroups[i]
    sub_tumor <- df_tumor[!is.na(df_tumor[[group_col]]) & df_tumor[[group_col]] == grp, ]
    if (nrow(sub_tumor) < 5) next

    combined <- rbind(sub_tumor, df_normal)
    y01 <- ifelse(combined$Group == "Tumor", 1, 0)

    res <- roc_with_ci(y01, combined$T_Score, n_boot = 2000, seed = 42 + i)

    roc_df <- data.frame(
      fpr      = 1 - res$roc_obj$specificities,
      tpr      = res$roc_obj$sensitivities,
      Subgroup = grp
    )
    roc_curves[[grp]] <- roc_df

    legend_entries <- c(legend_entries,
      sprintf("%s\nAUC=%.3f (95%% CI: %.3f-%.3f)\nn(tumor)=%d, n(normal)=%d",
              grp, res$auc, res$ci_lower, res$ci_upper,
              nrow(sub_tumor), nrow(df_normal)))
    colors_used <- c(colors_used, SUBGROUP_COLORS[i])
  }

  if (length(roc_curves) == 0) {
    cat("  No valid subgroups for", group_col, "\n")
    return(NULL)
  }

  roc_all <- do.call(rbind, roc_curves)
  roc_all$Subgroup <- factor(roc_all$Subgroup, levels = subgroups)

  # Map subgroups to legend entries
  legend_map <- setNames(legend_entries, names(roc_curves))
  roc_all$Legend <- legend_map[as.character(roc_all$Subgroup)]
  roc_all$Legend <- factor(roc_all$Legend, levels = legend_entries)

  color_map <- setNames(colors_used, legend_entries)

  p <- ggplot(roc_all, aes(x = fpr, y = tpr, color = Legend)) +
    geom_path(linewidth = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60",
                linewidth = 0.5) +
    scale_color_manual(values = color_map) +
    scale_x_continuous(expand = c(0.01, 0), limits = c(0, 1.01),
                       labels = scales::number_format(accuracy = 0.1)) +
    scale_y_continuous(expand = c(0.01, 0), limits = c(0, 1.01),
                       labels = scales::number_format(accuracy = 0.1)) +
    labs(title = title_text,
         x = "1 - Specificity (False Positive Rate)",
         y = "Sensitivity (True Positive Rate)",
         color = "") +
    pub_theme +
    theme(
      legend.position = c(0.68, 0.28),
      legend.background = element_rect(fill = alpha("white", 0.9),
                                       color = "grey70", linewidth = 0.3),
      legend.text = element_text(size = 8),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.6, "cm")
    )

  ggsave(filename, p, width = 9, height = 8, dpi = 300)
  print(p)
  p
}

# --- 4.1 ROC by Age Group ---
cat("--- 4.1 ROC by Age Group ---\n")
p_roc_age <- plot_subgroup_roc_ci(merged_df, "Age_Group",
  "ROC: Discriminatory Performance by Age Group",
  "results/subgroup_roc/ROC_by_Age.png")

# --- 4.2 ROC by Cancer Stage ---
cat("--- 4.2 ROC by Cancer Stage ---\n")
p_roc_stage <- plot_subgroup_roc_ci(merged_df, "Stage",
  "ROC: Discriminatory Performance by AJCC Cancer Stage",
  "results/subgroup_roc/ROC_by_Stage.png")

# --- 4.3 ROC by Cancer Subtype ---
cat("--- 4.3 ROC by Cancer Subtype ---\n")
p_roc_subtype <- plot_subgroup_roc_ci(merged_df, "Subtype",
  "ROC: Discriminatory Performance by Histological Subtype",
  "results/subgroup_roc/ROC_by_Subtype.png")


# ==============================================================================
# 5. T-SCORE BOXPLOTS BY SUBGROUP
#    (Tumor vs Normal comparison with Mann-Whitney U, exact p-values)
# ==============================================================================
cat("\n========== PHASE 5: T-SCORE BOXPLOTS BY SUBGROUP ==========\n")

# --- Enhanced boxplot function with exact p-values and sample sizes ---
plot_subgroup_boxplot <- function(df, subgroup_col, title_text, filename) {
  df_normal <- df[df$Group == "Normal", ]
  df_tumor  <- df[df$Group == "Tumor", ]

  valid_subgroups <- sort(unique(na.omit(df_tumor[[subgroup_col]])))
  if (is.factor(df_tumor[[subgroup_col]])) {
    valid_subgroups <- levels(df_tumor[[subgroup_col]])
    valid_subgroups <- valid_subgroups[valid_subgroups %in% unique(na.omit(df_tumor[[subgroup_col]]))]
  }

  plot_data <- data.frame()
  sample_size_labels <- list()

  for (grp in valid_subgroups) {
    sub_tumor <- df_tumor[!is.na(df_tumor[[subgroup_col]]) &
                          df_tumor[[subgroup_col]] == grp, ]
    if (nrow(sub_tumor) < 3) next

    # Use all normals for comparison in each subgroup
    tmp_tumor <- sub_tumor[, c("Group", "T_Score")]
    tmp_tumor$Subgroup <- grp

    tmp_normal <- df_normal[, c("Group", "T_Score")]
    tmp_normal$Subgroup <- grp

    plot_data <- rbind(plot_data, rbind(tmp_tumor, tmp_normal))

    # Store sample sizes for annotation
    sample_size_labels[[grp]] <- data.frame(
      Subgroup = c(grp, grp),
      Group    = c("Tumor", "Normal"),
      N        = c(nrow(sub_tumor), nrow(df_normal)),
      stringsAsFactors = FALSE
    )
  }

  if (nrow(plot_data) == 0) {
    cat("  No valid subgroups for", subgroup_col, "\n")
    return(NULL)
  }

  plot_data$Subgroup <- factor(plot_data$Subgroup, levels = valid_subgroups)
  sample_sizes <- do.call(rbind, sample_size_labels)

  # Compute exact p-values (Mann-Whitney U test)
  p_values <- data.frame()
  for (grp in valid_subgroups) {
    sub_data <- plot_data[plot_data$Subgroup == grp, ]
    t_scores <- sub_data$T_Score[sub_data$Group == "Tumor"]
    n_scores <- sub_data$T_Score[sub_data$Group == "Normal"]
    if (length(t_scores) >= 3 && length(n_scores) >= 3) {
      wt <- wilcox.test(t_scores, n_scores, exact = FALSE)
      p_values <- rbind(p_values, data.frame(
        Subgroup = grp,
        p = wt$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }

  # Format p-values
  p_values$p_label <- ifelse(p_values$p < 0.001,
                             format(p_values$p, scientific = TRUE, digits = 2),
                             sprintf("p = %.4f", p_values$p))

  # Build plot
  p <- ggplot(plot_data, aes(x = Subgroup, y = T_Score, fill = Group)) +
    stat_boxplot(geom = "errorbar", width = 0.2, position = position_dodge(0.75)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85,
                 position = position_dodge(0.75),
                 color = "black", linewidth = 0.4) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
                size = 0.6, alpha = 0.25, color = "black") +
    scale_fill_manual(values = DIAG_COLORS, name = "Diagnosis") +
    scale_y_continuous(limits = c(-0.05, 1.25), breaks = seq(0, 1, 0.2),
                       expand = c(0, 0)) +
    labs(
      title = title_text,
      x = "",
      y = "T-Score"
    ) +
    pub_theme +
    theme(
      axis.text.x = element_text(size = 10, face = "bold", angle = 15, hjust = 1),
      legend.position = "top",
      panel.grid = element_blank()
    )

  # Add comparison bars with exact p-values
  y_bar <- 1.08
  for (i in seq_len(nrow(p_values))) {
    grp <- p_values$Subgroup[i]
    x_pos <- which(valid_subgroups == grp)
    p <- p +
      annotate("segment", x = x_pos - 0.2, xend = x_pos + 0.2,
               y = y_bar, yend = y_bar, linewidth = 0.4) +
      annotate("text", x = x_pos, y = y_bar + 0.03,
               label = p_values$p_label[i], size = 3, fontface = "italic")
  }

  # Add sample size labels below
  n_label_data <- sample_sizes %>%
    mutate(
      x_pos = match(Subgroup, valid_subgroups),
      x_dodge = ifelse(Group == "Normal", x_pos - 0.19, x_pos + 0.19),
      label = paste0("n=", N)
    )

  p <- p +
    annotate("text",
             x = n_label_data$x_dodge,
             y = -0.04,
             label = n_label_data$label,
             size = 2.5, color = ifelse(n_label_data$Group == "Normal",
                                        DIAG_COLORS["Normal"],
                                        DIAG_COLORS["Tumor"]),
             fontface = "bold")

  ggsave(filename, p, width = max(7, length(valid_subgroups) * 1.8 + 2),
         height = 7, dpi = 300)
  print(p)
  p
}

# --- 5.1 T-Score Boxplot by Age Group ---
cat("--- 5.1 Boxplot by Age Group ---\n")
p_box_age <- plot_subgroup_boxplot(merged_df, "Age_Group",
  "T-Score: Tumor vs Normal by Age Group",
  "results/subgroup_box/Boxplot_TScore_Age.png")

# --- 5.2 T-Score Boxplot by Cancer Stage ---
cat("--- 5.2 Boxplot by Cancer Stage ---\n")
p_box_stage <- plot_subgroup_boxplot(merged_df, "Stage",
  "T-Score: Tumor vs Normal by AJCC Cancer Stage",
  "results/subgroup_box/Boxplot_TScore_Stage.png")

# --- 5.3 T-Score Boxplot by Cancer Subtype ---
cat("--- 5.3 Boxplot by Cancer Subtype ---\n")
p_box_subtype <- plot_subgroup_boxplot(merged_df, "Subtype",
  "T-Score: Tumor vs Normal by Histological Subtype",
  "results/subgroup_box/Boxplot_TScore_Subtype.png")

# --- 5.4 T-Score Boxplot by Dataset (Batch) ---
cat("--- 5.4 Boxplot by Dataset ---\n")
p_box_batch <- plot_subgroup_boxplot(merged_df, "Batch",
  "T-Score: Tumor vs Normal by Dataset",
  "results/subgroup_box/Boxplot_TScore_Batch.png")


# ==============================================================================
# 6. COMBINED FIGURE: ROC + BOXPLOT PANEL FOR EACH STRATIFICATION FACTOR
# ==============================================================================
cat("\n========== PHASE 6: COMBINED ROC + BOXPLOT PANELS ==========\n")

# Create side-by-side panels (ROC | Boxplot) for each stratification factor
create_combined_panel <- function(roc_plot, box_plot, factor_name, filename) {
  if (is.null(roc_plot) || is.null(box_plot)) {
    cat("  Skipping combined panel for", factor_name, "(missing plot)\n")
    return(NULL)
  }

  combined <- plot_grid(
    roc_plot + theme(legend.position = "bottom",
                     legend.text = element_text(size = 7)),
    box_plot + theme(legend.position = "bottom"),
    ncol = 2, labels = c("A", "B"), label_size = 16,
    rel_widths = c(1.1, 1)
  )

  title_grob <- ggdraw() +
    draw_label(paste0("Diagnostic Performance Stratified by ", factor_name),
               fontface = "bold", size = 14, hjust = 0.5)

  final <- plot_grid(title_grob, combined, ncol = 1, rel_heights = c(0.05, 1))

  ggsave(filename, final, width = 18, height = 8, dpi = 300)
  cat("  Combined panel saved:", filename, "\n")
  final
}

create_combined_panel(p_roc_age, p_box_age, "Age",
  "results/downstream/Combined_Age.png")
create_combined_panel(p_roc_stage, p_box_stage, "Cancer Stage",
  "results/downstream/Combined_Stage.png")
create_combined_panel(p_roc_subtype, p_box_subtype, "Histological Subtype",
  "results/downstream/Combined_Subtype.png")


# ==============================================================================
# 7. SUMMARY TABLE
# ==============================================================================
cat("\n========== PHASE 7: SUMMARY ==========\n")

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  DOWNSTREAM ANALYSIS SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("piRNA Signature (", length(top_feats), "features):\n")
cat("  ", paste(top_feats, collapse = ", "), "\n\n")

cat("Datasets merged: BRCA1 + yyfbatch1 + yyfbatch2\n")
cat("  Total samples:", nrow(merged_df), "\n")
cat("  Tumor:", sum(merged_df$Group == "Tumor"),
    "  Normal:", sum(merged_df$Group == "Normal"), "\n\n")

cat("Cox Regression (multivariate):\n")
if (exists("multi_cox_results")) {
  for (i in seq_len(nrow(multi_cox_results))) {
    cat(sprintf("  %s: HR=%.2f (%.2f-%.2f), p=%s\n",
                multi_cox_results$Variable[i],
                multi_cox_results$HR[i],
                multi_cox_results$Lower_95CI[i],
                multi_cox_results$Upper_95CI[i],
                multi_cox_results$P_text[i] %||% format(multi_cox_results$P_value[i], digits = 3)))
  }
}

cat("\nOutput files:\n")
cat("  results/cox_analysis/   - Cox regression tables + forest plot\n")
cat("  results/km_curves/      - Kaplan-Meier survival curves\n")
cat("  results/subgroup_roc/   - Subgroup ROC curves with 95% CI\n")
cat("  results/subgroup_box/   - T-Score boxplots by subgroup\n")
cat("  results/downstream/     - Combined ROC + Boxplot panels\n")

end_time_ds <- Sys.time()
cat("\nDownstream analysis runtime:",
    round(difftime(end_time_ds, start_time_ds, units = "mins"), 1), "minutes\n")
cat("\n*** DOWNSTREAM ANALYSIS COMPLETE ***\n")

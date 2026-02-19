# Copy-Paste Prompt: Breast Cancer miRNA Diagnosis Pipeline (R)

> Copy everything inside the code fence below and paste it directly into Claude Code or any AI assistant.

---

```
## ROLE & CONTEXT

You are a senior bioinformatics analyst specializing in non-coding RNA biomarker
discovery and machine learning-based cancer diagnostics. You write production-grade
R code with clear comments, reproducible seeds, and publication-quality ggplot2 figures.

## MY DATA

I have a breast cancer miRNA TPM (Transcripts Per Million) expression dataset:
- Format: CSV file, loaded as data frame `mirna_data`
- Rows = samples (patients), Columns = miRNAs + one "label" column
- Label: binary factor — "Cancer" vs "Normal" (or "Tumor" vs "Normal")
- Expression values: TPM (already normalized for sequencing depth, NOT log-transformed yet)
- [EDIT: replace with your actual dimensions, e.g., "800 samples x 2500 miRNAs"]

## OBJECTIVE

Build the best possible miRNA-based diagnostic classifier for breast cancer,
suitable for independent validation on an unseen cohort. The pipeline has
three phases executed sequentially in R:

  PHASE 1 → Multi-method feature selection (8 methods) → consensus miRNA set
  PHASE 2 → Build the optimal diagnosis model with hyperparameter tuning
  PHASE 3 → Independent validation on a held-out external dataset

---

## PHASE 1: COMPREHENSIVE FEATURE SELECTION (8 METHODS)

### 1.0 Preprocessing (before feature selection)

- Log2-transform TPM values: log2(TPM + 1)
- Remove miRNAs with zero expression in >80% of samples
- Remove miRNAs with near-zero variance (caret::nearZeroVar)
- Check and report class balance; if imbalanced (ratio > 2:1), note it for later
- Split data: 70% discovery set, 30% held-out independent validation set
  * Use createDataPartition() stratified by label
  * set.seed(2024)
  * The validation set is LOCKED — never touched until Phase 3
- All 8 feature selection methods below run ONLY on the 70% discovery set

### 1.1 Differential Expression Analysis (limma)

- Use limma with empirical Bayes moderation on log2-TPM
- Design matrix: ~label
- Criteria: adjusted p-value (BH) < 0.01 AND |log2 fold change| > 1.0
- Output: volcano plot with top miRNAs labeled, sorted results table
- Save selected miRNA list as `fs_limma`

### 1.2 Wilcoxon Rank-Sum Test (non-parametric)

- Per-miRNA Wilcoxon rank-sum test (Cancer vs Normal)
- BH-adjusted p-value < 0.01 AND |median fold change| > 1.5
- This captures non-linear differences limma might miss
- Save selected miRNA list as `fs_wilcoxon`

### 1.3 Random Forest Variable Importance

- Train randomForest (ntree=1000) on discovery set
- Extract importance: MeanDecreaseGini AND MeanDecreaseAccuracy
- Select top miRNAs: importance score > mean + 1*SD (or top 50, whichever is smaller)
- Plot: horizontal barplot of top 30 miRNAs by importance
- Save selected miRNA list as `fs_rf_importance`

### 1.4 LASSO Regression (L1 Regularization)

- Use glmnet with alpha=1 (pure LASSO), family="binomial"
- cv.glmnet with 10-fold CV to find optimal lambda (lambda.min and lambda.1se)
- Select miRNAs with non-zero coefficients at lambda.1se (more parsimonious)
- Plot: coefficient path plot + CV error plot
- Save selected miRNA list as `fs_lasso`

### 1.5 Elastic Net (L1+L2 Regularization)

- Use glmnet with alpha grid search: alpha = seq(0.1, 0.9, by=0.1)
- For each alpha, run cv.glmnet → select alpha with lowest CV error
- Extract miRNAs with non-zero coefficients at best alpha + lambda.1se
- Save selected miRNA list as `fs_elasticnet`

### 1.6 Recursive Feature Elimination (RFE)

- Use caret::rfe() with rfFuncs (Random Forest-based)
- Subset sizes: c(5, 10, 15, 20, 30, 50, 75, 100)
- 10-fold CV repeated 3 times
- Plot: number of features vs CV Accuracy
- Extract optimal feature subset
- Save selected miRNA list as `fs_rfe`

### 1.7 Boruta Feature Selection

- Use Boruta package (wrapper around Random Forest)
- maxRuns = 300, doTrace = 2
- Select "Confirmed" important features (reject "Tentative" and "Rejected")
- Plot: Boruta importance boxplot (green=confirmed, yellow=tentative, red=rejected)
- Save selected miRNA list as `fs_boruta`

### 1.8 Mutual Information / mRMR (minimum Redundancy Maximum Relevance)

- Use mRMRe or praznik package to compute mutual information
- Select top features via mRMR criterion (maximizes relevance, minimizes redundancy)
- Select top 50 miRNAs (or use elbow method if available)
- Save selected miRNA list as `fs_mrmr`

### 1.9 Consensus Feature Set Construction

After running all 8 methods, determine the final miRNA set using this strategy:

a) Create a FREQUENCY TABLE: for each miRNA, count how many of the 8 methods selected it
b) Create an UpSet plot (UpSetR package) showing method overlaps — this is KEY for the paper
c) Apply TWO selection strategies and compare:

   STRATEGY A — "Strict Intersection":
   Select miRNAs found by ≥ 6 out of 8 methods (≥75% consensus)
   If this yields < 5 miRNAs, relax to ≥ 5 methods
   If this yields > 50 miRNAs, tighten to ≥ 7 methods

   STRATEGY B — "Minimum Optimal Set":
   Rank miRNAs by selection frequency (descending)
   Perform forward stepwise selection: add miRNAs one by one in frequency order,
   evaluate 5-fold CV AUC after each addition using a Random Forest
   Stop when AUC improvement < 0.005 for 3 consecutive additions
   This gives the smallest set with near-maximal performance

d) Report both sets, recommend the better one (prefer smaller set if AUC difference < 0.01)
e) Final output:
   - Table: miRNA | selection_count | which_methods | mean_importance | individual_AUC
   - Heatmap: final miRNAs x samples, clustered, annotated by Cancer/Normal (use pheatmap)
   - Venn diagram or UpSet plot of all 8 methods
   - Print: "Final miRNA signature: N miRNAs selected by M/8 methods"

---

## PHASE 2: BUILD THE BEST DIAGNOSIS MODEL

### 2.0 Data Preparation

- Use the 70% discovery set with ONLY the consensus miRNAs from Phase 1
- Within the discovery set, use 10-fold CV repeated 5 times for model training
- Preprocessing within caret pipeline: center + scale (preProcess)
- If classes are imbalanced: use SMOTE inside CV folds (via trainControl sampling="smote")
  or use class weights — try both and compare

### 2.1 Train 6 Candidate Models (all via caret::train)

Configure trainControl:
  method = "repeatedcv", number = 10, repeats = 5
  classProbs = TRUE, summaryFunction = twoClassSummary
  savePredictions = "final", search = "grid"
  metric = "ROC"

MODEL 1 — Random Forest:
  method = "rf"
  tuneGrid: mtry = seq(2, min(ncol(features), 20), by=2)
  ntree = 1000

MODEL 2 — SVM (Radial Basis Function):
  method = "svmRadial"
  tuneGrid: C = c(0.001, 0.01, 0.1, 0.5, 1, 5, 10, 50, 100)
             sigma = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.5, 1)
  preProcess = c("center", "scale")

MODEL 3 — XGBoost:
  method = "xgbTree"
  tuneGrid: nrounds = c(50, 100, 200, 500)
            max_depth = c(2, 3, 4, 6)
            eta = c(0.01, 0.05, 0.1, 0.3)
            gamma = c(0, 0.1, 0.5)
            colsample_bytree = c(0.6, 0.8, 1.0)
            min_child_weight = c(1, 3, 5)
            subsample = c(0.7, 0.8, 1.0)
  (Use a random search of ~200 combinations if full grid is too large)

MODEL 4 — Elastic Net Logistic Regression:
  method = "glmnet"
  tuneGrid: alpha = seq(0, 1, by=0.1)
            lambda = 10^seq(-4, 0, length=50)

MODEL 5 — k-Nearest Neighbors:
  method = "knn"
  tuneGrid: k = seq(1, 41, by=2)
  preProcess = c("center", "scale")

MODEL 6 — Neural Network (single hidden layer):
  method = "nnet"
  tuneGrid: size = c(1, 3, 5, 7, 10)
            decay = c(0, 0.001, 0.01, 0.1, 1)
  MaxNWts = 5000, maxit = 500

### 2.2 Model Comparison & Selection

a) Collect resamples: cv_results <- resamples(list(RF=..., SVM=..., XGB=..., GLM=..., KNN=..., NNET=...))
b) Compare by ROC-AUC, Sensitivity, Specificity
c) Statistical test: diff(cv_results) to test if differences are significant
d) Visualize:
   - Boxplot of CV AUC distributions (all 6 models side by side)
   - Dot plot with confidence intervals
   - ROC curves on discovery test fold (overlaid, all 6 models)
e) Select the BEST model based on: highest median AUC, then highest Sensitivity as tiebreaker
f) For the best model, report:
   - Optimal hyperparameters
   - CV performance: AUC, Accuracy, Sensitivity, Specificity, PPV, NPV, F1
   - Confusion matrix heatmap

### 2.3 Model Refinement (Best Model Only)

After identifying the best model:
a) Fine-tune: narrow the hyperparameter grid around the best values, re-run CV
b) Calibration: plot calibration curve (predicted probability vs observed frequency)
   If poorly calibrated, apply Platt scaling or isotonic regression
c) Determine optimal probability threshold:
   - Default is 0.5, but find the threshold that maximizes Youden's J (Sensitivity+Specificity-1)
   - Also find the threshold that maximizes F1
   - Report both, recommend Youden's J for diagnostic use
d) Learning curve: plot training set size vs AUC to check if more data would help

---

## PHASE 3: INDEPENDENT VALIDATION

### 3.0 Preparation

- Take the 30% held-out validation set (never used in Phase 1 or 2)
- Apply the EXACT same preprocessing: log2(TPM+1), center+scale using
  training set parameters (NOT re-computed on validation set)
- Subset to the consensus miRNA set from Phase 1

### 3.1 Predict & Evaluate

Using the final tuned model from Phase 2:
a) Predict classes and probabilities on the validation set
b) Apply the optimal threshold from Phase 2.3
c) Compute and report:
   - AUC (with 95% CI via DeLong method, pROC::ci.auc)
   - Accuracy, Sensitivity, Specificity, PPV, NPV, F1
   - Confusion matrix heatmap
   - ROC curve with 95% CI band
   - Precision-Recall curve with AUC
d) Compare training CV performance vs validation performance in a table
   - Flag if validation AUC drops > 0.05 (potential overfitting warning)

### 3.2 Robustness Checks

a) Bootstrap validation: resample the validation set 1000 times
   - Compute AUC for each bootstrap → report median and 95% CI
b) Permutation test: shuffle labels 1000 times → compute null AUC distribution
   - Report empirical p-value: proportion of permuted AUC ≥ observed AUC
   - Plot: histogram of null distribution with observed AUC marked
c) Per-sample prediction confidence:
   - Plot histogram of predicted probabilities, colored by true label
   - Flag samples with probability between 0.4–0.6 as "uncertain"

### 3.3 Clinical Interpretation Summary

Print a final report:
- "Diagnostic miRNA Signature: [list of N miRNAs]"
- "Best Model: [model name] with [hyperparameters]"
- "Discovery Performance (10x5 CV): AUC = X.XX (95% CI: X.XX–X.XX)"
- "Validation Performance: AUC = X.XX (95% CI: X.XX–X.XX)"
- "Optimal Threshold: X.XX (Youden's J)"
- "At optimal threshold — Sensitivity: XX.X%, Specificity: XX.X%"
- "Bootstrap 95% CI for Validation AUC: [X.XX, X.XX]"
- "Permutation test p-value: X.XXXX"

---

## GLOBAL REQUIREMENTS

1. REPRODUCIBILITY
   - set.seed(2024) before EVERY stochastic operation
   - Print R sessionInfo() at the end
   - Record start and end time, print total runtime

2. CODE QUALITY
   - One clearly commented section per step
   - Progress messages with cat() for each major step
   - All objects stored in named lists for easy access
   - Wrap repeated operations in helper functions

3. VISUALIZATIONS (ggplot2, publication-quality)
   - theme_minimal(base_size=13) with bold centered titles
   - Color palette: c("#E31A1C","#1F78B4","#33A02C","#FF7F00","#6A3D9A","#B15928")
   - All plots saved to PNG (300 DPI, 8x6 inches) AND displayed on screen
   - Key plots: volcano, UpSet, heatmap, ROC, confusion matrix, CV boxplot,
     calibration curve, bootstrap histogram, prediction probability distribution

4. OUTPUT FILES
   Create output directory structure:
   results/
   ├── feature_selection/
   │   ├── fs_summary_table.csv        (all miRNAs, selection counts, importance)
   │   ├── consensus_mirnas.txt         (final miRNA list, one per line)
   │   ├── volcano_plot.png
   │   ├── upset_plot.png
   │   └── heatmap_consensus.png
   ├── models/
   │   ├── model_comparison.csv         (all 6 models performance)
   │   ├── best_model.rds               (saved final model object)
   │   ├── preprocessing_params.rds     (center/scale values)
   │   ├── cv_boxplot.png
   │   └── roc_all_models.png
   └── validation/
       ├── validation_results.csv
       ├── roc_validation.png
       ├── confusion_matrix.png
       ├── bootstrap_ci.png
       ├── permutation_test.png
       └── final_report.txt

5. PACKAGES (install check at top of script)
   Required: caret, randomForest, e1071, glmnet, xgboost, nnet, pROC,
             limma, Boruta, UpSetR, pheatmap, ggplot2, dplyr, tidyr,
             gridExtra, viridis, scales, smotefamily (or DMwR2), praznik (or mRMRe)

6. ERROR HANDLING
   - Wrap each feature selection method in tryCatch()
   - If any method fails, warn but continue with remaining methods
   - If consensus set is empty, fall back to LASSO-only features + top 20 RF importance

Please generate the complete R script now. Structure it clearly with section headers.
Prioritize correctness and biological interpretability over speed.
```

---

## How to Customize This Prompt

Before pasting, **search-and-replace** these placeholders with your actual values:

| Placeholder | Replace with | Example |
|---|---|---|
| `"Cancer" vs "Normal"` | Your actual label names | `"Tumor" vs "Adjacent"` |
| `800 samples x 2500 miRNAs` | Your actual dimensions | `150 samples x 1800 miRNAs` |
| `set.seed(2024)` | Any seed you prefer | `set.seed(42)` |
| `adj.p < 0.01` | Adjust if too strict/lenient | `adj.p < 0.05` |
| `log2FC > 1.0` | Adjust fold change cutoff | `log2FC > 0.5` |
| `≥ 6 out of 8` | Adjust consensus threshold | `≥ 5 out of 8` |

## Tips for Using This Prompt Effectively

1. **If the script is too long**, ask Claude to generate Phase 1, 2, and 3 separately
2. **If a method fails**, paste the error and say: "Method X.X failed with this error: [error]. Fix it."
3. **If you want to add a 9th feature selection method**, say: "Add Information Gain (FSelector package) as method 1.9 before the consensus step"
4. **If you have a separate validation cohort file**, replace Phase 3.0 with: "Load external validation data from `validation_cohort.csv` instead of using the 30% split"
5. **If runtime is too long**, say: "XGBoost grid is too large, reduce to random search with 50 combinations"

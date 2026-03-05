# Prompt Templates for Breast Cancer piRNA RNA-seq Analysis & ML Diagnosis in R

Copy-paste any of these prompts into Claude Code (or any AI assistant) to get help with your project.

---

## 1. Full Project Setup Prompt

```
I'm building a breast cancer diagnosis model using piRNA (Piwi-interacting RNA)
expression data from RNA-seq in R. The project involves:

- Data: piRNA expression matrix (rows = patient samples, columns = piRNA features),
  with a binary "label" column (e.g., "Cancer" vs "Normal", or "Tumor" vs "Healthy")
- Goal: classify breast cancer patients vs healthy controls using machine learning
- Language: R
- Key libraries: caret, randomForest, e1071, glmnet, xgboost, pROC, ggplot2

Please help me build a complete analysis pipeline that includes:

1. DATA LOADING & EXPLORATION
   - Load the piRNA expression matrix (CSV/TSV)
   - Summarize dimensions, class distribution, missing values
   - Visualize label distribution (barplot)

2. PREPROCESSING
   - Handle missing values (imputation or removal)
   - Normalize/transform expression data (log2, TMM, or variance stabilization)
   - Filter low-variance piRNAs
   - Check and address class imbalance (SMOTE or weighted models if needed)

3. FEATURE SELECTION
   - Differential expression analysis (limma or DESeq2-style)
   - Variance-based filtering
   - Recursive Feature Elimination (RFE) with caret
   - Correlation-based filtering to remove redundant features
   - Select top N informative piRNAs

4. MODEL TRAINING & COMPARISON (10-fold CV, 80/20 train/test split)
   - Random Forest (tuning: mtry, ntree)
   - Support Vector Machine with radial kernel (tuning: C, sigma)
   - Elastic Net / Logistic Regression (tuning: alpha, lambda)
   - XGBoost (tuning: nrounds, max_depth, eta)
   - Naive Bayes (baseline)

5. EVALUATION & VISUALIZATION
   - Confusion matrices (heatmaps)
   - ROC curves with AUC for all models on one plot
   - Precision-Recall curves (important if classes are imbalanced)
   - Model performance comparison barplot (Accuracy, Sensitivity, Specificity, AUC)
   - Cross-validation boxplots showing score distributions
   - Feature importance plot (top 20 piRNAs from best model)

6. RESULTS SUMMARY
   - Ranked model comparison table
   - Best model recommendation with justification
   - List of top biomarker piRNA candidates

Use set.seed(123) for reproducibility. Use ggplot2 with a clean minimal theme
for all plots. Store all models and results in named lists for easy comparison.
```

---

## 2. Data Preprocessing Only

```
I have a piRNA expression matrix for breast cancer RNA-seq analysis in R.
The data frame is called `pirna_data` with:
- Rows: patient samples
- Columns: piRNA expression values + a "label" column (Cancer/Normal)

Please write R code to:
1. Check dimensions, data types, and missing value summary
2. Log2-transform the expression values (add pseudocount of 1)
3. Remove piRNAs with near-zero variance using caret::nearZeroVar()
4. Normalize samples (quantile normalization or scale/center)
5. Check class balance and report if SMOTE/oversampling is needed
6. Produce a PCA plot colored by label to check sample separation
7. Produce a correlation heatmap of the top 50 most variable piRNAs

Use ggplot2 for all visualizations. Comment each step clearly.
```

---

## 3. Feature Selection / Differential Expression

```
I have a preprocessed piRNA expression matrix (`pirna_data`) for breast cancer
classification in R. The "label" column has two classes.

Please write R code to perform feature selection:

1. DIFFERENTIAL EXPRESSION
   - Use limma (with voom if counts) or Wilcoxon rank-sum test
   - Compute fold change and adjusted p-values for each piRNA
   - Volcano plot with labeled top hits
   - Filter significant piRNAs (adj.p < 0.05, |log2FC| > 1)

2. MACHINE LEARNING-BASED SELECTION
   - Recursive Feature Elimination (RFE) using caret with Random Forest
   - Plot number of features vs CV accuracy
   - Extract optimal feature subset

3. COMBINE & FINALIZE
   - Intersect DEG-based and RFE-based selected piRNAs
   - Report final feature set with a summary table
   - Heatmap of selected piRNAs (samples as columns, piRNAs as rows, annotated by label)

Use set.seed(123). Visualize with ggplot2 and pheatmap/ComplexHeatmap.
```

---

## 4. Model Training & Hyperparameter Tuning

```
I have a breast cancer piRNA expression dataset in R, already split into
train_x, train_y, test_x, test_y. Features are pre-selected piRNAs.

Please write R code to train and tune these models using caret with
10-fold cross-validation, optimizing for ROC:

1. Random Forest - tune mtry (2 to sqrt(p) in steps)
2. SVM (radial) - tune C (0.01, 0.1, 1, 10, 100) and sigma (0.001 to 1)
3. Elastic Net - tune alpha (0 to 1) and lambda (grid of 20 values)
4. XGBoost - tune nrounds, max_depth, eta, subsample
5. k-Nearest Neighbors - tune k (3 to 25, odd values)

For each model:
- Use trainControl with method="repeatedcv", number=10, repeats=3
- Enable classProbs=TRUE and summaryFunction=twoClassSummary
- Store predictions and probabilities on the test set
- Compute confusion matrix

After training all models, create:
- A summary data frame with Accuracy, Sensitivity, Specificity, AUC, F1 per model
- Identify and report the best model

Use set.seed(456) before each model. Print progress messages.
```

---

## 5. Visualization Suite

```
I have trained 5 ML models for breast cancer piRNA classification in R.
I have these objects ready:
- model_results: named list of confusionMatrix objects
- roc_data: named list of pROC roc objects
- rf_model: trained caret Random Forest (with variable importance)
- cv_results: caret resamples() object

Please write R code to generate publication-quality figures using ggplot2:

1. ROC CURVES - all models on one plot, with AUC in legend, diagonal reference line
2. PERFORMANCE BARPLOT - grouped bars for Accuracy/Sensitivity/Specificity/AUC
3. CONFUSION MATRIX HEATMAPS - one per model, arranged in a grid
4. FEATURE IMPORTANCE - horizontal bar chart of top 20 piRNAs from Random Forest
5. CV BOXPLOTS - boxplot of 10-fold CV AUC scores per model with jittered points
6. CALIBRATION CURVES - predicted probability vs observed frequency per model
7. DECISION CURVE ANALYSIS - net benefit across threshold probabilities

Settings:
- Use theme_minimal() with centered bold titles
- Color palette: viridis or custom (steelblue, tomato, forestgreen, darkorange, purple)
- Export-ready: ggsave() each plot at 300 DPI, 8x6 inches
- Arrange related plots with patchwork or gridExtra

Output each plot to both screen and PDF/PNG file.
```

---

## 6. Biomarker Discovery & Interpretation

```
I have a trained Random Forest and Elastic Net model for breast cancer piRNA
classification in R. I want to identify candidate piRNA biomarkers.

Please write R code to:

1. FEATURE IMPORTANCE ANALYSIS
   - Extract and rank variable importance from Random Forest (MeanDecreaseGini + MeanDecreaseAccuracy)
   - Extract non-zero coefficients from Elastic Net
   - Find consensus top piRNAs (present in both model rankings)
   - Create a combined importance score and ranked table

2. INDIVIDUAL piRNA ANALYSIS (for top 10 piRNAs)
   - Boxplot: expression in Cancer vs Normal (with p-values, Wilcoxon test)
   - Individual ROC curve with AUC for each piRNA as a standalone biomarker
   - Summary table: piRNA name, importance score, individual AUC, fold change, p-value

3. MULTI-MARKER PANEL
   - Train a logistic regression model using only the top 5 consensus piRNAs
   - Compare its ROC/AUC to the full model
   - Report if a small panel is sufficient for diagnosis

4. EXPORT RESULTS
   - Save biomarker table as CSV
   - Save all plots to a PDF report
```

---

## 7. Class Imbalance Handling

```
My breast cancer piRNA dataset is imbalanced (e.g., 70% Cancer, 30% Normal).
Please write R code in to address this:

1. Report current class ratio
2. Apply these strategies and compare model performance:
   a. SMOTE oversampling (using DMwR2 or smotefamily package)
   b. Random undersampling of majority class
   c. Class weights in caret (weights parameter or classwt for RF)
   d. Cost-sensitive learning (adjusted cutoff thresholds)
3. Train Random Forest under each strategy with 10-fold CV
4. Compare: ROC-AUC, PR-AUC, Sensitivity, Specificity, F1
5. Visualize with grouped barplot and PR curves
6. Recommend the best strategy for this dataset

Use Precision-Recall AUC as the primary metric (more informative than ROC for imbalanced data).
```

---

## 8. External Validation / New Sample Prediction

```
I have a trained and validated best model (saved as .rds) for breast cancer
piRNA classification. Now I want to:

1. SAVE THE FINAL MODEL
   - Save the best model object with saveRDS()
   - Save the preprocessing parameters (scaling center/scale values)
   - Save the selected feature list
   - Create a prediction wrapper function

2. LOAD & PREDICT ON NEW SAMPLES
   - Write a function predict_cancer(new_data, model_path) that:
     a. Loads the saved model
     b. Checks that required piRNA columns exist
     c. Applies the same preprocessing (centering, scaling)
     d. Returns: predicted class, probability of Cancer, confidence level
   - Handle edge cases (missing piRNAs, out-of-range values)

3. VALIDATION ON INDEPENDENT DATASET
   - Load an external validation cohort
   - Apply the prediction function
   - Report performance metrics and ROC curve
   - Compare training vs validation performance

Write clean, well-documented R code with error handling.
```

---

## 9. Complete Reproducible Pipeline (Single Script)

```
Please write a single, complete, reproducible R script for breast cancer
diagnosis using piRNA expression data. The script should run end-to-end
with minimal user modification (only the input file path).

Requirements:
- Input: CSV file with samples as rows, piRNAs as columns, "label" column
- Automatic package installation check (install if missing)
- Sections clearly separated with comments and cat() progress messages
- set.seed() at every stochastic step
- All results saved to an output directory (plots/, models/, tables/)
- Models: Random Forest, SVM, Elastic Net, XGBoost
- Evaluation: ROC, PR curves, confusion matrices, CV boxplots
- Feature importance and top piRNA biomarker table
- Final summary printed to console
- Total runtime printed at end

Structure the script as:
  0. Setup & package loading
  1. Data loading & exploration
  2. Preprocessing & QC
  3. Feature selection
  4. Train/test split
  5. Model training with CV
  6. Evaluation & visualization
  7. Biomarker analysis
  8. Save outputs
  9. Summary report

Target: publication-quality analysis suitable for a bioinformatics manuscript.
```

---

## 10. Quick Fix / Debug Prompts

```
# When something breaks, paste the error with context:

I'm running a piRNA breast cancer classification pipeline in R.
I got this error: [PASTE ERROR HERE]

The relevant code section is:
[PASTE CODE BLOCK]

My data looks like:
- Dimensions: [e.g., 200 samples x 1500 piRNAs]
- Label column: [e.g., factor with levels "Cancer" and "Normal"]
- R version: [e.g., 4.3.1]
- OS: [e.g., Ubuntu 22.04]

Please explain what's wrong and provide a fix.
```

---

## Tips for Best Results

1. **Be specific about your data**: Always mention dimensions, column names, class labels
2. **Share errors completely**: Include the full error message and traceback
3. **State your R version**: Some packages behave differently across versions
4. **Mention installed packages**: Avoids suggestions for packages you don't have
5. **Specify plot output**: Say if you need PNG, PDF, or screen display
6. **Ask for one thing at a time**: Break large requests into focused prompts
7. **Request reproducibility**: Always ask for `set.seed()` in stochastic steps

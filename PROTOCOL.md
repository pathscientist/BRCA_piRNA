# piRNA Breast Cancer Diagnostic Pipeline — Execution Protocol

## Overview

This pipeline builds a piRNA-based diagnostic classifier for breast cancer using multi-cohort data with ComBat batch correction, multi-method feature selection, and independent plasma validation.

**Data:**
| Cohort | Specimen | Role | Samples |
|--------|----------|------|---------|
| BRCA1 | Tissue (TCGA) | Training | ~194 (187 Tumor, 7 Normal) |
| yyfbatch1 | Plasma (in-house) | Training | 192 (106 Tumor, 86 Normal) |
| yyfbatch2 | Plasma (in-house) | Independent Validation | 212 (119 Tumor, 93 Normal) |

**Note:** Both yyfbatch1 and yyfbatch2 are **plasma** specimens. yyfbatch1 is included in training because BRCA1 alone has only 7 Normal samples, far too few for the model to learn Normal vs Tumor.

**Target:** ~10 piRNA signature (7–12 range).

---

## Prerequisites

- **R** version 4.0 or higher
- **Internet connection** for first run (auto-installs missing packages)
- Input files in `processed_results/`:
  - `BRCA1_processed_1.csv`
  - `yyfbatch1_processed.csv`
  - `yyfbatch2_processed.csv`

---

## Execution Order

Run the 5 scripts **sequentially** from the project root directory. Each script saves `.rds` files that the next one loads.

```
cd BRCA_piRNA/
```

### Script 01 — Data Loading & Batch Correction

```bash
Rscript 01_data_loading_batch_correction.R
```

**What it does:**
- Loads all cohort CSVs, harmonizes group labels (Tumor/Normal)
- Keeps all BRCA1 samples (no downsampling — only 7 normals)
- Finds common piRNAs across all cohorts
- Applies log2(TPM+1) transformation
- Runs ComBat batch correction on ALL cohorts together (BRCA1 + yyfbatch1 + yyfbatch2)
- Z-score normalizes
- Generates PCA plots (before/after ComBat)
- Splits: training (BRCA1 + yyfbatch1) vs validation (yyfbatch2)

**Check before proceeding:**
- Open `results/qc/pca_combat.png` — after ComBat, cohorts should overlap while Tumor/Normal still separate
- Console shows sample counts per cohort and class distribution

**Key outputs:**
```
results/models/combat_df_all.rds       # All cohorts combined
results/models/train_df.rds            # BRCA1 + yyfbatch1
results/models/valid_yyfbatch2.rds     # Independent hold-out
results/models/preprocessing_params.rds
results/qc/pca_combat.png
```

---

### Script 02 — Feature Selection (~10 piRNAs)

```bash
Rscript 02_feature_selection.R
```

**What it does:**
- Runs 8 feature selection methods on training data ONLY:
  1. limma (empirical Bayes DE)
  2. Wilcoxon rank-sum
  3. Random Forest importance
  4. LASSO (glmnet, alpha=1)
  5. Elastic Net (best alpha by CV)
  6. RFE (Recursive Feature Elimination)
  7. Boruta
  8. mRMR (praznik)
- Builds consensus via frequency threshold + forward stepwise
- Refines with forward/backward/swap optimization
- Final signature forced to 7–12 piRNAs
- Generates justification table, UpSet plot, volcano plot, incremental AUC curve, correlation heatmap

**Check before proceeding:**
- Console prints final piRNA signature (should be 7–12 piRNAs)
- Review `results/feature_selection/incremental_auc.png` — AUC should plateau around 8–10 features
- Review `results/feature_selection/fs_summary_table.csv` — each piRNA should have individual AUC > 0.55

**Key outputs:**
```
results/models/final_features.rds
results/feature_selection/consensus_pirnas.txt
results/feature_selection/fs_summary_table.csv
results/feature_selection/volcano_plot.png
results/feature_selection/rf_importance.png
results/feature_selection/upset_plot.png
results/feature_selection/incremental_auc.png
results/feature_selection/signature_correlation.png
```

---

### Script 03 — Model Training (6 Classifiers)

```bash
Rscript 03_model_training.R
```

**What it does:**
- Trains 6 classifiers using the consensus piRNAs:
  1. Random Forest
  2. SVM-RBF
  3. XGBoost
  4. Elastic Net (glmnet)
  5. KNN
  6. Neural Net (nnet)
- 10-fold CV × 5 repeats for each, with downsampling if class imbalance > 2:1
- Compares models by median AUC, selects best
- Fine-tunes best model (narrow grid around best hyperparameters)
- Determines optimal threshold: Youden's J + sensitivity ≥ 95% (screening)
- Generates calibration curve and learning curve

**Check before proceeding:**
- Console prints "Best model: [name], CV AUC = X.XX (95% CI: ...)"
- CV AUC should be > 0.85
- Review `results/models/cv_boxplot.png` for model comparison
- Review `results/models/calibration.png` for calibration

**Key outputs:**
```
results/models/final_model.rds
results/models/optimal_threshold.rds
results/models/all_candidate_models.rds
results/models/model_comparison_table.csv
results/models/cv_boxplot.png
results/models/roc_all_models.png
results/models/calibration.png
results/models/learning_curve.png
```

---

### Script 04 — Independent Validation

```bash
Rscript 04_independent_validation.R
```

**What it does:**
- Evaluates the final model on yyfbatch2 (plasma, NEVER seen during training)
- Computes AUC (DeLong 95% CI), sensitivity, specificity, PPV, NPV, F1
- Generates 3-panel ROC figure: Discovery CV / Validation / Overlay
- Per-batch ROC overlay across all cohorts
- Bootstrap AUC (1000 resamples) with 95% CI
- Permutation test (1000 label shuffles) with empirical p-value
- Performance drop assessment: Delta(train − validation)

**Check before proceeding:**
- Console prints validation AUC and confusion matrix stats
- Review `results/validation/Figure_ROC_validation.png`
- Check `results/validation/yyfbatch2_confusion.png`
- Delta(train–validation) ideally < 0.10

**Key outputs:**
```
results/validation/validation_results.csv
results/validation/Figure_ROC_validation.png / .pdf
results/validation/yyfbatch2_confusion.png
results/validation/yyfbatch2_roc.png
results/validation/yyfbatch2_prc.png
results/validation/yyfbatch2_probability_hist.png
results/validation/yyfbatch2_bootstrap.png
results/validation/yyfbatch2_permutation.png
results/validation/auc_comparison_bar.png
results/validation/roc_all_cohorts.png
results/validation/validation_details.rds
```

---

### Script 05 — Clinical Interpretation & Final Report

```bash
Rscript 05_clinical_interpretation_report.R
```

**What it does:**
- Generates signature characterization table (for paper Table 1): piRNA ID, N methods, direction, log2FC, adj.p, individual AUC, RF importance
- Creates summary heatmap across all cohorts (annotated by Group, Cohort, Specimen type)
- Produces final text report with all metrics
- Saves session info

**Key outputs:**
```
results/final_signature_table.csv
results/feature_selection/heatmap_all_cohorts.png
results/validation/final_report.txt
results/session_info.txt
```

---

### Script 06 — Function Prediction Figure (Pearson + Reactome)

```bash
Rscript 06_function_prediction_figure.R
```

**What it does:**
- Focuses on three piRNAs: `piR-hsa-41032`, `piR-hsa-1348371`, `piR-hsa-128633`
- Matches the TCGA piRNA cohort (`BRCA1_processed_1.csv`) to TCGA mRNA expression (`mRNA_expression/TCGA_BRCA_mRNA.csv`) by sample barcode
- Computes Pearson correlation between each piRNA and all genes; keeps top correlated genes (|r| > 0.3, FDR < 0.05)
- Runs Reactome pathway enrichment (`ReactomePA`) on the correlated genes
- Builds the two-panel figure:
  - **Panel A** — dot-heatmap of genes correlated with the piRNAs (color = Pearson r, size = −log10 p)
  - **Panel B** — piRNA → gene → Reactome pathway functional network
- Falls back to a curated gene set for any piRNA that is ~0 in the matched samples (e.g. low-abundance `piR-hsa-128633`); such rows are flagged `source = curated`

**Requirements:** `mRNA_expression/TCGA_BRCA_mRNA.csv` must be present. Standalone — does not depend on scripts 01–05.

**Key outputs:**
```
results/functional/Figure_function_prediction.png / .pdf
results/functional/piRNA_gene_correlations.csv
results/functional/piRNA_gene_pathway_edges.csv
results/functional/reactome_enrichment.csv
```

---

## Run All at Once

To run the entire pipeline in one go:

```bash
cd BRCA_piRNA/

Rscript 01_data_loading_batch_correction.R && \
Rscript 02_feature_selection.R && \
Rscript 03_model_training.R && \
Rscript 04_independent_validation.R && \
Rscript 05_clinical_interpretation_report.R
```

Or from within an R session:

```r
setwd("path/to/BRCA_piRNA")

source("01_data_loading_batch_correction.R")
source("02_feature_selection.R")
source("03_model_training.R")
source("04_independent_validation.R")
source("05_clinical_interpretation_report.R")
```

**Estimated total runtime:** 30–90 minutes depending on hardware (Step 03 is the longest due to 6 × 10×5 CV models).

---

## Output Directory Structure

```
results/
├── qc/
│   └── pca_combat.png                      # PCA before/after ComBat
├── feature_selection/
│   ├── volcano_plot.png                    # limma DE volcano
│   ├── rf_importance.png                   # RF importance barplot
│   ├── upset_plot.png                      # 8-method overlap
│   ├── incremental_auc.png                 # AUC vs number of piRNAs
│   ├── signature_correlation.png           # Pairwise correlation heatmap
│   ├── heatmap_all_cohorts.png             # Expression heatmap, all samples
│   ├── consensus_pirnas.txt                # Final piRNA list (one per line)
│   ├── fs_summary_table.csv                # Justification table
│   └── fs_frequency_table.csv              # Selection frequency
├── models/
│   ├── combat_df_all.rds                   # All cohorts, batch-corrected
│   ├── train_df.rds                        # Training data (BRCA1 + yyfbatch1)
│   ├── valid_yyfbatch2.rds                 # Validation hold-out
│   ├── final_model.rds                     # Trained classifier
│   ├── final_features.rds                  # Consensus piRNA names
│   ├── optimal_threshold.rds               # Youden's J threshold
│   ├── preprocessing_params.rds            # Z-score parameters
│   ├── all_candidate_models.rds            # All 6 trained models
│   ├── model_comparison_table.csv          # 6-model comparison
│   ├── cv_boxplot.png                      # CV AUC boxplot
│   ├── roc_all_models.png                  # Overlaid ROC, 6 models
│   ├── calibration.png                     # Calibration curve
│   └── learning_curve.png                  # AUC vs training size
├── validation/
│   ├── validation_results.csv              # Summary metrics table
│   ├── Figure_ROC_validation.png / .pdf    # Key multi-panel ROC figure
│   ├── yyfbatch2_confusion.png             # Confusion matrix
│   ├── yyfbatch2_roc.png                   # Individual ROC
│   ├── yyfbatch2_prc.png                   # Precision-Recall
│   ├── yyfbatch2_probability_hist.png      # Prediction probability histogram
│   ├── yyfbatch2_bootstrap.png             # Bootstrap AUC distribution
│   ├── yyfbatch2_permutation.png           # Permutation null distribution
│   ├── auc_comparison_bar.png              # AUC bar chart with CI
│   ├── roc_all_cohorts.png                 # Per-batch ROC overlay
│   ├── validation_details.rds              # Detailed results for script 05
│   └── final_report.txt                    # Final text report
├── final_signature_table.csv               # Signature characterization (Table 1)
└── session_info.txt                        # R session info
```

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| Package install fails | Run `install.packages("pkgname")` or `BiocManager::install("pkgname")` manually in R |
| `Error: No cohort files found` | Place CSV files in `processed_results/` directory |
| ComBat error | Check that all CSV files have a `Group` column and share common piRNA names |
| PCA shows no cohort separation after ComBat | This is correct — cohorts SHOULD overlap. Check that Tumor/Normal still separate |
| Feature selection finds < 7 piRNAs | Relax thresholds in script 02: `adj.p < 0.05` and `|log2FC| > 0.5` |
| CV AUC < 0.85 | Check PCA plot; may need more training data or relaxed feature count |
| Validation AUC much lower than CV AUC | Check PCA for residual batch effects; inspect probability histogram |
| Script fails mid-run | Re-run from the failed script only — all previous `.rds` files are already saved |
| `Error: file not found ... .rds` | Run the previous script first — each script depends on outputs from the one before it |

---

## Downstream Analysis (Optional)

After completing scripts 01–05, the following existing scripts can be run for additional analyses. They load objects from `results/models/`:

| Script | What it does | Extra data needed |
|--------|-------------|-------------------|
| `piRNA_downstream_analysis.R` | Cox regression, KM curves, subgroup ROC | Clinical data with OS_time/OS_status |
| `piRNA_functional_meta_analysis.R` | Expression heatmap, piRNA-mRNA correlation, pathway enrichment | mRNA expression matrix |
| `piRNA_network_analysis.R` | 8 network visualizations | None |
| `piRNA_advanced_analysis.R` | Nomogram, SHAP, immune, TMB, drug sensitivity | mRNA + clinical + TMB data |

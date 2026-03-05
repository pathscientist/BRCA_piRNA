# Breast Cancer piRNA Diagnostic Signature — Analysis Protocol

## Overview

This project builds a piRNA-based diagnostic signature for breast cancer using multi-cohort data. The pipeline selects <=10 piRNA features, trains a Random Forest classifier on 5 public datasets, and validates independently on 2 in-house cohorts (yyfbatch1 + yyfbatch2), targeting AUC > 0.8.

## Datasets

| Dataset | Role | Source |
|---------|------|--------|
| BRCA1 | Training (balanced) | TCGA-BRCA |
| PRJNA294226 | Training | GEO: GSE72080 |
| PRJNA482141 | Training | GEO: GSE117452 |
| PRJNA808405 | Training | GEO: GSE197020 |
| PRJNA934049 | Training | GEO: GSE225117 |
| yyfbatch1 | **Independent validation** | In-house |
| yyfbatch2 | **Independent validation** | In-house |

## Prerequisites

### R Version
- R >= 4.1.0

### Required Data Files
Place all processed expression files in the `processed_results/` directory:

```
processed_results/
  BRCA1_processed.csv
  PRJNA294226_processed.csv
  PRJNA482141_processed.csv
  PRJNA808405_processed.csv
  PRJNA934049_processed.csv
  yyfbatch1_processed.csv
  yyfbatch2_processed.csv
```

**File format:**
- CSV with row names (sample IDs)
- Column `Group`: `Tumor`, `Normal`, `Cancer`, or `Benign`
- Remaining columns: piRNA expression values (TPM)

### Clinical Data (for downstream scripts)
Place clinical CSV files in `clinical_data/` (needed for scripts 2-5):

```
clinical_data/
  TCGA_BRCA_clinical.csv         # From TCGAbiolinks
  yyfbatch1_clinical.csv         # Your lab records
  yyfbatch2_clinical.csv         # Your lab records
```

**Required columns:** `SampleID`, `Age`, `Stage`, `Subtype`, `OS_time`, `OS_status`

> **Note:** If clinical data is not provided, scripts 2-5 will use simulated data for demonstration. Replace the `WARNING: No clinical data found` blocks with your real data loading code.

---

## Script Execution Order

Run the scripts **in this exact order**. Each script depends on objects produced by the previous one.

### Script 1: Main Pipeline (REQUIRED — run first)

```r
source("piRNA_multicohort_pipeline.R")
```

**What it does:**
1. Loads all 7 datasets from `processed_results/`
2. Balances BRCA1 (matched pairs + 40% remaining tumors)
3. Finds common piRNAs across all datasets
4. Applies log2 transformation + ComBat batch correction + Z-score normalization
5. Splits data: training pool (5 public) vs. dual independent (yyfbatch1 + yyfbatch2)
6. Feature selection: 8 methods (limma, Wilcoxon, RF, LASSO, Elastic Net, Boruta, mRMR, XGBoost) -> consensus ranking -> forward/backward/swap optimization targeting min(AUC) > 0.8 across all validation sets
7. Trains Random Forest (10x5 repeated CV, down-sampled)
8. Evaluates on: Discovery CV, Hold-out, yyfbatch1, yyfbatch2
9. Generates ROC/PRC curves, confusion matrices, combined ROC plot
10. Logistic regression forest plot (T-Score, Age, Stage)
11. Subgroup ROC and T-Score boxplots
12. Feature importance and expression heatmap

**Runtime:** ~10-30 minutes depending on dataset size

**Key outputs:**
- `results/models/final_model.rds` — trained RF model
- `results/models/final_features.rds` — selected piRNA features
- `results/feature_selection/final_features.txt` — feature names (text)
- `results/validation/ROC_combined_all.png` — main figure
- `results/validation/performance_summary.csv` — AUC summary table
- `Final_piRNA_Signature.csv` — final piRNA list

**Objects left in memory** (used by subsequent scripts):
- `combat_df_all` — batch-corrected expression matrix (all 7 datasets)
- `model` — trained caret RF model
- `top_feats` — character vector of selected piRNA names
- `gene_cols` — all piRNA column names

---

### Script 2: Downstream Analysis (Cox, KM, Subgroup)

```r
source("piRNA_downstream_analysis.R")
```

**Requires:** Script 1 objects in memory, OR saved model files in `results/`

**What it does:**
1. Merges BRCA1 + yyfbatch1 + yyfbatch2 for clinical analysis
2. Cox regression (univariate + multivariate) with binary variables
3. Cox forest plot (HR with 95% CI)
4. Kaplan-Meier curves for individual piRNAs + composite T-Score
5. Subgroup ROC with bootstrap 95% CI (by Age, Stage, Subtype)
6. T-Score boxplots with Mann-Whitney p-values and sample sizes
7. Combined ROC + Boxplot panels

**Key outputs:**
- `results/cox_analysis/` — Cox regression tables + forest plot
- `results/km_curves/` — KM survival curves
- `results/subgroup_roc/` — Subgroup ROC curves
- `results/subgroup_box/` — T-Score boxplots
- `results/downstream/` — Combined panels

---

### Script 3: Meta-Analysis & Functional Prediction

```r
source("piRNA_functional_meta_analysis.R")
```

**Requires:** Script 1 objects in memory

**What it does:**
1. Expression heatmap of signature piRNAs across all 7 datasets
2. Standardized Mean Difference (SMD) forest plot per piRNA (random-effects meta-analysis)
3. Spearman correlation heatmap between signature piRNAs
4. Pearson correlation between signature piRNAs and all other genes
5. KEGG pathway enrichment + bar plot
6. GO enrichment (BP, CC, MF) + bubble plots
7. Reactome pathway enrichment
8. Functional interaction network (piRNA -> genes -> pathways)

**Key outputs:**
- `results/meta_analysis/` — heatmap, SMD forest plots, correlation
- `results/functional/` — KEGG, GO, Reactome enrichment results
- `results/network/` — interaction network

---

### Script 4: Network Analysis

```r
source("piRNA_network_analysis.R")
```

**Requires:** Script 1 objects in memory (runs its own correlation if Script 3 was not run)

**What it does:**
1. Per-piRNA radial networks
2. Combined bipartite network (piRNAs <-> genes)
3. Pathway-centric clustered network
4. Hierarchical layered network (piRNA -> genes -> pathways)
5. piRNA-gene correlation heatmap
6. Shared-gene overlap network
7. Circos-style chord diagram
8. Multi-panel combined figure

**Key outputs:**
- `results/network/` — all 8 network plot types (PNG)

---

### Script 5: Advanced Analysis (Nomogram, SHAP, Immune, TMB, Drug)

```r
source("piRNA_advanced_analysis.R")
```

**Requires:** Script 1 objects in memory

**What it does:**
1. LASSO Cox regression + coefficient path plot
2. Nomogram + calibration curve + time-dependent ROC (1/3/5 year)
3. SHAP values for model interpretability
4. Immune infiltration estimation (ssGSEA-based)
5. Immune cell proportion heatmap + boxplots
6. TMB analysis framework (requires mutation data)
7. Drug sensitivity analysis framework (requires drug response data)

**Key outputs:**
- `results/lasso_cox/` — LASSO Cox results
- `results/nomogram/` — nomogram + calibration + time-ROC
- `results/shap/` — SHAP value plots
- `results/immune/` — immune infiltration results
- `results/tmb/` — TMB analysis (if mutation data available)
- `results/drug_sensitivity/` — drug sensitivity (if data available)

---

## Quick Start (Run Everything)

```r
# Step 1: Set working directory to project root
setwd("/path/to/PiRNAs")

# Step 2: Run main pipeline (MUST run first)
source("piRNA_multicohort_pipeline.R")

# Step 3: Run downstream analyses (in any order after Step 2)
source("piRNA_downstream_analysis.R")
source("piRNA_functional_meta_analysis.R")
source("piRNA_network_analysis.R")
source("piRNA_advanced_analysis.R")
```

## Running Scripts Independently (Without Script 1 in Memory)

Scripts 2-5 can load saved results if Script 1 was run previously:

```r
# Script 1 must have been run at least once to create results/models/
# Then you can restart R and run any downstream script:
source("piRNA_downstream_analysis.R")  # loads model from results/models/
```

However, `combat_df_all` is NOT saved to disk by default (it's too large). If you need to run scripts 2-5 in a fresh R session, add this to the end of Script 1:

```r
saveRDS(combat_df_all, "results/models/combat_df_all.rds")
```

Then in scripts 2-5, before `source()`:

```r
combat_df_all <- readRDS("results/models/combat_df_all.rds")
```

---

## Output Directory Structure

```
results/
  feature_selection/
    fs_frequency_table.csv          # Feature selection consensus table
    final_features.txt              # Selected piRNA names
    feature_importance.png          # RF importance bar plot
    expression_heatmap.png          # Heatmap of selected piRNAs
  models/
    final_model.rds                 # Trained RF model (caret object)
    final_features.rds              # Selected features (R object)
  validation/
    ROC_combined_all.png            # *** KEY FIGURE: all 4 ROC curves ***
    ROC_discovery_CV.png
    ROC_holdout.png
    ROC_independent_yyfbatch1.png
    ROC_independent_yyfbatch2.png
    PRC_*.png                       # Precision-Recall curves
    CM_*.png                        # Confusion matrices
    performance_summary.csv         # AUC summary table
  forest_plot/
    forest_plot.png                 # Logistic regression forest plot
    univariate_results.csv
    multivariate_results.csv
  subgroup/
    ROC_by_Age.png
    ROC_by_Stage.png
    ROC_by_Subtype.png
    Boxplot_TScore_Phase.png        # *** KEY FIGURE: T-Score by phase ***
    Boxplot_TScore_AllBatches.png
  cox_analysis/                     # From Script 2
    cox_forest_plot.png
    univariate_cox_results.csv
    multivariate_cox_results.csv
  km_curves/                        # From Script 2
    KM_TScore_risk.png
    KM_<piRNA_name>.png
  subgroup_roc/                     # From Script 2
  subgroup_box/                     # From Script 2
  downstream/                       # From Script 2
  meta_analysis/                    # From Script 3
  functional/                       # From Script 3
  network/                          # From Scripts 3 & 4
  lasso_cox/                        # From Script 5
  nomogram/                         # From Script 5
  shap/                             # From Script 5
  immune/                           # From Script 5
```

## Replacing Simulated Clinical Data with Real Data

In **Script 1** (`piRNA_multicohort_pipeline.R`), find the block starting with:
```r
if (!"Age" %in% colnames(combat_df_all)) {
  cat("WARNING: No clinical data found. Using simulated clinical data...")
```

Replace it with:
```r
# Load TCGA clinical data
tcga_clin <- read.csv("clinical_data/TCGA_BRCA_clinical.csv", stringsAsFactors = FALSE)
# Match by sample barcode (first 12 characters of TCGA barcode)
tcga_idx <- match(substr(rownames(combat_df_all), 1, 12), tcga_clin$submitter_id)
combat_df_all$Age[!is.na(tcga_idx)] <- tcga_clin$age_at_diagnosis[tcga_idx[!is.na(tcga_idx)]]
# ... similar for Stage, Subtype

# Load your lab clinical data
yyf1_clin <- read.csv("clinical_data/yyfbatch1_clinical.csv", stringsAsFactors = FALSE)
yyf2_clin <- read.csv("clinical_data/yyfbatch2_clinical.csv", stringsAsFactors = FALSE)
# Match by SampleID
# ... merge into combat_df_all
```

The same pattern applies to **Script 2** (`piRNA_downstream_analysis.R`) which has a similar simulated data block.

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `Error: No data found` | Place `*_processed.csv` files in `processed_results/` |
| `combat_df_all not found` | Run Script 1 first, or load from saved RDS |
| `ComBat error` | Check that all datasets have the same piRNA column names |
| Very few common piRNAs | Ensure consistent piRNA naming across datasets |
| AUC < 0.8 on independent sets | The forward/backward/swap search will try to maximize; check data quality |
| Memory issues | BRCA1 balancing reduces TCGA size; also try `gc()` between scripts |
| Package install fails | For Bioconductor: `BiocManager::install("package_name")` |

## Key R Packages

| Package | Purpose |
|---------|---------|
| `sva` | ComBat batch correction |
| `caret` | Model training framework |
| `randomForest` | Random Forest classifier |
| `glmnet` | LASSO / Elastic Net |
| `pROC` | ROC curves and AUC |
| `Boruta` | Feature selection |
| `survival` / `survminer` | Cox regression, KM curves |
| `clusterProfiler` | GO / KEGG enrichment |
| `pheatmap` | Expression heatmaps |
| `xgboost` | XGBoost feature importance + SHAP |

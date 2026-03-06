# Breast Cancer piRNA Diagnostic Signature — Analysis Protocol

## Overview

This project builds a piRNA-based diagnostic and prognostic signature for breast cancer using multi-cohort data. The pipeline selects <=10 piRNA features, trains a Random Forest classifier on 5 public datasets, validates independently on 2 in-house cohorts (yyfbatch1 + yyfbatch2), and performs downstream survival, immune, and drug-sensitivity analyses.

---

## Required Data Files (Checklist)

**Before running anything, make sure you have ALL of the files below.**
If a file is missing, the corresponding script will fail (no simulated fallbacks).

### 1. piRNA Expression Matrices

Place in `processed_results/`:

| File | Dataset | Role |
|------|---------|------|
| `BRCA1_processed.csv` | TCGA-BRCA | Training (balanced) |
| `PRJNA294226_processed.csv` | GEO: GSE72080 | Training |
| `PRJNA482141_processed.csv` | GEO: GSE117452 | Training |
| `PRJNA808405_processed.csv` | GEO: GSE197020 | Training |
| `PRJNA934049_processed.csv` | GEO: GSE225117 | Training |
| `yyfbatch1_processed.csv` | In-house | Independent validation |
| `yyfbatch2_processed.csv` | In-house | Independent validation |

**Format:**
- CSV with row names = sample IDs
- Column `Group`: `Tumor` or `Normal` (or `Cancer`/`Benign`)
- Remaining columns: piRNA expression values (TPM)
- All files must share a common set of piRNA column names

### 2. Clinical / Phenotype Data

Place in `clinical_data/`:

| File | Used by | Required columns |
|------|---------|-----------------|
| `TCGA_BRCA_clinical.csv` | Scripts 1-6 | `SampleID`, `Age`, `Stage`, `Subtype`, `OS_time`, `OS_status` |
| `yyfbatch1_clinical.csv` | Scripts 1-6 | `SampleID`, `Age`, `Stage`, `Subtype`, `OS_time`, `OS_status` |
| `yyfbatch2_clinical.csv` | Scripts 1-6 | `SampleID`, `Age`, `Stage`, `Subtype`, `OS_time`, `OS_status` |

**Column details:**

| Column | Type | Description |
|--------|------|-------------|
| `SampleID` | character | Must match row names in the expression CSVs |
| `Age` | numeric | Age at diagnosis (years) |
| `Stage` | character | `Stage I`, `Stage II`, `Stage III`, or `Stage IV` |
| `Subtype` | character | e.g. `Luminal A`, `Luminal B`, `HER2+`, `Triple-negative` |
| `OS_time` | numeric | Overall survival time (months) |
| `OS_status` | integer | 0 = censored/alive, 1 = dead/event |

> **How to get TCGA clinical data:**
> ```r
> library(TCGAbiolinks)
> clin <- GDCquery_clinic("TCGA-BRCA", type = "clinical")
> # Keep: submitter_id, age_at_diagnosis, ajcc_pathologic_stage,
> #       paper_BRCA_Subtype_PAM50, days_to_death, days_to_last_follow_up, vital_status
> # Rename to SampleID, Age, Stage, Subtype, OS_time, OS_status
> ```

### 3. Matched mRNA Expression (for functional prediction)

| File | Used by | Description |
|------|---------|-------------|
| `mRNA_expression/TCGA_BRCA_mRNA.csv` | Script 3 (functional), Script 5 (immune ssGSEA) | Gene-level expression matrix (rows = samples, columns = gene symbols, values = TPM/FPKM) |

**Format:** Same sample IDs as piRNA expression. Columns are HGNC gene symbols.

Without this file:
- Script 3: piRNA-gene correlation, KEGG/GO/Reactome enrichment will be skipped
- Script 5: immune infiltration will not run ssGSEA (no immune marker genes available)

### 4. TMB Data (for TMB analysis in Script 5)

| File | Used by | Required columns |
|------|---------|-----------------|
| `clinical_data/tmb_data.csv` | Script 5 (Phase 5) | `SampleID`, `TMB` |

- `TMB`: mutations per megabase (numeric)
- Can be derived from MAF files using `maftools::tmb()`

### 5. Drug Sensitivity Data (for drug analysis in Script 5)

| File | Used by | Required columns |
|------|---------|-----------------|
| `clinical_data/drug_ic50.csv` | Script 5 (Phase 6) | `SampleID`, one column per drug (IC50 values) |

- Generate with `oncoPredict` or `pRRophetic` from matched mRNA expression
- Example drug columns: `Tamoxifen`, `Paclitaxel`, `Doxorubicin`, `Cisplatin`, etc.

### Summary: What do I need?

| Priority | File(s) | Scripts that need it |
|----------|---------|---------------------|
| **REQUIRED** | 7 x `*_processed.csv` in `processed_results/` | All scripts |
| **REQUIRED** | 3 x clinical CSVs in `clinical_data/` | Scripts 1-6 (Cox, KM, nomogram, subgroup, prognostic figure) |
| Recommended | `mRNA_expression/TCGA_BRCA_mRNA.csv` | Script 3 (KEGG/GO/Reactome), Script 5 (immune ssGSEA) |
| Optional | `clinical_data/tmb_data.csv` | Script 5 (TMB analysis) |
| Optional | `clinical_data/drug_ic50.csv` | Script 5 (drug sensitivity) |

---

## Prerequisites

### R Version
- R >= 4.1.0

### Key R Packages

| Package | Purpose |
|---------|---------|
| `sva` | ComBat batch correction |
| `caret` | Model training framework |
| `randomForest` | Random Forest classifier |
| `glmnet` | LASSO / Elastic Net / Ridge |
| `pROC` | ROC curves and AUC |
| `Boruta` | Feature selection |
| `survival` / `survminer` | Cox regression, KM curves |
| `timeROC` | Time-dependent ROC |
| `rms` | Nomogram + calibration |
| `gbm` | Gradient Boosting Machines |
| `xgboost` | XGBoost + SHAP |
| `SHAPforxgboost` | SHAP value plots |
| `GSVA` / `GSEABase` | ssGSEA immune estimation |
| `clusterProfiler` | GO / KEGG enrichment |
| `ReactomePA` | Reactome enrichment |
| `metafor` / `meta` | Meta-analysis |
| `pheatmap` | Expression heatmaps |
| `igraph` / `ggraph` | Network visualization |

---

## Script Execution Order

Run the scripts **in this exact order**. Each script depends on objects produced by the previous one.

### Script 1: Main Pipeline (REQUIRED — run first)

```r
source("piRNA_multicohort_pipeline.R")
```

**What it does:**
1. Loads all 7 datasets from `processed_results/`
2. Loads clinical data from `clinical_data/` and merges by SampleID
3. Balances BRCA1 (matched pairs + 40% remaining tumors)
4. Finds common piRNAs across all datasets
5. Applies log2 transformation + ComBat batch correction + Z-score normalization
6. Splits data: training pool (5 public) vs. dual independent (yyfbatch1 + yyfbatch2)
7. Feature selection: 8 methods -> consensus ranking -> forward/backward/swap optimization
8. Trains Random Forest (10x5 repeated CV, down-sampled)
9. Evaluates on: Discovery CV, Hold-out, yyfbatch1, yyfbatch2
10. Generates ROC/PRC curves, confusion matrices, combined ROC plot
11. Logistic regression forest plot (T-Score, Age, Stage)
12. Subgroup ROC and T-Score boxplots
13. Feature importance and expression heatmap

**Data required:** expression CSVs + clinical CSVs

**Key outputs:**
- `results/models/final_model.rds` — trained RF model
- `results/models/final_features.rds` — selected piRNA features
- `results/feature_selection/final_features.txt` — feature names (text)
- `results/validation/ROC_combined_all.png` — main figure
- `results/validation/performance_summary.csv` — AUC summary table

**Objects left in memory** (used by subsequent scripts):
- `combat_df_all` — batch-corrected expression matrix with clinical columns
- `model` — trained caret RF model
- `top_feats` — character vector of selected piRNA names
- `gene_cols` — all piRNA column names

---

### Script 2: Downstream Analysis (Cox, KM, Subgroup)

```r
source("piRNA_downstream_analysis.R")
```

**Requires:** Script 1 objects + clinical data already merged in `combat_df_all`

**What it does:**
1. Merges BRCA1 + yyfbatch1 + yyfbatch2 for clinical analysis
2. Cox regression (univariate + multivariate) with binary variables
3. Cox forest plot (HR with 95% CI)
4. Kaplan-Meier survival curves for individual piRNAs + composite T-Score
5. Subgroup ROC curves (Age, Stage, Subtype) with 95% CI
6. T-Score boxplots with Mann-Whitney p-values

**Data required:** `OS_time`, `OS_status`, `Age`, `Stage`, `Subtype` columns in `combat_df_all`

**Key outputs:**
- `results/cox_analysis/` — Cox regression tables + forest plot
- `results/km_curves/` — KM survival curves
- `results/subgroup_roc/` — Subgroup ROC curves
- `results/downstream/` — Combined panels

---

### Script 3: Meta-Analysis & Functional Prediction

```r
source("piRNA_functional_meta_analysis.R")
```

**Requires:** Script 1 objects + mRNA expression data (for Part B)

**What it does:**
1. Expression heatmap of signature piRNAs across all 7 datasets
2. Standardized Mean Difference (SMD) forest plot per piRNA
3. Spearman correlation heatmap between signature piRNAs
4. Pearson correlation between signature piRNAs and mRNA genes
5. KEGG pathway enrichment + bar plot
6. GO enrichment (BP, CC, MF) + bubble plots
7. Reactome pathway enrichment
8. Functional interaction network (piRNA -> genes -> pathways)

**Data required:** `mRNA_expression/TCGA_BRCA_mRNA.csv` for items 4-8

**Key outputs:**
- `results/meta_analysis/` — heatmap, SMD forest plots, correlation
- `results/functional/` — KEGG, GO, Reactome enrichment results
- `results/network/` — interaction network

---

### Script 4: Network Analysis

```r
source("piRNA_network_analysis.R")
```

**Requires:** Script 1 objects (runs its own correlation if Script 3 was not run)

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

**Requires:** Script 1 objects + clinical data + optional TMB & drug data

**What it does:**
1. LASSO Cox regression + coefficient path plot
2. Nomogram + calibration curve + time-dependent ROC (1/3/5 year)
3. SHAP values for model interpretability
4. Immune infiltration estimation (ssGSEA — requires mRNA data)
5. Immune cell proportion heatmap + boxplots
6. TMB analysis (requires `clinical_data/tmb_data.csv`)
7. Drug sensitivity analysis (requires `clinical_data/drug_ic50.csv`)

**Data required:**
- Clinical columns (`OS_time`, `OS_status`, `Age`, `Stage`) — **required**
- `mRNA_expression/TCGA_BRCA_mRNA.csv` — for immune ssGSEA
- `clinical_data/tmb_data.csv` — for TMB analysis
- `clinical_data/drug_ic50.csv` — for drug sensitivity

**Key outputs:**
- `results/lasso_cox/` — LASSO Cox coefficient path + CV curve
- `results/nomogram/` — nomogram + calibration + time-ROC
- `results/shap/` — SHAP value plots
- `results/immune/` — immune infiltration results
- `results/tmb/` — TMB analysis
- `results/drug_sensitivity/` — drug sensitivity

---

### Script 6: Prognostic Model Figure (NEW)

```r
source("pipeline/12_prognostic_model_figure.R")
```

**Requires:** Script 1 objects + clinical data (OS_time, OS_status)

**What it does:**
1. **Panel A** — Univariate Cox forest plot of core piRNA features (HR + 95% CI)
2. **Panel B** — ML pipeline comparison heatmap: 8 feature-selection methods x 7 classifiers = 56 pipelines, ranked by C-index across training and validation cohorts
3. **Panel C** — Kaplan-Meier survival curve for training cohort (best pipeline)
4. **Panel D** — Kaplan-Meier survival curve for validation cohort (best pipeline)
5. **Panel E** — Time-dependent ROC at 1/3/5 yr, training cohort
6. **Panel F** — Time-dependent ROC at 1/3/5 yr, validation cohort
7. Automatically selects the best-performing pipeline for panels C-F
8. Assembles all 6 panels into a composite figure

**Feature-selection methods:** Lasso, Ridge, ElasticNet, StepCox (forward/backward/both), RSF

**Classifiers:** CoxPH, Lasso, Ridge, ElasticNet, GBM, RSF, SVM

**Data required:**
- `OS_time` and `OS_status` columns in `combat_df_all` — **required**
- Clinical data must be real (survival endpoints); this script has no simulated fallback

**Key outputs:**
- `results/prognostic_figure/Figure_prognostic_model.png` — composite 6-panel figure
- `results/prognostic_figure/Figure_prognostic_model.pdf` — vector PDF
- `results/prognostic_figure/pipeline_comparison.csv` — all 56 pipeline C-index results
- Individual panel PNGs in `results/prognostic_figure/`

---

## Quick Start (Run Everything)

```r
# Step 1: Set working directory to project root
setwd("/path/to/BRCA_piRNA")

# Step 2: Run main pipeline (MUST run first)
source("piRNA_multicohort_pipeline.R")

# Step 3: Run downstream analyses (in any order after Step 2)
source("piRNA_downstream_analysis.R")
source("piRNA_functional_meta_analysis.R")
source("piRNA_network_analysis.R")
source("piRNA_advanced_analysis.R")

# Step 4: Run prognostic model figure
source("pipeline/12_prognostic_model_figure.R")
```

## Running Scripts Independently (Without Script 1 in Memory)

Scripts 2-6 can load saved results if Script 1 was run previously.
However, `combat_df_all` is NOT saved to disk by default. Add this to the end of Script 1:

```r
saveRDS(combat_df_all, "results/models/combat_df_all.rds")
```

Then in scripts 2-6, before `source()`:

```r
combat_df_all <- readRDS("results/models/combat_df_all.rds")
```

---

## Directory Structure

```
BRCA_piRNA/
│
├── processed_results/                    # INPUT: piRNA expression matrices
│   ├── BRCA1_processed.csv
│   ├── PRJNA294226_processed.csv
│   ├── PRJNA482141_processed.csv
│   ├── PRJNA808405_processed.csv
│   ├── PRJNA934049_processed.csv
│   ├── yyfbatch1_processed.csv
│   └── yyfbatch2_processed.csv
│
├── clinical_data/                        # INPUT: clinical/phenotype data
│   ├── TCGA_BRCA_clinical.csv           # SampleID, Age, Stage, Subtype, OS_time, OS_status
│   ├── yyfbatch1_clinical.csv
│   ├── yyfbatch2_clinical.csv
│   ├── tmb_data.csv                     # (optional) SampleID, TMB
│   └── drug_ic50.csv                    # (optional) SampleID, drug IC50 columns
│
├── mRNA_expression/                      # INPUT: matched mRNA data
│   └── TCGA_BRCA_mRNA.csv              # (recommended) gene expression matrix
│
├── piRNA_multicohort_pipeline.R          # Script 1: main pipeline
├── piRNA_downstream_analysis.R           # Script 2: Cox, KM, subgroup
├── piRNA_functional_meta_analysis.R      # Script 3: meta-analysis, enrichment
├── piRNA_network_analysis.R              # Script 4: network visualization
├── piRNA_advanced_analysis.R             # Script 5: nomogram, SHAP, immune, TMB, drug
├── pipeline/
│   ├── 00_config.R ... 11_select_best_independent_cohort.R
│   └── 12_prognostic_model_figure.R      # Script 6: prognostic model figure
│
└── results/                              # OUTPUT: all generated results
    ├── feature_selection/
    ├── models/
    ├── validation/
    ├── forest_plot/
    ├── subgroup/
    ├── cox_analysis/
    ├── km_curves/
    ├── subgroup_roc/
    ├── downstream/
    ├── meta_analysis/
    ├── functional/
    ├── network/
    ├── lasso_cox/
    ├── nomogram/
    ├── shap/
    ├── immune/
    ├── tmb/
    ├── drug_sensitivity/
    └── prognostic_figure/                # NEW: Script 6 outputs
```

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `Error: No data found` | Place `*_processed.csv` files in `processed_results/` |
| `combat_df_all not found` | Run Script 1 first, or load from saved RDS |
| Missing `OS_time`/`OS_status` | Provide clinical CSVs with survival data |
| Missing `Age`/`Stage`/`Subtype` | Provide clinical CSVs with demographic columns |
| `ComBat error` | Check that all datasets have the same piRNA column names |
| Very few common piRNAs | Ensure consistent piRNA naming across datasets |
| AUC < 0.8 on independent sets | Check data quality; the optimization will try to maximize |
| Immune ssGSEA fails | Provide matched mRNA expression data |
| TMB/Drug sections skip | Provide `tmb_data.csv` / `drug_ic50.csv` |
| Memory issues | BRCA1 balancing reduces TCGA size; also try `gc()` between scripts |
| Package install fails | For Bioconductor: `BiocManager::install("package_name")` |

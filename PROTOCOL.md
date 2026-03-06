# Breast Cancer piRNA Diagnostic Signature — Analysis Protocol

## Overview

This project builds a piRNA-based diagnostic and prognostic signature for breast cancer using multi-cohort data. The pipeline selects <=10 piRNA features, trains a Random Forest classifier on 5 public datasets, validates independently on 2 in-house cohorts (yyfbatch1 + yyfbatch2), and performs downstream survival, immune, and drug-sensitivity analyses.

---

## Pipeline Workflow

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                         DATA PREPARATION (Step 0)                               │
│                                                                                 │
│  ┌──────────────────────┐   ┌──────────────────────┐   ┌────────────────────┐   │
│  │ 5 Public Datasets    │   │ yyfbatch1 (n≈192)    │   │ yyfbatch2 (n≈192)  │   │
│  │ BRCA1, PRJNA294226,  │   │ In-house validation  │   │ In-house validation│   │
│  │ PRJNA482141,         │   │ Survival: 31% avail  │   │ Survival: 0%       │   │
│  │ PRJNA808405,         │   └──────────┬───────────┘   └─────────┬──────────┘   │
│  │ PRJNA934049          │              │                         │               │
│  └──────────┬───────────┘              │                         │               │
│             │                          │                         │               │
│  ┌──────────┴──────────────────────────┴─────────────────────────┴──────────┐    │
│  │  Reformat clinical data: SampleID, Age, Stage, Subtype, OS_time/status  │    │
│  │  TCGA via TCGAbiolinks  ·  yyfbatch: rename AJCC_Stage → Stage          │    │
│  └─────────────────────────────────┬───────────────────────────────────────-┘    │
└────────────────────────────────────┼────────────────────────────────────────────-┘
                                     ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│              STEP 1: DIAGNOSTIC MODEL  (piRNA_multicohort_pipeline.R)            │
│                                                                                 │
│  7 datasets ──► Common piRNAs ──► log2 ──► ComBat ──► Z-score                  │
│                                                                                 │
│  ┌─────────────────────────────┐    ┌──────────────────────────────────────┐     │
│  │ TRAINING (5 public)         │    │ VALIDATION (independent)             │     │
│  │                             │    │                                      │     │
│  │ 8 Feature Selection Methods │    │  yyfbatch1: Dx validation            │     │
│  │ → Consensus (≤10 piRNAs)    │    │  yyfbatch2: Dx validation            │     │
│  │ → Forward/Backward/Swap     │    │                                      │     │
│  │ → Random Forest (10×5 CV)   │    │  ROC, PRC, Confusion Matrices       │     │
│  └──────────────┬──────────────┘    └──────────────────┬───────────────────┘     │
│                 │                                      │                         │
│                 └──────────┬───────────────────────────-┘                         │
│                            ▼                                                     │
│              final_model.rds + final_features.rds                               │
│              combat_df_all (in memory for downstream)                            │
└────────────────────────────┼────────────────────────────────────────────────────-┘
                             │
              ┌──────────────┼──────────────┐
              ▼              ▼              ▼
┌─────────────────┐ ┌───────────────┐ ┌────────────────┐
│ STEP 2          │ │ STEP 3        │ │ STEP 4         │
│ Downstream      │ │ Meta-Analysis │ │ Network        │
│ Analysis        │ │ & Functional  │ │ Analysis       │
│                 │ │               │ │                │
│ Cox regression  │ │ SMD forest    │ │ 8 network      │
│  (TCGA only)    │ │ Expression    │ │ visualizations │
│ KM curves       │ │   heatmap     │ │ (radial,       │
│  (TCGA only)    │ │ piRNA-mRNA    │ │  bipartite,    │
│ Subgroup ROC    │ │   correlation │ │  chord, etc.)  │
│ T-Score boxplot │ │ KEGG/GO/      │ │                │
│                 │ │  Reactome     │ │                │
│                 │ │ (needs mRNA)  │ │                │
└────────┬────────┘ └───────┬───────┘ └───────┬────────┘
         │                  │                 │
         └──────────────────┼─────────────────┘
                            ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│              STEP 5: ADVANCED ANALYSIS  (piRNA_advanced_analysis.R)              │
│                                                                                 │
│  ┌─────────────────┐ ┌──────────────┐ ┌──────────┐ ┌──────────┐ ┌───────────┐  │
│  │ LASSO-Cox       │ │ Nomogram +   │ │ SHAP     │ │ Immune   │ │ Drug      │  │
│  │ (TCGA only)     │ │ Calibration  │ │ (XGBoost)│ │ ssGSEA   │ │ Sensitiv. │  │
│  │ Coef path plot  │ │ + tdROC      │ │          │ │ (needs   │ │ (needs    │  │
│  │                 │ │ (TCGA only)  │ │          │ │  mRNA)   │ │  mRNA +   │  │
│  │                 │ │              │ │          │ │          │ │ oncoPredict│  │
│  └─────────────────┘ └──────────────┘ └──────────┘ └──────────┘ └───────────┘  │
│                                                                                 │
│  ┌────────────────────────────────┐                                              │
│  │ TMB Analysis (optional)       │                                              │
│  │ Needs: clinical_data/tmb.csv  │                                              │
│  └────────────────────────────────┘                                              │
└────────────────────────────────────────┬────────────────────────────────────────-┘
                                         ▼
┌─────────────────────────────────────────────────────────────────────────────────┐
│         STEP 6: PROGNOSTIC MODEL FIGURE  (pipeline/12_prognostic_model_figure.R) │
│                                                                                 │
│  6-panel composite (TCGA survival data only):                                   │
│  ┌───────┐ ┌───────────────┐ ┌───────────┐ ┌───────────┐ ┌──────┐ ┌──────┐     │
│  │ A:    │ │ B: ML Pipeline│ │ C: KM     │ │ D: KM     │ │ E:   │ │ F:   │     │
│  │ Uni-  │ │ Heatmap       │ │ (train)   │ │ (valid)   │ │ tdROC│ │ tdROC│     │
│  │ variate│ │ 7 FS × 7 ML  │ │           │ │           │ │ train│ │ valid│     │
│  │ Cox   │ │ C-index       │ │           │ │           │ │      │ │      │     │
│  └───────┘ └───────────────┘ └───────────┘ └───────────┘ └──────┘ └──────┘     │
│                                                                                 │
│  Output: Figure_prognostic_model.{png,pdf}                                      │
└─────────────────────────────────────────────────────────────────────────────────┘

  ╔═══════════════════════════════════════════════════════════════════╗
  ║  SURVIVAL DATA SCOPE:                                           ║
  ║  All survival analyses (Cox, KM, Nomogram, LASSO-Cox, tdROC,    ║
  ║  prognostic figure) use TCGA-BRCA data ONLY.                    ║
  ║  yyfbatch1: 69% missing  ·  yyfbatch2: 100% missing            ║
  ║                                                                 ║
  ║  DIAGNOSTIC MODEL: Uses Group labels only (all 7 cohorts).      ║
  ╚═══════════════════════════════════════════════════════════════════╝
```

---

## Data Readiness Checklist

### Status of All Required Input Files

**Legend:** READY = file exists and is usable; NEEDS PREP = file exists but needs reformatting; MISSING = not yet available; COMPUTED = generated by the pipeline itself.

#### 1. piRNA Expression Matrices (`processed_results/`)

| File | Status | Notes |
|------|--------|-------|
| `BRCA1_processed.csv` | READY (upload full file) | Training. User has locally, only a subset (`BRCA1_processed_1.csv`) is in repo |
| `PRJNA294226_processed.csv` | READY (upload) | Training. User has locally |
| `PRJNA482141_processed.csv` | READY (upload) | Training. User has locally |
| `PRJNA808405_processed.csv` | READY (upload) | Training. User has locally |
| `PRJNA934049_processed.csv` | READY (upload) | Training. User has locally |
| `yyfbatch1_processed.csv` | READY | In-house validation |
| `yyfbatch2_processed.csv` | READY | In-house validation |

**Format:** CSV with row names = sample IDs, column `Group` = `Tumor`/`Normal`/`Cancer`/`Benign`, remaining columns = piRNA TPM values. All files must share common piRNA column names.

#### 2. Clinical / Phenotype Data (`clinical_data/`)

The main pipeline (`piRNA_multicohort_pipeline.R`, line 865) expects clinical files with these **exact column names**: `SampleID`, `Age`, `Stage`, `Subtype`, `OS_time`, `OS_status`.

| File | Status | Action Needed |
|------|--------|---------------|
| `TCGA_BRCA_clinical.csv` | NEEDS PREP | Current file has raw TCGA column names (`submitter_id.samples`, `tumor_stage.diagnoses`, etc.) and is missing `Age`, `Subtype`, `OS_time`, `OS_status`. **Must be reformatted** — see Step 0 below |
| `BRCA_pheno.csv` | READY | Contains sample-to-group mapping. Used during clinical merge |
| `yyfbatch1_clinical_clean.csv` | NEEDS PREP | Has `AJCC_Stage` instead of `Stage`; `OS_time`/`OS_status` are 69% empty; `Subtype` 45% empty |
| `yyfbatch2_clinical_clean.csv` | NEEDS PREP | Has `AJCC_Stage` instead of `Stage`; `OS_time`/`OS_status` are **100% empty**; `Subtype` 44% empty |

**Critical column gap summary:**

| Column | TCGA | yyfbatch1 | yyfbatch2 |
|--------|------|-----------|-----------|
| SampleID | Missing (has `submitter_id.samples`) | OK | OK |
| Age | Missing (need TCGAbiolinks) | OK (0% missing) | OK (1% missing) |
| Stage | Missing (has `tumor_stage.diagnoses`) | Partial (has `AJCC_Stage`, 51% missing) | Partial (has `AJCC_Stage`, 44% missing) |
| Subtype | Missing entirely | 45% missing | 44% missing |
| OS_time | Missing (has `days_to_death`/`days_to_last_follow_up`) | 69% missing | **100% missing** |
| OS_status | Missing (has `vital_status.diagnoses`) | 69% missing | **100% missing** |

#### 3. mRNA Expression (`mRNA_expression/`)

| File | Status | Used by |
|------|--------|---------|
| `TCGA_BRCA_mRNA.csv` | MISSING | Script 3 (piRNA-gene correlation, KEGG/GO/Reactome enrichment), Script 5 (immune ssGSEA) |

**How to get it:**
```r
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-BRCA", data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", workflow.type = "STAR - Counts")
GDCdownload(query)
data <- GDCprepare(query)
# Extract TPM matrix, transpose to samples x genes, save as CSV
```

**Without this file:** piRNA-gene correlation, pathway enrichment (KEGG/GO/Reactome), network analysis, and immune ssGSEA will all be skipped. The diagnostic pipeline (Script 1) still works.

#### 4. TMB Data (`clinical_data/tmb_data.csv`)

| File | Status | Used by |
|------|--------|---------|
| `tmb_data.csv` | MISSING (optional) | Script 5, Phase 5 |

**How to get it:**
```r
library(TCGAbiolinks)
library(maftools)
query <- GDCquery(project = "TCGA-BRCA", data.category = "Simple Nucleotide Variation",
                  data.type = "Masked Somatic Mutation", access = "open")
GDCdownload(query)
maf <- GDCprepare(query)
maf_obj <- read.maf(maf)
tmb_df <- tmb(maf_obj)  # mutations per megabase
# Save with columns: SampleID, TMB
```

#### 5. Drug Sensitivity Data (`clinical_data/drug_ic50.csv`)

| File | Status | Used by |
|------|--------|---------|
| `drug_ic50.csv` | COMPUTED (not a download) | Script 5, Phase 6 |

**This is NOT available from TCGA.** IC50 values are predicted using cell line pharmacogenomics training data:
```r
library(oncoPredict)
# Download GDSC2 training data (bundled with oncoPredict)
# Requires mRNA expression matrix (TCGA_BRCA_mRNA.csv)
calcPhenotype(trainingExprData = GDSC2_Expr, trainingPvalData = GDSC2_Res,
              testExprData = your_mRNA_matrix, batchCorrect = "eb")
# Output: predicted IC50 per drug per sample
```

The current pipeline has a **placeholder stub** for this step. It will run once you provide the mRNA expression data and install `oncoPredict`.

---

## Step 0: Prepare Clinical Data (MUST DO BEFORE RUNNING)

The uploaded `TCGA_BRCA_clinical.csv` needs reformatting to match what the pipeline expects. Run this preparation script:

```r
# === Reformat TCGA_BRCA_clinical.csv ===
library(TCGAbiolinks)

# 1. Download full TCGA-BRCA clinical data (includes Age, Subtype)
clin <- GDCquery_clinic("TCGA-BRCA", type = "clinical")

# 2. Build the standardized clinical file
tcga_clin <- data.frame(
  SampleID    = paste0(clin$submitter_id, "-01A"),  # match expression row names
  Age         = floor(clin$age_at_diagnosis / 365.25),
  Stage       = clin$ajcc_pathologic_stage,
  Subtype     = clin$paper_BRCA_Subtype_PAM50,
  OS_time     = ifelse(!is.na(clin$days_to_death),
                       clin$days_to_death / 30.44,    # convert days -> months
                       clin$days_to_last_follow_up / 30.44),
  OS_status   = ifelse(clin$vital_status == "Dead", 1, 0),
  stringsAsFactors = FALSE
)

write.csv(tcga_clin, "clinical_data/TCGA_BRCA_clinical.csv", row.names = FALSE)

# === Fix yyfbatch column names ===
for (batch in c("yyfbatch1", "yyfbatch2")) {
  path <- paste0("clinical_data/", batch, "_clinical_clean.csv")
  df <- read.csv(path, stringsAsFactors = FALSE)
  if ("AJCC_Stage" %in% colnames(df) && !"Stage" %in% colnames(df)) {
    df$Stage <- df$AJCC_Stage
  }
  write.csv(df, path, row.names = FALSE)
}
```

**Note on survival data:** yyfbatch2 has 0% survival data and yyfbatch1 has only ~31% (69% missing). Survival analyses (Cox regression, KM curves, prognostic model figure) will effectively only use TCGA-BRCA data. The diagnostic model (Script 1) does NOT require survival data — it works with Group labels only.

---

## Script Execution Order

Run in this exact order. Each script depends on the previous one.

### Step 1: Main Pipeline (REQUIRED)

```r
source("piRNA_multicohort_pipeline.R")
```

**Reads:**
- `processed_results/{BRCA1,PRJNA294226,PRJNA482141,PRJNA808405,PRJNA934049,yyfbatch1,yyfbatch2}_processed.csv`
- `clinical_data/TCGA_BRCA_clinical.csv` (expects: `SampleID`, `Age`, `Stage`, `Subtype`, `OS_time`, `OS_status`)
- `clinical_data/yyfbatch1_clinical_clean.csv` (same columns)
- `clinical_data/yyfbatch2_clinical_clean.csv` (same columns)

**Does:**
1. Loads 7 expression datasets, recodes Group labels
2. Balances BRCA1 (matched tumor-normal pairs + 40% excess tumors)
3. Finds common piRNAs across all datasets
4. log2 transform + ComBat batch correction + Z-score normalization
5. Splits: 5 public datasets (training) vs. yyfbatch1 + yyfbatch2 (independent validation)
6. Feature selection: 8 methods -> consensus -> forward/backward/swap optimization (<=10 piRNAs)
7. Random Forest training (10x5 repeated CV, down-sampled)
8. Evaluation: Discovery CV, Hold-out, yyfbatch1, yyfbatch2
9. Merges clinical data by SampleID matching
10. ROC/PRC curves, confusion matrices, combined ROC plot
11. Logistic regression forest plot, subgroup ROC, T-Score boxplots, heatmap

**Outputs:** `results/models/final_model.rds`, `results/models/final_features.rds`, `results/validation/` plots

**Objects in memory** (needed by later scripts): `combat_df_all`, `model`, `top_feats`, `gene_cols`

---

### Step 2: Downstream Analysis (Cox, KM, Subgroup)

```r
source("piRNA_downstream_analysis.R")
```

**Requires:** Step 1 objects. Survival columns (`OS_time`, `OS_status`) must be present.

**Does:** Cox regression (univariate + multivariate), KM survival curves, subgroup ROC, T-Score boxplots

**Outputs:** `results/cox_analysis/`, `results/km_curves/`, `results/subgroup_roc/`

---

### Step 3: Meta-Analysis & Functional Prediction

```r
source("piRNA_functional_meta_analysis.R")
```

**Requires:** Step 1 objects + `mRNA_expression/TCGA_BRCA_mRNA.csv` (for piRNA-gene correlation and enrichment)

**Does:** Expression heatmap, SMD forest plot, correlation heatmap, piRNA-mRNA correlation, KEGG/GO/Reactome enrichment, functional network

**Outputs:** `results/meta_analysis/`, `results/functional/`, `results/network/`

---

### Step 4: Network Analysis

```r
source("piRNA_network_analysis.R")
```

**Requires:** Step 1 objects (runs own correlation if Step 3 was skipped)

**Does:** 8 types of network visualizations (radial, bipartite, pathway-centric, hierarchical, chord diagram, etc.)

**Outputs:** `results/network/`

---

### Step 5: Advanced Analysis (Nomogram, SHAP, Immune, TMB, Drug)

```r
source("piRNA_advanced_analysis.R")
```

**Requires:** Step 1 objects + survival data + optional extras

**Does:**
1. LASSO Cox regression + coefficient path plot
2. Nomogram + calibration + time-dependent ROC (1/3/5 yr)
3. SHAP values (XGBoost model interpretability)
4. Immune infiltration (ssGSEA) — **needs `mRNA_expression/TCGA_BRCA_mRNA.csv`**
5. TMB analysis — **needs `clinical_data/tmb_data.csv`**
6. Drug sensitivity — **needs mRNA data + `oncoPredict` package** (no separate IC50 file needed)

**Outputs:** `results/lasso_cox/`, `results/nomogram/`, `results/shap/`, `results/immune/`, `results/tmb/`, `results/drug_sensitivity/`

---

### Step 6: Prognostic Model Figure

```r
source("pipeline/12_prognostic_model_figure.R")
```

**Requires:** Step 1 objects + real survival data (`OS_time`, `OS_status`). No simulated fallback.

**Does:** 6-panel composite figure:
- A: Univariate Cox forest plot
- B: ML pipeline heatmap (7 FS methods x 7 classifiers = 49 pipelines, C-index)
- C-D: KM survival curves (training + validation)
- E-F: Time-dependent ROC (1/3/5 yr, training + validation)

**Outputs:** `results/prognostic_figure/Figure_prognostic_model.{png,pdf}`, `pipeline_comparison.csv`

---

## Quick Start

```r
setwd("/path/to/BRCA_piRNA")

# Step 0: Prepare clinical data (run once)
source("scripts/prepare_clinical_data.R")  # or run the code in Step 0 above

# Step 1: Main pipeline (MUST run first)
source("piRNA_multicohort_pipeline.R")

# Save combat_df_all so later scripts can reload it
saveRDS(combat_df_all, "results/models/combat_df_all.rds")

# Steps 2-4: Can run in any order after Step 1
source("piRNA_downstream_analysis.R")
source("piRNA_functional_meta_analysis.R")   # needs mRNA data for full output
source("piRNA_network_analysis.R")

# Step 5: Advanced analyses
source("piRNA_advanced_analysis.R")          # needs mRNA + optionally TMB data

# Step 6: Prognostic model figure
source("pipeline/12_prognostic_model_figure.R")
```

---

## What Still Needs to Be Done (Action Items)

### Must do before running

1. **Upload the 5 missing expression files** to `processed_results/`:
   - `BRCA1_processed.csv` (full version, not the subset)
   - `PRJNA294226_processed.csv`
   - `PRJNA482141_processed.csv`
   - `PRJNA808405_processed.csv`
   - `PRJNA934049_processed.csv`

2. **Reformat `TCGA_BRCA_clinical.csv`** using TCGAbiolinks to add `SampleID`, `Age`, `Subtype`, `OS_time`, `OS_status` (see Step 0 above)

3. **Rename `AJCC_Stage` -> `Stage`** in yyfbatch clinical files (see Step 0)

### Recommended (for full downstream analysis)

4. **Download TCGA-BRCA mRNA expression** and place in `mRNA_expression/TCGA_BRCA_mRNA.csv`
   - Enables: piRNA-gene correlation, KEGG/GO/Reactome enrichment, immune ssGSEA, drug sensitivity prediction

### Optional

5. **Download TCGA-BRCA TMB data** and place in `clinical_data/tmb_data.csv`
   - Enables: TMB-risk group association analysis

6. **Drug IC50**: No file to download. Will be computed automatically by `oncoPredict` once mRNA data is available

---

## Directory Structure

```
BRCA_piRNA/
|
+-- processed_results/                    # INPUT: piRNA expression matrices
|   +-- BRCA1_processed.csv             # TCGA-BRCA (upload full version)
|   +-- PRJNA294226_processed.csv       # (upload)
|   +-- PRJNA482141_processed.csv       # (upload)
|   +-- PRJNA808405_processed.csv       # (upload)
|   +-- PRJNA934049_processed.csv       # (upload)
|   +-- yyfbatch1_processed.csv         # READY
|   +-- yyfbatch2_processed.csv         # READY
|
+-- clinical_data/                        # INPUT: clinical/phenotype data
|   +-- TCGA_BRCA_clinical.csv          # NEEDS PREP (reformat columns)
|   +-- BRCA_pheno.csv                  # READY (sample-group mapping)
|   +-- yyfbatch1_clinical_clean.csv    # NEEDS PREP (rename AJCC_Stage)
|   +-- yyfbatch2_clinical_clean.csv    # NEEDS PREP (rename AJCC_Stage, no survival)
|   +-- tmb_data.csv                    # (optional, download from TCGA)
|
+-- mRNA_expression/                      # INPUT: matched mRNA data
|   +-- TCGA_BRCA_mRNA.csv             # (recommended, download from TCGA)
|
+-- piRNA_multicohort_pipeline.R          # Step 1: main diagnostic pipeline
+-- piRNA_downstream_analysis.R           # Step 2: Cox, KM, subgroup
+-- piRNA_functional_meta_analysis.R      # Step 3: meta-analysis, enrichment
+-- piRNA_network_analysis.R              # Step 4: network visualization
+-- piRNA_advanced_analysis.R             # Step 5: nomogram, SHAP, immune, TMB, drug
+-- pipeline/
|   +-- 00_config.R ... 11_select_best_independent_cohort.R
|   +-- 12_prognostic_model_figure.R      # Step 6: prognostic model figure
|
+-- results/                              # OUTPUT: all generated results
    +-- models/, validation/, feature_selection/
    +-- cox_analysis/, km_curves/, subgroup_roc/, downstream/
    +-- meta_analysis/, functional/, network/
    +-- lasso_cox/, nomogram/, shap/
    +-- immune/, tmb/, drug_sensitivity/
    +-- prognostic_figure/
```

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `Error: No data found` | Place all 7 `*_processed.csv` files in `processed_results/` |
| `combat_df_all not found` | Run Step 1 first, or `combat_df_all <- readRDS("results/models/combat_df_all.rds")` |
| `OS_time`/`OS_status` missing | Reformat TCGA clinical (Step 0). Only TCGA has usable survival data |
| `Age`/`Stage`/`Subtype` missing | Reformat clinical files (Step 0). Subtype partially empty in yyfbatch — OK |
| `ComBat error` | All 7 datasets must share common piRNA column names |
| mRNA / immune / enrichment skipped | Provide `mRNA_expression/TCGA_BRCA_mRNA.csv` |
| TMB section skipped | Provide `clinical_data/tmb_data.csv` |
| Drug sensitivity placeholder | Install `oncoPredict` + provide mRNA data |
| Package install fails | For Bioconductor: `BiocManager::install("package_name")` |

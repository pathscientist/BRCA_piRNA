# MINI Workflow for Breast Cancer piRNA Diagnosis

## Overview

This MINI workflow is a simplified study design modeled after a discovery → hold-out validation → independent validation framework similar in spirit to Feng et al. It uses:

- **Five public datasets** plus **yyfbatch2** as the model-development pool
- **yyfbatch1** as the **primary in-house independent validation cohort**
- **Random Forest** as the **only machine-learning model**
- **TCGA-BRCA / BRCA1 only** for downstream survival analysis

## Cohort assignment

### Public tissue datasets
- BRCA1 (TCGA-BRCA)
- PRJNA482141

### Public plasma datasets
- PRJNA294226
- PRJNA808405
- PRJNA934049

### In-house plasma datasets
- yyfbatch1
- yyfbatch2

## MINI study design

### Primary design
- **Development pool for 70/30 split:**
  - BRCA1
  - PRJNA482141
  - PRJNA294226
  - PRJNA808405
  - PRJNA934049
  - yyfbatch2
- **Primary independent validation cohort:**
  - yyfbatch1

### Batch-effect harmonization
ComBat batch correction is performed **together across all public datasets plus yyfbatch1 and yyfbatch2**, followed by Z-score normalization.

This means:
- **yyfbatch1 is included in batch-effect harmonization**
- **yyfbatch1 is not included in the 70/30 development split**
- **yyfbatch2 is included in the development pool and enters the 70/30 split**

## Workflow steps

### Step 1. Data harmonization
1. Load all seven datasets
2. Recode group labels to unified diagnostic labels
3. Intersect common piRNAs across cohorts
4. Apply log2 transformation
5. Perform ComBat batch correction across:
   - BRCA1
   - PRJNA482141
   - PRJNA294226
   - PRJNA808405
   - PRJNA934049
   - yyfbatch1
   - yyfbatch2
6. Perform Z-score normalization

### Step 2. Build development dataset
Use the following cohorts as the development pool:
- BRCA1
- PRJNA482141
- PRJNA294226
- PRJNA808405
- PRJNA934049
- yyfbatch2

### Step 3. Split development data
Perform a **70% / 30% split** of the development pool.

Recommended split rule:
- Stratify by **Dataset** and **Group**
- Preserve tumor/normal balance where possible

Result:
- **70% Discovery / Training set**
- **30% Hold-out validation set**

### Step 4. Model construction
Use **Random Forest only**.

- No multi-model comparison in the MINI workflow
- Random Forest is trained on the 70% discovery/training subset
- Feature selection can still be constrained to a compact final signature, but the prediction model remains Random Forest only

### Step 5. Validation structure
Evaluate the model on:
1. Discovery/training performance
2. 30% hold-out validation performance
3. **yyfbatch1 independent validation performance**

## Required evaluation outputs

### Main diagnostic evaluation
Report:
- AUROC
- AUPRC
- Sensitivity
- Specificity
- Accuracy
- Balanced accuracy
- PPV
- NPV
- F1 score
- Threshold used

### Combined ROC and PRC figures
Create combined overlay plots showing multiple groups together.

#### Main figure panels
- Discovery/training
- Hold-out validation
- yyfbatch1 independent validation
- Tissue pooled samples
- Plasma pooled samples

#### Supplementary figure panels
- Individual dataset ROC overlays
- Individual dataset PRC overlays

### Grouped reporting
Report diagnostic performance for:
- Discovery/training
- Hold-out validation
- yyfbatch1 independent validation
- Tissue pooled samples
- Plasma pooled samples
- Individual datasets

### Confusion matrix
- **Optional / supplementary only**

## Tissue and plasma grouping

### Tissue datasets
- BRCA1
- PRJNA482141

### Plasma datasets
- PRJNA294226
- PRJNA808405
- PRJNA934049
- yyfbatch1
- yyfbatch2

For grouped evaluation, show:
- Tissue-pooled ROC/PRC
- Plasma-pooled ROC/PRC
- Independent validation ROC/PRC for yyfbatch1

## Survival analysis
Use **TCGA-BRCA / BRCA1 only**.

Do **not** use yyfbatch1 or yyfbatch2 for survival analysis.

### Allowed downstream survival outputs
- Cox regression
- Kaplan–Meier survival analysis
- Time-dependent ROC
- Nomogram / prognostic figure if needed

## Suggested interpretation
This MINI workflow is designed to be close to a journal-style study design:

- Public + in-house development pool for discovery and hold-out validation
- One in-house cohort reserved as the primary independent endpoint cohort
- One consistent classifier (Random Forest)
- Compact grouped ROC/PRC reporting
- TCGA-only survival analysis

## Summary statement

**Final MINI design:**
- **yyfbatch2 is included in the 70/30 development split**
- **yyfbatch1 is the primary independent validation cohort**
- **ComBat harmonization includes all public datasets plus yyfbatch1 and yyfbatch2**
- **Random Forest is the only machine-learning model**
- **TCGA-BRCA / BRCA1 only is used for survival analysis**

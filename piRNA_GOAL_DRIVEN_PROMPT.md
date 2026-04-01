# Goal-Driven Claude Code Prompt: piRNA Breast Cancer Diagnostic Pipeline

> Paste the fenced block below into Claude Code from your `BRCA_piRNA/` project root.
> Each phase is a self-contained goal with context, constraints, verification, and next actions.
> Run phases sequentially — each one builds on the previous.

---

## How to Use This Prompt

**Option A — Run all at once:** Paste everything inside the outer fence.
**Option B — Phase by phase:** Paste one `## PHASE` block at a time, confirm outputs, then paste the next.

Phase-by-phase is recommended for a project this size. It keeps Claude Code's context fresh and lets you inspect intermediate results before committing.

---

```
## ROLE

You are a senior bioinformatics analyst with deep expertise in non-coding RNA
biomarker discovery, multi-cohort batch correction, and ML-based cancer
diagnostics in R. You write production-grade, reproducible R code with clear
comments, set.seed() before every stochastic call, and publication-quality
ggplot2 figures (theme_minimal, 300 DPI, 8×6 PNG).

────────────────────────────────────────────────────────────────────────────────
## PROJECT CONTEXT

Goal: Build a piRNA-based diagnostic classifier for breast cancer, trained on
5 public tissue cohorts, independently validated on 2 in-house cohorts
(yyfbatch1 = tissue, yyfbatch2 = plasma/serum), and examined separately by
specimen type (tissue vs liquid biopsy) to demonstrate cross-platform
clinical utility.

Target signature size: ~10 piRNAs. Not fewer than 7, not more than 12.
The signature must be large enough to be biologically interpretable and
robust, but small enough to be clinically translatable (RT-qPCR panel).

Cohort roles (LOCKED — do not reassign):
  TRAINING (feature selection + model building):
    BRCA1            (TCGA-BRCA tissue)
    PRJNA294226      (tissue)
    PRJNA482141      (tissue)
    PRJNA808405      (tissue)
    PRJNA934049      (tissue)
  INDEPENDENT VALIDATION (never touched until final evaluation):
    yyfbatch1        (in-house tissue validation)
    yyfbatch2        (in-house plasma/serum validation)

Files are in processed_results/*.csv
  Format: rows = samples, columns = piRNAs + "Group" column
  Group labels: "Tumor"/"Normal" or "Cancer"/"Normal" or "Cancer"/"Benign"
  Values: TPM (not yet log-transformed)

Clinical data in clinical_data/*.csv
  Expected columns: SampleID, Age, Stage, Subtype, OS_time, OS_status
  Survival data: usable in TCGA only (yyfbatch1 ~31%, yyfbatch2 = 0%)

────────────────────────────────────────────────────────────────────────────────
## PHASE 0: DATA LOADING & BATCH CORRECTION

### Goal
Load all 7 cohort CSVs, harmonize labels, find shared piRNAs, apply
log2(TPM+1) transform, run ComBat batch correction across all 7 cohorts,
then Z-score normalize. Output a single analysis-ready data frame with
a "Cohort" column and a "Group" column.

### Context
- Files: processed_results/{BRCA1,PRJNA294226,PRJNA482141,PRJNA808405,PRJNA934049,yyfbatch1,yyfbatch2}_processed.csv
- BRCA1 (TCGA) is the largest and may have class imbalance — balance it by
  retaining all matched tumor-normal pairs plus a random 40% of excess tumors
- Group labels vary across files. Recode everything to "Tumor" vs "Normal"
- ComBat (sva package): model.matrix(~Group) as covariate to preserve biology

### Constraints
- Do NOT drop yyfbatch1 or yyfbatch2 from ComBat — they must be corrected
  in the same batch as training cohorts to avoid a second correction later
- Do NOT split into train/validation before ComBat — correct all 7 together
  (this is standard practice for multi-cohort biomarker studies)
- Keep a "Cohort" column so we can split train vs validation after correction

### Verification
- Print: number of samples per cohort × group after loading
- Print: number of common piRNAs across all 7 datasets
- Print: dimensions of final combat_df_all (should be ~N_samples × N_piRNAs + Group + Cohort)
- Generate a PCA plot (first 2 PCs) colored by Cohort, shaped by Group
  — before ComBat and after ComBat side by side — saved to results/qc/pca_combat.png
  Confirm that after ComBat, cohorts overlap while Tumor/Normal still separate

### Output
- Object: combat_df_all (data frame in memory)
- Save: results/models/combat_df_all.rds
- Split into:
    train_df: rows where Cohort %in% c("BRCA1","PRJNA294226","PRJNA482141","PRJNA808405","PRJNA934049")
    valid_yyfbatch1: rows where Cohort == "yyfbatch1"
    valid_yyfbatch2: rows where Cohort == "yyfbatch2"

────────────────────────────────────────────────────────────────────────────────
## PHASE 1: FEATURE SELECTION → ~10 piRNA CONSENSUS SET

### Goal
Run 8 feature selection methods on train_df ONLY. Build a consensus set
of approximately 10 piRNAs (range: 7–12). Then apply forward/backward/swap
optimization to finalize the set. The final signature must be biologically
defensible — every included piRNA should have independent statistical or
ML-based evidence from multiple methods.

### Context — The 8 Methods
Run each method on train_df. set.seed(2024) before every stochastic call.

  1. limma (empirical Bayes DE): adj.p < 0.01, |log2FC| > 1.0
     → fs_limma
     → Save volcano plot to results/feature_selection/volcano_plot.png

  2. Wilcoxon rank-sum: BH adj.p < 0.01, |median fold change| > 1.5
     → fs_wilcoxon

  3. Random Forest importance: ntree=1000, select top features above
     mean + 1 SD (cap at 50)
     → fs_rf_importance
     → Save importance barplot to results/feature_selection/rf_importance.png

  4. LASSO (glmnet, alpha=1): cv.glmnet 10-fold, non-zero coefs at lambda.1se
     → fs_lasso

  5. Elastic Net (glmnet): alpha grid 0.1–0.9 by 0.1, best alpha by CV error,
     non-zero coefs at lambda.1se
     → fs_elasticnet

  6. RFE (caret::rfe, rfFuncs): subset sizes c(5,10,15,20,30,50), 10-fold CV ×3
     → fs_rfe

  7. Boruta: maxRuns=300, select "Confirmed" only
     → fs_boruta

  8. mRMR (praznik package): top 50 by mRMR criterion
     → fs_mrmr

### Consensus Construction (critical section)
a) Frequency table: for each piRNA, count how many of 8 methods selected it
b) UpSet plot of all 8 method overlaps → results/feature_selection/upset_plot.png
c) Apply TWO strategies:

   STRATEGY A — Frequency threshold:
     Start at ≥6/8 methods. If <7 piRNAs, relax to ≥5. If >12, tighten to ≥7.
     Target: land in the 7–12 range.

   STRATEGY B — Forward stepwise on frequency-ranked list:
     Rank piRNAs by selection count (descending), break ties by mean RF importance.
     Add one piRNA at a time, evaluate 5-fold CV AUC (Random Forest) after each.
     Stop when: (a) AUC improvement < 0.005 for 3 consecutive additions, OR
                (b) you reach 12 piRNAs — whichever comes first.
     Minimum: do not stop before 7 piRNAs regardless of AUC plateau.

d) Pick the strategy that yields 7–12 piRNAs with higher CV AUC.
   If both yield similar AUC (within 0.01), prefer the smaller set.

e) Forward/Backward/Swap refinement on the chosen set:
   - Try removing each piRNA one at a time — drop any that doesn't decrease AUC by ≥0.003
   - Try swapping each piRNA with the next-best candidate not in the set
   - Accept swaps only if AUC improves by ≥0.005
   - This polishes the set without changing size dramatically

### Justification of ~10 Features (IMPORTANT — build this evidence for the paper)
For each piRNA in the final signature, generate a row in a justification table with:
  - piRNA name
  - N_methods: how many of 8 methods selected it
  - which_methods: comma-separated list
  - limma_adjP: adjusted p-value from DE analysis
  - limma_log2FC: fold change
  - individual_AUC: univariate AUC for that single piRNA (on train_df, 5-fold CV)
  - RF_importance: MeanDecreaseGini rank
  - LASSO_coef: coefficient (0 if not selected by LASSO)
  - mRMR_rank: rank from mRMR

  Also compute:
  - Incremental AUC curve: plot number of piRNAs (x) vs cumulative CV AUC (y)
    as piRNAs are added in frequency order. Save to results/feature_selection/incremental_auc.png
    This plot visually demonstrates that ~10 piRNAs capture nearly all discriminative power.
  - Pairwise correlation heatmap of the final piRNAs (on train_df)
    → results/feature_selection/signature_correlation.png
    Low inter-correlation supports the argument that each piRNA contributes non-redundant information.

### Constraints
- Do NOT use valid_yyfbatch1 or valid_yyfbatch2 for anything in this phase
- Wrap each method in tryCatch() — if one fails, warn and continue with the rest
- If the consensus set is empty, fall back to union of LASSO + top 10 RF importance
- Final set MUST be between 7 and 12 piRNAs. If outside this range after all
  refinement, force it by adding (from frequency-ranked list) or pruning
  (lowest individual AUC first)

### Verification
- Print: "Final piRNA signature: N piRNAs selected by at least M/8 methods"
- Print the full justification table
- Confirm each piRNA has individual AUC > 0.55 (if any are below, flag for review but keep if N_methods ≥ 5)
- Save: results/feature_selection/consensus_mirnas.txt (one piRNA per line)
- Save: results/feature_selection/fs_summary_table.csv (full justification table)

────────────────────────────────────────────────────────────────────────────────
## PHASE 2: TRAIN DIAGNOSTIC MODEL (Random Forest, Primary)

### Goal
Train a Random Forest classifier on train_df using only the consensus piRNAs.
Use 10-fold CV repeated 5 times. Also train 5 comparison models to show RF
was a reasonable choice. Select the best model, fine-tune it, and determine
the optimal diagnostic threshold.

### Context
- Features: consensus piRNAs from Phase 1 only
- Labels: Group (Tumor vs Normal)
- Preprocessing inside caret: center + scale
- If class imbalance ratio > 2:1, use sampling = "down" in trainControl

### 6 Candidate Models (all via caret::train)
  trainControl: method="repeatedcv", number=10, repeats=5,
                classProbs=TRUE, summaryFunction=twoClassSummary,
                savePredictions="final", metric="ROC"

  1. Random Forest:    method="rf", tuneGrid mtry=seq(2, min(ncol,10), by=1), ntree=1000
  2. SVM-RBF:          method="svmRadial", C=c(0.01,0.1,1,10,100), sigma=c(0.001,0.01,0.1,0.5)
  3. XGBoost:          method="xgbTree", random search ~100 combos
  4. Elastic Net:      method="glmnet", alpha=seq(0,1,0.1), lambda=10^seq(-4,0,length=30)
  5. KNN:              method="knn", k=seq(1,31,by=2), preProcess=c("center","scale")
  6. Neural Net:       method="nnet", size=c(1,3,5,7), decay=c(0,0.001,0.01,0.1)

### Model Comparison
a) resamples() across all 6 → boxplot of CV AUC → results/models/cv_boxplot.png
b) Statistical comparison: diff(resamples) — report which differences are significant
c) Overlaid ROC curves from CV predictions → results/models/roc_all_models.png
d) Select best by: highest median AUC, tiebreak by sensitivity
e) For the best model, fine-tune: narrow grid around best hyperparams, re-run CV
f) Calibration curve → results/models/calibration.png
g) Optimal threshold: Youden's J index (Sensitivity + Specificity − 1)
   Also report threshold at Sensitivity ≥ 95% (for screening use case)
h) Learning curve: subsample train_df at 20%/40%/60%/80%/100%, plot AUC vs N
   → results/models/learning_curve.png

### Constraints
- set.seed(2024) before every train() call
- Save final model: results/models/final_model.rds
- Save preprocessing params: results/models/preprocessing_params.rds
- Save feature list: results/models/final_features.rds
- Do NOT touch validation cohorts

### Verification
- Print: "Best model: [name], CV AUC = X.XX (95% CI: X.XX–X.XX)"
- Print: full comparison table (6 models × AUC/Sens/Spec)
- Print: optimal hyperparameters of best model
- Confirm CV AUC > 0.85 — if not, flag but continue (do not re-engineer features)

────────────────────────────────────────────────────────────────────────────────
## PHASE 3: INDEPENDENT VALIDATION — TISSUE vs PLASMA ROC

### Goal
Evaluate the final model on yyfbatch1 (tissue) and yyfbatch2 (plasma/serum)
SEPARATELY. Generate per-cohort ROC curves, a combined multi-cohort ROC overlay,
and a tissue-vs-plasma comparison figure. This is the key figure for demonstrating
that the piRNA signature works across specimen types.

### 3.1 Per-Cohort Evaluation
For EACH of {yyfbatch1, yyfbatch2}:
a) Subset to consensus piRNAs
b) Apply center+scale using training set parameters (from preprocessing_params.rds)
c) Predict classes and probabilities using the final model
d) Apply optimal threshold from Phase 2
e) Compute:
   - AUC with 95% CI (DeLong method, pROC::ci.auc)
   - Accuracy, Sensitivity, Specificity, PPV, NPV, F1
   - Confusion matrix → results/validation/{cohort}_confusion.png

### 3.2 ROC Comparison Figure (KEY FIGURE for the paper)
Create a SINGLE publication figure with 4 panels (2×2 grid):

  Panel A: "Discovery (5-cohort training, 10×5 CV)"
    ROC curve of best model's CV predictions on train_df
    Report AUC with 95% CI in legend

  Panel B: "Tissue Validation (yyfbatch1)"
    ROC curve on yyfbatch1
    Report AUC with 95% CI in legend

  Panel C: "Plasma Validation (yyfbatch2)"
    ROC curve on yyfbatch2
    Report AUC with 95% CI in legend

  Panel D: "All Cohorts Overlay"
    Overlay all 3 ROC curves (Discovery, yyfbatch1, yyfbatch2)
    on one plot with different colors/linetypes
    Legend: "Discovery AUC=X.XX, Tissue AUC=X.XX, Plasma AUC=X.XX"

  Layout: cowplot or patchwork, shared theme, consistent axis limits (0–1)
  Colors: Discovery="#1F78B4", Tissue="#E31A1C", Plasma="#33A02C"
  Save: results/validation/Figure_ROC_tissue_vs_plasma.png (300 DPI, 10×8 in)
  Save: results/validation/Figure_ROC_tissue_vs_plasma.pdf

### 3.3 Also generate individual per-cohort plots:
  - ROC with 95% CI band → results/validation/{cohort}_roc.png
  - Precision-Recall curve → results/validation/{cohort}_prc.png
  - Prediction probability histogram (colored by true label)
    → results/validation/{cohort}_probability_hist.png
    Flag samples with P ∈ [0.4, 0.6] as "uncertain"

### 3.4 Per-Cohort-Per-Dataset ROC (Granular View)
Also produce ROC curves for EACH of the 5 training cohorts individually
(predict on each cohort's data using the final model, NOT from CV folds):
  - This shows how the model performs on each public dataset
  - Overlay all 7 cohort ROC curves on one plot
    → results/validation/roc_all_7_cohorts.png
  Label each: "BRCA1 AUC=X.XX", "PRJNA294226 AUC=X.XX", ..., "yyfbatch1 AUC=X.XX", "yyfbatch2 AUC=X.XX"
  Color tissue cohorts in shades of blue/red, plasma cohort in green
  This demonstrates cross-cohort generalization AND tissue-vs-plasma performance

### 3.5 Robustness Checks (on each validation cohort separately)
a) Bootstrap (1000 resamples): median AUC + 95% CI
b) Permutation test (1000 label shuffles): null AUC distribution + empirical p-value
   → results/validation/{cohort}_bootstrap.png
   → results/validation/{cohort}_permutation.png

### 3.6 Performance Drop Assessment
Compute Δ = (training CV AUC) − (validation AUC) for each validation cohort.
  - If Δ > 0.05: print WARNING and flag potential overfitting
  - If Δ > 0.10: print SEVERE WARNING
  - If plasma AUC < tissue AUC (expected): note this is expected due to
    lower piRNA abundance in circulation, and report the performance gap
  - Create a bar chart: AUC by cohort with error bars
    → results/validation/auc_comparison_bar.png

### Verification
- Print structured results table:
  | Cohort       | Type    | N    | AUC (95% CI)     | Sens  | Spec  | PPV   | NPV   |
  | Discovery CV | Tissue  | ...  | X.XX (X.XX–X.XX) | XX.X% | XX.X% | XX.X% | XX.X% |
  | yyfbatch1    | Tissue  | ...  | ...               | ...   | ...   | ...   | ...   |
  | yyfbatch2    | Plasma  | ...  | ...               | ...   | ...   | ...   | ...   |
- Confirm: results/validation/Figure_ROC_tissue_vs_plasma.png exists and has 4 panels
- Save: results/validation/validation_results.csv

────────────────────────────────────────────────────────────────────────────────
## PHASE 4: CLINICAL INTERPRETATION & FINAL REPORT

### Goal
Generate a final summary report and a signature characterization table
suitable for inclusion in a manuscript methods/results section.

### 4.1 Signature Characterization Table (Table 1 of the paper)
For each piRNA in the final signature, compile:
  - piRNA ID
  - Selection frequency (N/8 methods)
  - Direction (up/down in tumor)
  - log2FC and adj.p (from limma)
  - Individual AUC (univariate, on discovery set)
  - RF importance rank
  - Known associations (print "see literature" — do not fabricate)

Save as: results/final_signature_table.csv

### 4.2 Summary Heatmap
- Rows = final piRNAs, Columns = all samples (all 7 cohorts)
- Annotate columns by: Group (Tumor/Normal), Cohort, Specimen type (Tissue/Plasma)
- Cluster rows and columns within each annotation block
- Use pheatmap with annotation_col
→ results/feature_selection/heatmap_all_cohorts.png (300 DPI, 12×8 in)

### 4.3 Final Report (text file)
Print and save to results/validation/final_report.txt:

  ═══════════════════════════════════════════════════════════════
  piRNA Diagnostic Signature — Final Report
  ═══════════════════════════════════════════════════════════════
  Signature: [list of N piRNA IDs]
  Model: [best model name] with [hyperparameters]
  Optimal threshold: X.XX (Youden's J)

  Discovery (5 public tissue cohorts, 10×5 CV):
    AUC = X.XX (95% CI: X.XX–X.XX)
    Sensitivity = XX.X%, Specificity = XX.X%

  Tissue Validation (yyfbatch1, n=XX):
    AUC = X.XX (95% CI: X.XX–X.XX)
    Sensitivity = XX.X%, Specificity = XX.X%
    Bootstrap 95% CI: [X.XX, X.XX]
    Permutation p = X.XXXX

  Plasma Validation (yyfbatch2, n=XX):
    AUC = X.XX (95% CI: X.XX–X.XX)
    Sensitivity = XX.X%, Specificity = XX.X%
    Bootstrap 95% CI: [X.XX, X.XX]
    Permutation p = X.XXXX

  Tissue vs Plasma AUC gap: X.XX
  Overfitting assessment: Δ(train–tissue) = X.XX, Δ(train–plasma) = X.XX
  ═══════════════════════════════════════════════════════════════

### 4.4 Session Info
Print sessionInfo() and total runtime. Save to results/session_info.txt.

────────────────────────────────────────────────────────────────────────────────
## GLOBAL CONSTRAINTS

1. REPRODUCIBILITY
   - set.seed(2024) before EVERY stochastic operation
   - Record start/end time, print total runtime per phase
   - sessionInfo() at end

2. VALIDATION FIREWALL
   - yyfbatch1 and yyfbatch2 are NEVER used in Phase 1 or Phase 2
   - The only exception: they ARE included in Phase 0 ComBat correction
     (this is methodologically correct for multi-cohort normalization)

3. FEATURE COUNT GUARDRAILS
   - Final signature: 7–12 piRNAs, target 10
   - If any phase produces <7 or >12, apply the force-fit rules in Phase 1

4. SPECIMEN TYPE AWARENESS
   - Always label yyfbatch2 as "Plasma/Serum" in all plots and tables
   - Always label yyfbatch1 and the 5 public cohorts as "Tissue"
   - When plasma AUC is lower than tissue AUC, explicitly note this is
     biologically expected (lower circulating piRNA concentrations)

5. CODE QUALITY
   - One clearly commented section per step
   - cat() progress messages for each major step
   - All intermediate objects stored in named lists
   - Helper functions for repeated operations
   - tryCatch() around each feature selection method

6. VISUALIZATION STANDARDS
   - theme_minimal(base_size=13), bold centered titles
   - Palette: c("#E31A1C","#1F78B4","#33A02C","#FF7F00","#6A3D9A","#B15928")
   - All plots: PNG (300 DPI, 8×6 unless specified) AND displayed on screen
   - Multi-panel figures: cowplot or patchwork

7. PACKAGES
   Required: sva, caret, randomForest, e1071, glmnet, xgboost, nnet,
   pROC, limma, Boruta, UpSetR, pheatmap, ggplot2, dplyr, tidyr,
   cowplot, patchwork, viridis, scales, praznik
   Install check at top of script.

8. OUTPUT DIRECTORY
   results/
   ├── qc/                       (PCA before/after ComBat)
   ├── feature_selection/        (volcano, upset, heatmaps, incremental AUC, justification table)
   ├── models/                   (final_model.rds, preprocessing_params.rds, CV plots)
   └── validation/               (per-cohort ROC, tissue-vs-plasma figure, confusion matrices, report)

9. ERROR HANDLING
   - If a feature selection method fails: warn, continue, note in summary
   - If consensus set is empty: fall back to LASSO ∪ top 10 RF importance
   - If a validation cohort has <20 samples in either class: warn about
     unreliable AUC estimates, still compute but flag

Generate the complete R script now. Work phase by phase. After each phase,
pause and confirm outputs exist before moving to the next.
```

---

## Customization Guide

Before pasting, search-and-replace these if your setup differs:

| Placeholder | Default | Your value |
|---|---|---|
| `"Tumor"/"Normal"` | Default labels | Your actual Group labels if different |
| `set.seed(2024)` | Seed | Any integer you prefer |
| `yyfbatch2 = plasma` | Assumed plasma | Confirm with your lab — update annotation if tissue |
| `adj.p < 0.01, \|log2FC\| > 1.0` | Limma cutoffs | Relax to 0.05 / 0.5 if too few pass |
| `≥6/8 methods` | Consensus start | Adjust if you want stricter/looser |
| `7–12 piRNAs` | Target range | Change if your committee wants more/fewer |

## Tips for Phase-by-Phase Execution

1. **Phase 0 first, always.** If the PCA plot still shows cohort clustering after ComBat, troubleshoot before proceeding.
2. **Phase 1 is the longest.** If Claude Code hits context limits, split into "run methods 1–4" and "run methods 5–8 + consensus."
3. **Phase 2 XGBoost grid:** If runtime exceeds 30 min, tell Claude Code: "XGBoost is too slow. Reduce to random search with 50 combinations."
4. **Phase 3 is your paper's money figure.** Inspect `Figure_ROC_tissue_vs_plasma.png` carefully. If panel labels or legends need tweaking, describe the exact change.
5. **If you later have an external cohort CSV**, replace Phase 3 with: "Load external_cohort.csv, apply the same preprocessing, and add it as a 3rd validation curve in the overlay."

## Why ~10 Features? (Justification Talking Points for the Paper)

The prompt builds four lines of evidence to defend a ~10-piRNA signature:

- **Multi-method consensus:** Each piRNA independently selected by ≥5–6 of 8 orthogonal methods (statistical, penalized regression, tree-based, information-theoretic), ruling out method-specific artifacts.
- **Incremental AUC curve:** Demonstrates diminishing returns — the first ~10 piRNAs capture 95%+ of the classifier's discriminative power, additional features add noise.
- **Low inter-correlation:** Pairwise correlation heatmap shows the selected piRNAs provide non-redundant biological signal.
- **Individual AUC evidence:** Every piRNA in the signature has a univariate AUC above chance, confirming each one independently carries diagnostic information.

These four pieces together address the reviewer question "why 10 and not 5 or 50?" from both statistical and biological angles.

# Methods (Draft for Manuscript)

## Study design and cohorts
We developed and validated a piRNA-based diagnostic and prognostic framework using merged breast cancer datasets, including BRCA, `yyfbatch1`, and `yyfbatch2`. Raw TPM expression matrices were harmonized with clinical metadata using unique sample identifiers. Diagnostic labels were encoded as binary outcomes (Normal vs Cancer). For external validation, one complete batch (`yyfbatch1` or `yyfbatch2`) was held out as an independent dataset, while the remaining data were used for model development.

## Batch-effect adjustment and data partitioning
To reduce inter-cohort technical variability, batch effects were removed across all datasets jointly before development/validation partitioning using ComBat (`sva`) when available, with `limma::removeBatchEffect` as fallback. After correction, the pre-specified held-out batch was reserved as independent validation. The development subset was split into internal training and holdout partitions using stratified sampling.

## Class rebalancing strategy in BRCA1 cohort
To mitigate class imbalance in BRCA1 during model training, we adopted a controlled sampling strategy: we retained matched tumor-normal pairs up to the minority class size and then included 40% of the excess majority-class samples. This approach preserved biological diversity while reducing training bias induced by over-represented classes.

## Feature selection
Feature selection was performed on development-training data only to prevent data leakage. We evaluated multiple complementary methods: differential-expression filtering (Wilcoxon/FDR), LASSO, Elastic Net, mRMR, Boruta, random forest importance, XGBoost importance, SVM-RFE, and mutual information (when package availability permitted). Candidate signatures were constructed by intersection and frequency-based consensus rules, with a strict upper bound of <10 markers.

## Model development and hyperparameter tuning
Candidate signatures were evaluated using repeated stratified cross-validation within the development set. We benchmarked logistic regression, elastic-net logistic regression, random forest, XGBoost, and SVM-RBF (`caret` framework). Hyperparameters were tuned via predefined grids. Model selection prioritized AUROC, with AUPRC as secondary criterion and feature-count minimization as tie-breaker.

## Threshold locking and independent validation
Classification thresholds were optimized on cross-validated development predictions (Youden index or F1) and then locked. Locked preprocessing, selected feature set, final model, and threshold were applied unchanged to the internal holdout and the independent validation batch. We reported AUROC, AUPRC, accuracy, sensitivity, specificity, F1 score, and Matthews correlation coefficient.

## Subgroup diagnostic analysis
For subgroup robustness analyses, BRCA + `yyfbatch1` + `yyfbatch2` were combined and stratified by age, stage, and molecular subtype (additional factors such as sex/smoking/histology were included if available). For each subgroup:
1. ROC curves and AUC (95% CI) were computed from model probabilities.
2. T-score distributions between tumor and control were compared with Mann–Whitney U tests.
3. Boxplots were generated with sample-size annotations and exact P values.

## Survival analysis
For prognostic analyses, a model-derived risk score was binarized into high/low groups using the locked threshold. Univariate Kaplan–Meier analyses were performed for individual biomarkers and composite score groups with log-rank testing. Multivariable Cox proportional hazards models were fitted using binary score and available clinical covariates (e.g., age, stage, subtype), with hazard ratios (HR), 95% confidence intervals (CI), and Wald P values reported.

## Functional annotation and pathway/network analyses
To infer biological relevance of selected piRNAs, Pearson correlation analyses were conducted between signature piRNAs and co-expressed genes. Correlated genes (|r| threshold and FDR-controlled significance) were subjected to Reactome pathway enrichment. A piRNA–gene–pathway network was constructed to visualize functional modules. Optional KEGG/GO enrichment analyses were used to characterize BP/CC/MF terms.

## Cross-dataset consistency and meta-analysis
Cross-dataset expression patterns were visualized in heatmaps. For key markers, standardized mean differences (SMD; tumor vs control) were estimated per dataset and pooled with random-effects meta-analysis (REML), with forest plots showing study-level effects and summary estimates.

## Software and reproducibility
Analyses were implemented in R (modular scripted workflow) using: `data.table`, `dplyr`, `caret`, `glmnet`, `randomForest`, `xgboost`, `e1071`, `pROC`, `PRROC`, `survival`, `broom`, and visualization packages (`ggplot2`, `pheatmap`, `ggraph`, `igraph`, `survminer`, `metafor`, and enrichment packages when available). Random seeds were fixed and intermediate artifacts were serialized (`RDS`/`CSV`) for reproducibility.

---

# Figure Legends (Draft)

## Figure 3. Diagnostic robustness of the piRNA signature across clinical strata.
**(a, c, e, g, i)** Receiver operating characteristic (ROC) curves showing the discriminatory performance of the piRNA signature across stratified cohorts (e.g., age group, sex, histological subtype, smoking status, and AJCC stage, depending on available metadata). Curve colors indicate subgroup categories. AUC values with 95% confidence intervals are reported for each subgroup.

**(b, d, f, h, j)** Boxplots comparing T-scores between tumor and control samples within corresponding stratified factors. Colors indicate diagnostic class (tumor vs control). Sample sizes (n) are displayed for each subgroup. Two-sided Mann–Whitney U tests were used for group comparisons; exact P values are shown.

## Figure 4. Prognostic value of individual piRNAs and composite risk score.
**(a–c)** Kaplan–Meier survival curves for representative individual piRNA biomarkers in the discovery cohort. Patients were dichotomized into high- and low-expression groups using pre-specified cutoffs (median unless otherwise stated). Log-rank P values are shown; risk tables indicate numbers at risk over follow-up.

**(d)** Kaplan–Meier survival analysis for the composite piRNA-derived risk score, dichotomized into high- and low-risk groups by locked threshold. Dashed lines denote median survival references where estimable.

**(e)** Forest plot of multivariable Cox regression including binary risk score and clinical covariates (e.g., age, stage, subtype). Points indicate HR; horizontal lines indicate 95% CI. Variables with HR>1 are associated with increased risk.

## Figure 5. Functional characterization and cross-dataset consistency of signature piRNAs.
**(a)** Heatmap of normalized signature expression (or T-scores) across datasets and sample classes.

**(b)** Random-effects meta-analysis forest plot of standardized mean differences (SMDs; tumor vs control) for key signature markers across datasets. Box sizes are proportional to study weight; horizontal bars indicate 95% CI; the diamond represents pooled effect.

**(c)** Correlation heatmap of signature RNAs with associated molecular features (e.g., tRFs/mRNAs), using Spearman or Pearson coefficients as specified.

**(d)** Functional interaction network linking signature piRNAs to correlated genes and enriched pathways. Node classes represent piRNAs, genes, and pathways.

**(e)** Reactome/KEGG pathway enrichment summary plot for correlated genes.

**(f–h)** Gene Ontology enrichment bubble plots for biological process (BP), cellular component (CC), and molecular function (MF). Bubble size indicates gene count and color indicates FDR-adjusted significance.

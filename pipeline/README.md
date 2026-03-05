# Breast Cancer miRNA TPM ML Pipeline (R)

This folder contains a leakage-aware template pipeline for:

1. multi-method feature selection (8+ methods),
2. model tuning and comparison,
3. locked independent validation,
4. subgroup ROC/boxplot analysis and binary-score Cox regression,
5. functional downstream analyses (Pearson/Reactome/network, KM, meta-analysis templates).

## Script order

- `00_config.R` - configuration, package checks, helper metrics, threshold optimization, feature cap (`<10`).
- `01_load_qc.R` - data loading, merge by sample ID, global batch correction across all datasets, locked validation-batch split (`yyfbatch1` or `yyfbatch2`), BRCA1 imbalance rebalancing.
- `02_feature_selection.R` - DE, LASSO, Elastic Net, mRMR, Boruta, RF, XGBoost, SVM-RFE, MI.
- `03_model_training.R` - repeated stratified CV training and tuning for logistic / EN / RF / XGB / SVM-RBF.
- `04_model_selection.R` - best model selection using AUROC primary, AUPRC secondary, feature-count tie-break.
- `05_external_validation.R` - locked threshold + locked feature set evaluation on holdout/external sets.
- `06_reports_plots.R` - model ranking plot, ROC plots, final report with `sessionInfo()`.
- `07_clinical_subgroup_analysis.R` - merge BRCA + yyfbatch1 + yyfbatch2, subgroup ROC (age/stage/subtype), tumor-vs-normal T-score boxplots, binary-score Cox regression and forest plot.
- `08_downstream_biology_survival_meta.R` - templates for Pearson-correlated genes, Reactome enrichment, pathway network, stratified ROC+Mann-Whitney boxplots, KM curves, multivariable Cox, and SMD meta-analysis forest plots.

## Notes

- Fill placeholders in `00_config.R` before running.
- To avoid leakage, run feature selection and tuning only on training splits.
- Global batch correction is applied **before** selecting the locked independent validation batch.
- External validation batch is selected by `config$input$validation_batch` (`yyfbatch1` or `yyfbatch2`).
- BRCA1 balancing strategy keeps matched tumor/normal pairs and adds 40% of the excess class in training.
- Optional methods are skipped with warnings when packages are unavailable.

# Breast Cancer piRNA TPM ML Pipeline (R)

This folder contains a leakage-aware template pipeline for:

1. multi-method feature selection (8+ methods),
2. model tuning and comparison,
3. locked independent validation,
4. subgroup ROC/boxplot analysis and binary-score Cox regression,
5. functional downstream analyses (Pearson/Reactome/network, KM, meta-analysis templates),
6. dedicated gene-network centrality/module analysis and Cytoscape export,
7. advanced clinical/omics extensions (univariate Cox, LASSO-Cox, nomogram, immune infiltration, TMB, drug-sensitivity stub).

## Script order

- `00_config.R` - configuration, package checks, helper metrics, threshold optimization, feature cap (`<10`), and explicit default dataset file names (`*_processed.csv`).
- `00b_data_preflight.R` - one-command real-data preflight that validates file existence, required columns, label mapping settings, and writes pass/fail reports before running modeling scripts.
- `01_load_qc.R` - load and merge multiple processed piRNA dataset files, optional label mapping join, global batch correction, locked validation-batch split (`yyfbatch1` or `yyfbatch2`), BRCA1 imbalance rebalancing.
- `02_feature_selection.R` - DE, LASSO, Elastic Net, mRMR, Boruta, RF, XGBoost, SVM-RFE, MI.
- `03_model_training.R` - repeated stratified CV training and tuning for logistic / EN / RF / XGB / SVM-RBF.
- `04_model_selection.R` - best model selection using AUROC primary, AUPRC secondary, feature-count tie-break.
- `05_external_validation.R` - locked threshold + locked feature set evaluation on holdout/external sets.
- `06_reports_plots.R` - model ranking plot, ROC plots, final report with `sessionInfo()`.
- `07_clinical_subgroup_analysis.R` - merge BRCA + yyfbatch1 + yyfbatch2, subgroup ROC (age/stage/subtype), tumor-vs-normal T-score boxplots, binary-score Cox regression and forest plot.
- `08_downstream_biology_survival_meta.R` - templates for Pearson-correlated genes, Reactome enrichment, pathway network, stratified ROC+Mann-Whitney boxplots, KM curves, multivariable Cox, and SMD meta-analysis forest plots.
- `09_gene_network_analysis.R` - builds piRNA-gene/gene-gene networks from correlation tables, computes centrality and modules, exports Cytoscape-ready node/edge tables, and saves a network overview plot.
- `10_advanced_clinical_omics_analysis.R` - extended module for univariate Cox screening, LASSO-Cox risk score, KM + time-ROC, optional nomogram, immune infiltration scoring, TMB association, and drug-sensitivity placeholder outputs.
- `11_select_best_independent_cohort.R` - runs end-to-end training twice (yyfbatch1/yyfbatch2 as external cohort), compares external AUC, enforces AUC >= 0.80 and feature count <= 10, and exports best-cohort summary artifacts.

## Notes

- Fill placeholders in `00_config.R` before running.
- To avoid leakage, run feature selection and tuning only on training splits.
- Global batch correction is applied **before** selecting the locked independent validation batch.
- External validation batch is selected by `config$input$validation_batch` (`yyfbatch1` or `yyfbatch2`).
- BRCA1 balancing strategy keeps matched tumor/normal pairs and adds 40% of the excess class in training.
- Optional methods are skipped with warnings when packages are unavailable.

- Input files are configured in `config$input$dataset_files` (e.g., `BRCA1_processed.csv`, `PRJNA*_processed.csv`, `yyfbatch*_processed.csv`).

- Run preflight first: `Rscript pipeline/00b_data_preflight.R`.

- Independent-cohort comparison outputs: `artifacts/batch_compare_summary.csv`, `artifacts/batch_compare_best.csv`, `artifacts/batch_compare_summary.md`.

- Full command protocol: `pipeline/RSCRIPT_PROTOCOL.md`.

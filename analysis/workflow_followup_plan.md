# Follow-up analysis plan (matching workflow-style figure)

This extension maps to the requested workflow pattern:

1. **Univariate Cox screening** for candidate piRNAs.
2. **LASSO-Cox signature** construction and risk-score calculation.
3. **Kaplan-Meier + time-dependent ROC** performance evaluation.
4. **Nomogram** (optional, when `rms` is installed and survival fields are complete).
5. **Immune infiltration and immune-cell proportion** estimation (optional `immunedeconv`).
6. **TMB association** test between high/low risk groups.
7. **Drug sensitivity** placeholder output with notes for environment-specific prediction tools.

## Script
Run `pipeline/10_advanced_clinical_omics_analysis.R` after generating:
- `artifacts/models/best_model_artifact.rds`
- `artifacts/reports/merged_subgroup_scored_data.rds`

## Primary outputs
- `artifacts/reports/univariate_cox_signature.csv`
- `artifacts/reports/lasso_cox_coefficients.csv`
- `artifacts/reports/lasso_cox_km.png`
- `artifacts/reports/lasso_cox_timeROC_auc.csv`
- `artifacts/reports/nomogram.png` (optional)
- `artifacts/reports/immune_infiltration_scores.csv` (optional)
- `artifacts/reports/tmb_riskgroup_wilcox.csv` (if `TMB` exists)
- `artifacts/reports/drug_sensitivity_note.txt`

## Required columns
- Survival: `OS_time`, `OS_event`
- Label: `label` with Normal/Cancer levels
- Optional: `TMB`, `age`, `stage`, `subtype`

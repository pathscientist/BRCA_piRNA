# Rscript Protocol (One-command execution guide)

This protocol explains exactly how to run the pipeline with `Rscript` and where to check outputs.

## 1) Prerequisites

- Install R (version 4.2+ recommended).
- Ensure `Rscript` is available in your shell.
- Set your working directory to the repository root (`PiRNAs`).

Quick check:

```bash
Rscript --version
```

If that command fails, install R first.

## 2) Prepare required input files

Place your processed files in `processed_results/` (or update paths in `pipeline/00_config.R`):

- `processed_results/BRCA1_processed.csv`
- `processed_results/PRJNA294226_processed.csv`
- `processed_results/PRJNA482141_processed.csv`
- `processed_results/PRJNA808405_processed.csv`
- `processed_results/PRJNA934049_processed.csv`
- `processed_results/yyfbatch1_processed.csv`
- `processed_results/yyfbatch2_processed.csv`

Required columns per file:

- `sample_id`
- `label` (`Normal` / `Cancer`)
- piRNA expression columns

Recommended for downstream modules:

- subgroup: `age`, `stage`, `subtype`
- survival: `OS_time`, `OS_event`

## 3) Configure pipeline

Edit `pipeline/00_config.R`:

- `config$input$dataset_files` -> real file paths
- `config$input$validation_batch` -> `yyfbatch1` or `yyfbatch2`
- `config$input$batch_col`, `cohort_col`, `label_col`, `sample_id_col`
- `config$model$max_features` (default already `<= 10`)

## 4) Run preflight (must run first)

```bash
Rscript pipeline/00b_data_preflight.R
```

Check reports:

- `artifacts/reports/data_preflight_report.csv`
- `artifacts/reports/data_preflight_report.md`

If preflight fails, fix missing files/columns before continuing.

## 5) Standard full run (manual step-by-step)

```bash
Rscript pipeline/01_load_qc.R
Rscript pipeline/02_feature_selection.R
Rscript pipeline/03_model_training.R
Rscript pipeline/04_model_selection.R
Rscript pipeline/05_external_validation.R
Rscript pipeline/06_reports_plots.R
```

Optional downstream:

```bash
Rscript pipeline/07_clinical_subgroup_analysis.R
Rscript pipeline/08_downstream_biology_survival_meta.R
Rscript pipeline/09_gene_network_analysis.R
Rscript pipeline/10_advanced_clinical_omics_analysis.R
```

## 6) One-command best independent cohort selection

This runs both `yyfbatch1` and `yyfbatch2` as independent cohorts and chooses the better one.

```bash
Rscript pipeline/11_select_best_independent_cohort.R
```

Outputs:

- `artifacts/batch_compare_summary.csv`
- `artifacts/batch_compare_best.csv`
- `artifacts/batch_compare_summary.md`

Hard guards in script 11:

- best external AUC must be `>= 0.80`
- selected feature count must be `<= 10`

## 7) Troubleshooting

- **`Rscript: command not found`**: install R and reopen terminal.
- **Package missing error**: install required package in R, then re-run preflight.
- **Column missing error**: correct column names in source CSVs or `00_config.R` mapping.
- **AUC guard failure (`<0.80`)**: review feature selection thresholds, batch strategy, and cohort quality; then re-run script 11.

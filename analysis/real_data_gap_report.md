# Real-data gap report (what is currently missing)

This report answers: **what real data are still missing to produce real results with the current pipeline**.

## Repository status check
A repository scan shows scripts/docs only; no `processed_results/` dataset files were found in the repo snapshot.

## Missing files to run end-to-end with real data

You need these real input files available at runtime (or update paths in `config$input$dataset_files`):

1. `processed_results/BRCA1_processed.csv`
2. `processed_results/PRJNA294226_processed.csv`
3. `processed_results/PRJNA482141_processed.csv`
4. `processed_results/PRJNA808405_processed.csv`
5. `processed_results/PRJNA934049_processed.csv`
6. `processed_results/yyfbatch1_processed.csv`
7. `processed_results/yyfbatch2_processed.csv`

## Required columns inside each processed file

Minimum for base modeling:
- `sample_id`
- `label` (`Normal`/`Cancer`)
- piRNA expression columns

Recommended to avoid implicit auto-filling and to support stratified analyses:
- `batch`
- `cohort`
- `age`
- `stage`
- `subtype`
- `OS_time`
- `OS_event`

Optional (only if needed for specific modules):
- `TMB` (for TMB association in script 10)

## Optional additional file

Only if your processed CSV files do **not** contain `label`:
- `labels_mapping.csv` with columns `sample_id`, `label`
- and set `config$input$labels_path` to this mapping file.

## Extra requirements for advanced modules (real results)

For `pipeline/10_advanced_clinical_omics_analysis.R` real outputs:
- survival columns (`OS_time`, `OS_event`) must be complete enough for Cox/KM
- optional package-dependent components require installation (`timeROC`, `rms`, `survminer`, `immunedeconv`)

For `pipeline/09_gene_network_analysis.R`:
- `artifacts/reports/pearson_correlated_genes.csv` must exist from upstream correlation step.

## Generated artifacts that must exist before advanced scripts

- `artifacts/models/best_model_artifact.rds`
- `artifacts/reports/merged_subgroup_scored_data.rds`
- `artifacts/reports/pearson_correlated_genes.csv` (for network script)


# Quick check of uploaded legacy script (`pirna_rf_model (4).r`)

## What I checked
- Verified the script exists and reviewed top sections for data assumptions, model flow, and dependency usage.
- Confirmed it assumes an in-memory object `pirna_data` and does not include strict file/column preflight.

## Key compatibility notes with current pipeline
1. The script expects `pirna_data` to be preloaded and `label` to exist.
2. It uses positional feature selection (`[, -1]`) which is fragile if `label` is not column 1.
3. It uses `gather()` but does not load `tidyr` explicitly.
4. It does not have file-level checks for missing required columns.

## Action added now
- Added `pipeline/00b_data_preflight.R` to run a one-command file/column validation before pipeline execution.
- This checker writes:
  - `artifacts/reports/data_preflight_report.csv`
  - `artifacts/reports/data_preflight_report.md`

## One-command usage
```bash
Rscript pipeline/00b_data_preflight.R
```

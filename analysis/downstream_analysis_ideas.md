# Downstream analysis ideas after selecting key piRNA features

This note maps your target figure/caption ideas to concrete analyses and output files.

## 1) Function prediction via Pearson + Reactome

- Compute Pearson correlations between each selected piRNA and mRNA expression across matched samples.
- Keep gene pairs by `|r| >= threshold` and FDR cutoff.
- Run Reactome enrichment on correlated genes.
- Draw a tri-partite network (`piRNA -> gene -> pathway`) with pathway-colored modules.

Suggested outputs:
- `pearson_correlated_genes.csv`
- `reactome_enrichment.csv`
- `pirna_gene_pathway_network.png`

## 2) Stratified diagnostic robustness (ROC + boxplots)

For age/sex/histology/smoking/AJCC stage (or available strata):
- Draw ROC curves by subgroup with AUC + 95% CI.
- Draw T-score tumor-vs-control boxplots within each subgroup.
- Use Mann–Whitney U test and annotate exact p values.
- Print subgroup sample sizes under each panel.

Suggested outputs:
- `strat_<factor>_auc_ci.csv`
- `strat_<factor>_roc.png`
- `strat_<factor>_mann_whitney.csv`
- `strat_<factor>_tscore_box.png`

## 3) Prognostic analysis (KM + Cox)

- For each piRNA biomarker, define high/low by median split (or pre-registered cutoff), then run Kaplan–Meier + log-rank.
- For composite risk score, convert to binary risk groups (high/low) using training-derived threshold.
- Plot KM curves with at-risk tables.
- Run multivariable Cox with binary score + clinical covariates; report HR/95%CI/p-values.

Suggested outputs:
- `km_<feature>.png`
- `km_composite_score.png`
- `multivariable_cox_results.csv`
- `cox_binary_score_forest.png`

## 4) Cross-dataset consistency and meta-analysis

- Draw heatmap of signature expression/T-score across datasets.
- For each key piRNA/tRF, run per-dataset SMD (tumor-control) and random-effects meta-analysis forest plot.
- Build Spearman correlation heatmap between signature RNAs and derived fragments/targets.

Suggested outputs:
- `signature_expression_heatmap.png`
- `meta_forest_<feature>.png`
- `signature_spearman_heatmap.png`

## 5) Functional biology enrichment (KEGG/GO)

- Build functional interaction network between signatures and correlated mRNAs.
- Perform KEGG enrichment on target/correlated gene lists.
- Perform GO BP/CC/MF enrichment and visualize as bubble plots.

Suggested outputs:
- `kegg_enrichment_bar.png`
- `go_bp_bubble.png`
- `go_cc_bubble.png`
- `go_mf_bubble.png`
- `functional_interaction_network.png`

## 6) Suggested report structure (for manuscript-ready figures)

1. Feature discovery and compact signature (<10 markers)
2. Internal/external diagnostic performance
3. Stratified ROC and subgroup T-score comparisons
4. Survival impact (single markers + composite score)
5. Multivariable Cox independence
6. Functional interpretation (Pearson/Reactome/KEGG/GO/network)
7. Cross-dataset consistency + meta-analysis

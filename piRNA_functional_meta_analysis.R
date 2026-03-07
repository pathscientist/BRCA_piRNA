################################################################################
#                                                                              #
#   Breast Cancer piRNA — Meta-Analysis & Functional Prediction                #
#                                                                              #
#   This script runs AFTER piRNA_multicohort_pipeline.R and adds:              #
#                                                                              #
#   PART A — Meta-Analysis Across Datasets                                     #
#     1. Expression heatmap of signature piRNAs across all datasets            #
#        (T-score color scale: orange→purple)                                  #
#     2. SMD forest plot per piRNA across datasets                             #
#        (individual effect + overall diamond from random-effects meta)        #
#     3. Spearman correlation heatmap between signature piRNAs                 #
#                                                                              #
#   PART B — Functional Prediction                                             #
#     4. Pearson's correlation between signature piRNAs and all genes          #
#        (requires mRNA expression data from the same cohort)                  #
#     5. KEGG pathway enrichment analysis + bar plot                           #
#     6. GO enrichment analysis (BP, CC, MF) + bubble plots                   #
#     7. Reactome pathway enrichment                                           #
#     8. Functional interaction network (piRNA → genes → pathways)             #
#                                                                              #
#   Requires objects from main pipeline:                                       #
#     - combat_df_all, model, top_feats, ready_list / clean_list               #
#   OR loads from saved results if run independently                           #
#                                                                              #
################################################################################

start_time_fm <- Sys.time()

# ==============================================================================
# 0. PACKAGES
# ==============================================================================

# CRAN packages
cran_pkgs <- c(
  "ggplot2", "dplyr", "tidyr", "pheatmap", "RColorBrewer",
  "meta", "igraph", "ggraph", "ggrepel", "scales",
  "gridExtra", "cowplot", "reshape2", "corrplot", "ggExtra"
)

# Bioconductor packages
bioc_pkgs <- c(
  "clusterProfiler", "org.Hs.eg.db", "DOSE",
  "enrichplot", "ReactomePA", "AnnotationDbi"
)

# Install CRAN packages
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE, quiet = TRUE)
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", quiet = TRUE)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(RColorBrewer)
  library(meta)
  library(igraph)
  library(ggraph)
  library(ggrepel)
  library(scales)
  library(gridExtra)
  library(cowplot)
  library(reshape2)
  library(corrplot)
  library(ggExtra)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)
  library(ReactomePA)
  library(AnnotationDbi)
})

cat("All meta-analysis & functional analysis packages loaded.\n")

# ==============================================================================
# 0.1 SETTINGS
# ==============================================================================
SEED <- 2024
set.seed(SEED)

pub_theme <- theme_bw() +
  theme(
    panel.grid.major  = element_line(color = "grey92"),
    panel.grid.minor  = element_blank(),
    axis.text         = element_text(size = 11, color = "black"),
    axis.title        = element_text(size = 13, face = "bold"),
    plot.title        = element_text(size = 14, hjust = 0.5, face = "bold"),
    plot.subtitle     = element_text(size = 10, hjust = 0.5, color = "grey40"),
    panel.border      = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.text       = element_text(size = 10),
    strip.text        = element_text(size = 11, face = "bold")
  )

dir.create("results/meta_analysis",   recursive = TRUE, showWarnings = FALSE)
dir.create("results/functional",      recursive = TRUE, showWarnings = FALSE)
dir.create("results/network",         recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 0.2 LOAD PIPELINE RESULTS
# ==============================================================================
cat("\n========== LOADING PIPELINE RESULTS ==========\n")

if (!exists("model")) {
  if (file.exists("results/models/final_model.rds")) {
    model <- readRDS("results/models/final_model.rds")
    cat("  Loaded model.\n")
  } else {
    stop("Run piRNA_multicohort_pipeline.R first.")
  }
}

if (!exists("top_feats")) {
  if (file.exists("results/models/final_features.rds")) {
    top_feats <- readRDS("results/models/final_features.rds")
  } else if (file.exists("results/feature_selection/final_features.txt")) {
    top_feats <- readLines("results/feature_selection/final_features.txt")
    top_feats <- top_feats[nchar(top_feats) > 0]
  }
}

if (!exists("combat_df_all")) {
  stop("combat_df_all not found. Please run the main pipeline first or load it.")
}

cat("piRNA signature:", paste(top_feats, collapse = ", "), "\n")
cat("Total samples:", nrow(combat_df_all), "\n")

# Extract dataset-level info
gene_cols <- setdiff(colnames(combat_df_all),
                     c("Group", "Batch", "T_Score", "T_Score_binary", "T_Score_all",
                       "Age", "Stage", "Subtype", "Age_binary", "Stage_binary",
                       "Age_Group", "Outcome01", "OS_time", "OS_status"))


# ══════════════════════════════════════════════════════════════════════════════
#                    PART A: META-ANALYSIS ACROSS DATASETS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat(paste(rep("═", 70), collapse = ""), "\n")
cat("  PART A: META-ANALYSIS ACROSS DATASETS\n")
cat(paste(rep("═", 70), collapse = ""), "\n")


# ==============================================================================
# A1. EXPRESSION HEATMAP OF SIGNATURE piRNAs ACROSS DATASETS
#     Color scale: T-score values 1.0 (orange) → -1.0 (purple)
# ==============================================================================
cat("\n========== A1: Expression Heatmap Across Datasets ==========\n")

# Compute mean T-score (standardized expression) per piRNA per dataset
datasets <- unique(combat_df_all$Batch)
heatmap_data <- matrix(NA, nrow = length(top_feats), ncol = length(datasets),
                       dimnames = list(top_feats, datasets))

for (ds in datasets) {
  ds_data <- combat_df_all[combat_df_all$Batch == ds, ]
  for (feat in top_feats) {
    if (feat %in% colnames(ds_data)) {
      tumor_vals  <- ds_data[[feat]][ds_data$Group == "Tumor"]
      normal_vals <- ds_data[[feat]][ds_data$Group == "Normal"]
      # Compute standardized mean difference (Cohen's d approximation)
      if (length(tumor_vals) >= 2 && length(normal_vals) >= 2) {
        pooled_sd <- sqrt((var(tumor_vals) * (length(tumor_vals) - 1) +
                           var(normal_vals) * (length(normal_vals) - 1)) /
                          (length(tumor_vals) + length(normal_vals) - 2))
        if (pooled_sd > 0) {
          heatmap_data[feat, ds] <- (mean(tumor_vals) - mean(normal_vals)) / pooled_sd
        }
      }
    }
  }
}

# Clip to [-1, 1] for color scale
heatmap_data_clipped <- pmin(pmax(heatmap_data, -1), 1)

# Annotation: sample sizes per dataset
ds_info <- data.frame(
  Dataset = datasets,
  N_Tumor  = sapply(datasets, function(d)
    sum(combat_df_all$Group[combat_df_all$Batch == d] == "Tumor")),
  N_Normal = sapply(datasets, function(d)
    sum(combat_df_all$Group[combat_df_all$Batch == d] == "Normal"))
)
rownames(ds_info) <- datasets

annotation_col <- data.frame(
  `N(Tumor)`  = ds_info$N_Tumor,
  `N(Normal)` = ds_info$N_Normal,
  row.names = datasets,
  check.names = FALSE
)

# Color palette: purple → white → orange
heatmap_colors <- colorRampPalette(c("#7B3294", "#F7F7F7", "#E66101"))(100)

png("results/meta_analysis/expression_heatmap_datasets.png",
    width = 10, height = 6, units = "in", res = 300)
pheatmap(
  heatmap_data_clipped,
  color = heatmap_colors,
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  main = "Expression Levels of Signature piRNAs Across Datasets",
  fontsize = 10,
  fontsize_row = 11,
  fontsize_col = 10,
  border_color = "grey80",
  cellwidth = 45,
  cellheight = 25,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  legend_breaks = c(-1, -0.5, 0, 0.5, 1),
  legend_labels = c("-1.0", "-0.5", "0", "0.5", "1.0")
)
dev.off()

cat("  Expression heatmap saved.\n")
cat("  Dataset summary:\n")
print(ds_info, row.names = FALSE)


# ==============================================================================
# A2. SMD FOREST PLOT (Random-Effects Meta-Analysis per piRNA)
#     Individual dataset = gray box (sized by N), blue CI lines
#     Overall = red diamond
# ==============================================================================
cat("\n========== A2: SMD Forest Plot (Meta-Analysis) ==========\n")

# Compute SMD + SE for each piRNA × dataset
compute_smd <- function(tumor_vals, normal_vals) {
  n1 <- length(tumor_vals)
  n2 <- length(normal_vals)
  if (n1 < 3 || n2 < 3) return(list(smd = NA, se = NA, n1 = n1, n2 = n2))

  m1 <- mean(tumor_vals); m2 <- mean(normal_vals)
  s1 <- sd(tumor_vals);   s2 <- sd(normal_vals)
  pooled_sd <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))

  if (pooled_sd == 0) return(list(smd = 0, se = NA, n1 = n1, n2 = n2))

  # Hedges' g (bias-corrected SMD)
  d <- (m1 - m2) / pooled_sd
  J <- 1 - 3 / (4 * (n1 + n2 - 2) - 1)  # correction factor
  g <- d * J
  se_g <- sqrt((n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2)))

  list(smd = g, se = se_g, n1 = n1, n2 = n2)
}

# Run meta-analysis for each piRNA
meta_results_all <- list()

for (feat in top_feats) {
  smd_list <- list()
  for (ds in datasets) {
    ds_data <- combat_df_all[combat_df_all$Batch == ds, ]
    if (!feat %in% colnames(ds_data)) next

    tumor_vals  <- ds_data[[feat]][ds_data$Group == "Tumor"]
    normal_vals <- ds_data[[feat]][ds_data$Group == "Normal"]
    res <- compute_smd(tumor_vals, normal_vals)

    if (!is.na(res$smd) && !is.na(res$se)) {
      smd_list[[ds]] <- data.frame(
        Dataset = ds, SMD = res$smd, SE = res$se,
        N_Tumor = res$n1, N_Normal = res$n2,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(smd_list) < 2) next
  smd_df <- do.call(rbind, smd_list)

  # Random-effects meta-analysis
  meta_obj <- tryCatch({
    metagen(TE = smd_df$SMD, seTE = smd_df$SE,
            studlab = smd_df$Dataset,
            sm = "SMD", method.tau = "REML",
            common = FALSE, random = TRUE)
  }, error = function(e) NULL)

  if (!is.null(meta_obj)) {
    meta_results_all[[feat]] <- list(smd_df = smd_df, meta = meta_obj)
    cat(sprintf("  %s: Overall SMD=%.3f (%.3f, %.3f), p=%.3e, I2=%.1f%%\n",
                feat, meta_obj$TE.random,
                meta_obj$lower.random, meta_obj$upper.random,
                meta_obj$pval.random, meta_obj$I2 * 100))
  }
}

# --- Custom ggplot2-based SMD Forest Plot ---
cat("\nGenerating SMD forest plots...\n")

for (feat in names(meta_results_all)) {
  res <- meta_results_all[[feat]]
  smd_df <- res$smd_df
  meta_obj <- res$meta

  # Build data for plot
  plot_df <- smd_df %>%
    mutate(
      Lower = SMD - 1.96 * SE,
      Upper = SMD + 1.96 * SE,
      Weight = 1 / SE^2,
      Type = "Study"
    )

  # Add overall row
  overall_row <- data.frame(
    Dataset  = "Overall",
    SMD      = meta_obj$TE.random,
    SE       = meta_obj$seTE.random,
    N_Tumor  = sum(smd_df$N_Tumor),
    N_Normal = sum(smd_df$N_Normal),
    Lower    = meta_obj$lower.random,
    Upper    = meta_obj$upper.random,
    Weight   = NA,
    Type     = "Summary",
    stringsAsFactors = FALSE
  )
  plot_df <- rbind(plot_df, overall_row)

  # Normalize weights for box sizing
  study_weights <- plot_df$Weight[plot_df$Type == "Study"]
  plot_df$RelWeight <- NA
  plot_df$RelWeight[plot_df$Type == "Study"] <-
    (study_weights / max(study_weights, na.rm = TRUE)) * 4 + 1

  # Y position
  plot_df$y_pos <- rev(seq_len(nrow(plot_df)))

  # Build annotation text
  plot_df$CI_text <- sprintf("%.2f [%.2f, %.2f]", plot_df$SMD, plot_df$Lower, plot_df$Upper)
  plot_df$N_text <- paste0(plot_df$N_Tumor, "/", plot_df$N_Normal)

  p_meta <- ggplot(plot_df, aes(x = SMD, y = y_pos)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    # Study-level: gray squares + blue CI lines
    geom_errorbarh(
      data = filter(plot_df, Type == "Study"),
      aes(xmin = Lower, xmax = Upper),
      height = 0.2, color = "#377EB8", linewidth = 0.7
    ) +
    geom_point(
      data = filter(plot_df, Type == "Study"),
      aes(size = RelWeight),
      shape = 15, color = "grey50"
    ) +
    # Summary: red diamond
    geom_point(
      data = filter(plot_df, Type == "Summary"),
      shape = 23, size = 5, fill = "#E41A1C", color = "#E41A1C"
    ) +
    geom_errorbarh(
      data = filter(plot_df, Type == "Summary"),
      aes(xmin = Lower, xmax = Upper),
      height = 0.3, color = "#E41A1C", linewidth = 1
    ) +
    # Annotation: CI text
    geom_text(aes(x = max(plot_df$Upper, na.rm = TRUE) + 0.5,
                  label = CI_text),
              hjust = 0, size = 3) +
    scale_y_continuous(
      breaks = plot_df$y_pos,
      labels = plot_df$Dataset,
      expand = expansion(mult = c(0.08, 0.08))
    ) +
    scale_size_identity() +
    labs(
      title = paste0("Meta-Analysis: ", feat),
      subtitle = sprintf("Random-effects model (I\u00B2 = %.1f%%, p = %.2e)",
                         meta_obj$I2 * 100, meta_obj$pval.random),
      x = "Standardized Mean Difference (Hedges' g)",
      y = ""
    ) +
    pub_theme +
    theme(
      plot.margin = ggplot2::margin(10, 100, 10, 10, unit = "pt"),
      axis.text.y = element_text(size = 10)
    )

  fname <- gsub("[^a-zA-Z0-9_.-]", "_", feat)
  ggsave(paste0("results/meta_analysis/SMD_forest_", fname, ".png"),
         p_meta, width = 11, height = max(4, nrow(plot_df) * 0.6 + 2), dpi = 300)
}

cat("  SMD forest plots saved.\n")

# Save meta-analysis summary table
meta_summary <- do.call(rbind, lapply(names(meta_results_all), function(feat) {
  m <- meta_results_all[[feat]]$meta
  data.frame(
    piRNA       = feat,
    SMD_overall = m$TE.random,
    Lower_95CI  = m$lower.random,
    Upper_95CI  = m$upper.random,
    P_value     = m$pval.random,
    I2          = m$I2,
    Tau2        = m$tau2,
    N_datasets  = m$k,
    stringsAsFactors = FALSE
  )
}))
write.csv(meta_summary, "results/meta_analysis/meta_analysis_summary.csv", row.names = FALSE)
cat("  Meta-analysis summary table saved.\n")


# ==============================================================================
# A3. SPEARMAN CORRELATION HEATMAP BETWEEN SIGNATURE piRNAs
#     Color scale: -1 (red) → 1 (dark blue)
# ==============================================================================
cat("\n========== A3: Spearman Correlation Heatmap ==========\n")

# Compute Spearman correlation among signature piRNAs across all samples
sig_expr <- combat_df_all[, top_feats, drop = FALSE]
cor_matrix <- cor(sig_expr, method = "spearman", use = "pairwise.complete.obs")

# Also compute per-dataset correlation and average
cor_per_dataset <- list()
for (ds in datasets) {
  ds_data <- combat_df_all[combat_df_all$Batch == ds, top_feats, drop = FALSE]
  if (nrow(ds_data) >= 10) {
    cor_per_dataset[[ds]] <- cor(ds_data, method = "spearman",
                                 use = "pairwise.complete.obs")
  }
}

# Average correlation across datasets
if (length(cor_per_dataset) > 1) {
  avg_cor <- Reduce("+", cor_per_dataset) / length(cor_per_dataset)
  # Replace NaN/NA with 0 to avoid hclust failure
  avg_cor[is.na(avg_cor)] <- 0
} else {
  avg_cor <- cor_matrix
}

# Color scale: red (-1) → white (0) → dark blue (1)
cor_colors <- colorRampPalette(c("#B2182B", "#F4A582", "#FDDBC7", "#F7F7F7",
                                  "#D1E5F0", "#92C5DE", "#2166AC"))(100)

# Plot 1: Overall correlation
png("results/meta_analysis/spearman_correlation_all.png",
    width = 8, height = 7, units = "in", res = 300)
pheatmap(
  cor_matrix,
  color = cor_colors,
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Spearman Correlation: Signature piRNAs (All Samples)",
  fontsize = 11,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  cellwidth = 50,
  cellheight = 50,
  border_color = "grey80"
)
dev.off()

# Plot 2: Averaged across datasets
png("results/meta_analysis/spearman_correlation_averaged.png",
    width = 8, height = 7, units = "in", res = 300)
pheatmap(
  avg_cor,
  color = cor_colors,
  breaks = seq(-1, 1, length.out = 101),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Spearman Correlation: Aggregated Across All Datasets",
  fontsize = 11,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  cellwidth = 50,
  cellheight = 50,
  border_color = "grey80"
)
dev.off()

cat("  Spearman correlation heatmaps saved.\n")


# ══════════════════════════════════════════════════════════════════════════════
#                    PART B: FUNCTIONAL PREDICTION
#     Predict piRNA functions via piRNA–mRNA correlation & pathway enrichment
# ══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat(paste(rep("\u2550", 70), collapse = ""), "\n")
cat("  PART B: FUNCTIONAL PREDICTION (piRNA\u2013mRNA Correlation)\n")
cat(paste(rep("\u2550", 70), collapse = ""), "\n")


# ==============================================================================
# B1. LOAD mRNA EXPRESSION DATA (TCGA-BRCA)
# ==============================================================================
cat("\n========== B1: Loading mRNA Expression Data ==========\n")

mRNA_df <- NULL
mRNA_cache <- "processed_results/TCGA_BRCA_mRNA_log2tpm.rds"

# --- Option 1: Load from cached RDS ---
if (file.exists(mRNA_cache)) {
  cat("  Loading cached mRNA data:", mRNA_cache, "\n")
  mRNA_df <- readRDS(mRNA_cache)
}

# --- Option 2: Load from CSV/TXT ---
if (is.null(mRNA_df)) {
  mRNA_paths <- c(
    "mRNA_expression/TCGA_BRCA_mRNA.csv",
    "data/mRNA_expression.csv", "data/TCGA_BRCA_mRNA.csv",
    "data/mRNA_expression.txt", "data/TCGA_BRCA_mRNA_tpm.csv"
  )
  for (fp in mRNA_paths) {
    if (file.exists(fp)) {
      cat("  Loading mRNA data from:", fp, "\n")
      mRNA_df <- tryCatch({
        if (grepl("\\.csv$", fp))
          read.csv(fp, row.names = 1, check.names = FALSE)
        else
          read.table(fp, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
      }, error = function(e) NULL)
      if (!is.null(mRNA_df)) break
    }
  }
}

# --- Option 3: Download TCGA-BRCA mRNA via TCGAbiolinks ---
if (is.null(mRNA_df)) {
  cat("  No local mRNA file found. Downloading TCGA-BRCA via TCGAbiolinks...\n")

  for (pkg in c("TCGAbiolinks", "SummarizedExperiment")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      BiocManager::install(pkg, quiet = TRUE, update = FALSE)
  }
  library(TCGAbiolinks)
  library(SummarizedExperiment)

  mRNA_df <- tryCatch({
    query <- GDCquery(
      project = "TCGA-BRCA",
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts"
    )
    GDCdownload(query, method = "api", directory = "GDCdata")
    se <- GDCprepare(query, directory = "GDCdata")

    # Extract TPM, map to gene symbols, log2-transform
    tpm_mat <- assay(se, "tpm_unstrand")
    gene_sym <- rowData(se)$gene_name
    keep <- !is.na(gene_sym) & gene_sym != "" & !duplicated(gene_sym)
    tpm_mat <- tpm_mat[keep, ]
    rownames(tpm_mat) <- gene_sym[keep]
    tpm_log <- log2(tpm_mat + 1)

    # Transpose: samples in rows, genes in columns
    df_out <- as.data.frame(t(tpm_log))
    dir.create(dirname(mRNA_cache), showWarnings = FALSE, recursive = TRUE)
    saveRDS(df_out, mRNA_cache)
    cat("  TCGA-BRCA mRNA data downloaded and cached.\n")
    df_out
  }, error = function(e) {
    cat("  TCGAbiolinks download failed:", conditionMessage(e), "\n")
    NULL
  })
}

# --- Validate ---
run_partB <- !is.null(mRNA_df) && ncol(mRNA_df) >= 100
if (!run_partB) {
  cat("\n  *** mRNA expression data not available. ***\n")
  cat("  To enable functional prediction, provide mRNA expression data:\n")
  cat("    - Place a CSV at data/mRNA_expression.csv\n")
  cat("    - Rows = samples (TCGA barcodes), Columns = gene symbols\n")
  cat("  Or install TCGAbiolinks for automatic download.\n")
} else {
  cat("  mRNA data:", nrow(mRNA_df), "samples x", ncol(mRNA_df), "genes\n")
}


# ==============================================================================
# B2. piRNA–mRNA PEARSON CORRELATION
# ==============================================================================
if (run_partB) {

cat("\n========== B2: piRNA\u2013mRNA Pearson Correlation ==========\n")

# Match samples between piRNA and mRNA data using TCGA barcodes
# TCGA barcodes: first 15 characters identify the sample
pirna_samples <- rownames(combat_df_all[combat_df_all$Group == "Tumor" &
                                         combat_df_all$Batch == "BRCA1", ])
mrna_samples  <- rownames(mRNA_df)

# Try matching by truncated TCGA barcodes (first 15 chars)
pirna_short <- substr(pirna_samples, 1, 15)
mrna_short  <- substr(mrna_samples, 1, 15)

common_short <- intersect(pirna_short, mrna_short)

if (length(common_short) < 20) {
  # Try exact rowname match
  common_exact <- intersect(pirna_samples, mrna_samples)
  if (length(common_exact) >= 20) {
    pirna_matched <- common_exact
    mrna_matched  <- common_exact
    cat("  Matched by exact rownames:", length(common_exact), "samples\n")
  } else {
    cat("  WARNING: Only", max(length(common_short), length(common_exact)),
        "samples matched. Using all tumor samples with index matching.\n")
    # Fallback: use all tumor samples by row order if names don't match
    tumor_idx <- combat_df_all$Group == "Tumor"
    n_use <- min(sum(tumor_idx), nrow(mRNA_df))
    pirna_matched <- rownames(combat_df_all[tumor_idx, ])[1:n_use]
    mrna_matched  <- rownames(mRNA_df)[1:n_use]
  }
} else {
  # Build matched index
  pirna_idx <- match(common_short, pirna_short)
  mrna_idx  <- match(common_short, mrna_short)
  pirna_matched <- pirna_samples[pirna_idx]
  mrna_matched  <- mrna_samples[mrna_idx]
  cat("  Matched by TCGA barcodes:", length(common_short), "samples\n")
}

# Extract matched expression matrices
expr_pirna <- as.matrix(combat_df_all[pirna_matched, top_feats, drop = FALSE])
expr_mrna  <- as.matrix(mRNA_df[mrna_matched, , drop = FALSE])

# Remove mRNA genes with zero variance
mrna_var <- apply(expr_mrna, 2, var, na.rm = TRUE)
expr_mrna <- expr_mrna[, mrna_var > 0, drop = FALSE]
cat("  Computing Pearson correlations:", length(top_feats), "piRNAs x",
    ncol(expr_mrna), "mRNAs...\n")

# Compute correlations for each signature piRNA vs all mRNAs
cor_results <- list()
for (feat in top_feats) {
  pirna_vec <- expr_pirna[, feat]
  cors <- apply(expr_mrna, 2, function(mrna_vec) {
    complete <- complete.cases(pirna_vec, mrna_vec)
    if (sum(complete) < 20) return(c(r = NA, p = NA))
    ct <- cor.test(pirna_vec[complete], mrna_vec[complete], method = "pearson")
    c(r = ct$estimate, p = ct$p.value)
  })
  cor_df <- data.frame(
    piRNA     = feat,
    Gene      = colnames(cors),
    Pearson_r = cors["r.cor", ],
    P_value   = cors["p", ],
    stringsAsFactors = FALSE
  )
  cor_df <- cor_df[!is.na(cor_df$P_value), ]
  cor_df$FDR <- p.adjust(cor_df$P_value, method = "BH")
  cor_df$Significant <- abs(cor_df$Pearson_r) > 0.3 & cor_df$FDR < 0.05
  cor_results[[feat]] <- cor_df
}

all_cors <- do.call(rbind, cor_results)
rownames(all_cors) <- NULL
sig_cors <- all_cors[all_cors$Significant, ]

cat(sprintf("  Total significant piRNA\u2013mRNA pairs: %d (|r|>0.3, FDR<0.05)\n",
            nrow(sig_cors)))
for (feat in top_feats) {
  n_sig <- sum(sig_cors$piRNA == feat)
  cat(sprintf("    %s: %d correlated mRNAs\n", feat, n_sig))
}

# Save tables
dir.create("results/functional", showWarnings = FALSE, recursive = TRUE)
write.csv(all_cors, "results/functional/pearson_piRNA_mRNA_all.csv", row.names = FALSE)
write.csv(sig_cors, "results/functional/pearson_piRNA_mRNA_significant.csv", row.names = FALSE)


# ==============================================================================
# B2a. DOT PLOT: piRNA–mRNA Correlation (Figure A style)
#      Y-axis = Gene, X-axis = piRNA, size = -log10(FDR), color = Pearson r
# ==============================================================================
cat("\n  Generating piRNA\u2013mRNA correlation dot plot...\n")

# Select top correlated genes per piRNA for the plot
top_n_per_pirna <- 20
dot_df <- sig_cors %>%
  group_by(piRNA) %>%
  arrange(desc(abs(Pearson_r))) %>%
  slice_head(n = top_n_per_pirna) %>%
  ungroup()

if (nrow(dot_df) > 0) {
  # Cap -log10(FDR) for display
  dot_df$neg_log10_fdr <- pmin(-log10(dot_df$FDR), 20)

  # Order genes by average |r| across piRNAs
  gene_order <- dot_df %>%
    group_by(Gene) %>%
    summarise(mean_r = mean(abs(Pearson_r)), .groups = "drop") %>%
    arrange(desc(mean_r)) %>%
    pull(Gene)
  dot_df$Gene <- factor(dot_df$Gene, levels = rev(gene_order))

  p_dotplot <- ggplot(dot_df, aes(x = piRNA, y = Gene)) +
    geom_point(aes(size = neg_log10_fdr, color = Pearson_r), alpha = 0.85) +
    scale_size_continuous(
      range = c(1.5, 7), name = expression(-log[10](FDR)),
      breaks = c(5, 10, 15)
    ) +
    scale_color_gradient2(
      low = "#B2182B", mid = "#F7F7F7", high = "#2166AC", midpoint = 0,
      name = "Pearson\ncorrelation",
      limits = c(-1, 1)
    ) +
    labs(
      title = "Genes Correlated with Signature piRNAs",
      subtitle = "Pearson correlation (FDR < 0.05, |r| > 0.3)",
      x = "", y = ""
    ) +
    pub_theme +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      axis.text.y = element_text(size = 8),
      panel.grid.major = element_line(color = "grey92"),
      legend.position = "right"
    )

  n_genes <- length(unique(dot_df$Gene))
  ggsave("results/functional/piRNA_mRNA_correlation_dotplot.png",
         p_dotplot, width = max(6, length(top_feats) * 1.5 + 4),
         height = max(8, n_genes * 0.22 + 3), dpi = 300,
         limitsize = FALSE)
  cat("  Dot plot saved: results/functional/piRNA_mRNA_correlation_dotplot.png\n")
}


# ==============================================================================
# B3. PATHWAY ENRICHMENT (KEGG, GO, Reactome)
#     Use correlated mRNA genes to predict piRNA functions
# ==============================================================================
cat("\n========== B3: Pathway Enrichment of Correlated Genes ==========\n")

correlated_genes <- unique(sig_cors$Gene)
cat("  Unique correlated mRNA genes:", length(correlated_genes), "\n")

# Map gene symbols to Entrez IDs
entrez_map <- tryCatch({
  bitr(correlated_genes, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)
}, error = function(e) {
  cat("  SYMBOL mapping failed, trying ALIAS...\n")
  tryCatch(
    bitr(correlated_genes, fromType = "ALIAS", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db),
    error = function(e2) NULL
  )
})

run_enrichment <- !is.null(entrez_map) && nrow(entrez_map) >= 5

if (run_enrichment) {
  entrez_ids <- unique(entrez_map$ENTREZID)
  cat("  Mapped", length(entrez_ids), "genes to Entrez IDs.\n")

  # --- KEGG Enrichment ---
  cat("\n  Running KEGG enrichment...\n")
  kegg_res <- tryCatch({
    enrichKEGG(gene = entrez_ids, organism = "hsa",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2,
               minGSSize = 5, maxGSSize = 500)
  }, error = function(e) {
    cat("    KEGG failed:", conditionMessage(e), "\n"); NULL
  })

  if (!is.null(kegg_res) && nrow(kegg_res@result[kegg_res@result$p.adjust < 0.05, ]) > 0) {
    kegg_sig <- kegg_res@result[kegg_res@result$p.adjust < 0.05, ]
    write.csv(kegg_sig, "results/functional/kegg_enrichment.csv", row.names = FALSE)

    n_show <- min(20, nrow(kegg_sig))
    kegg_plot_df <- head(kegg_sig, n_show) %>%
      mutate(Count = as.numeric(Count),
             Description = factor(Description, levels = rev(Description)))

    p_kegg <- ggplot(kegg_plot_df, aes(x = Count, y = Description, fill = p.adjust)) +
      geom_col(width = 0.7) +
      scale_fill_gradient(low = "#C62828", high = "#FFCDD2",
                          name = "FDR") +
      labs(title = "KEGG Pathway Enrichment",
           subtitle = paste0("Based on ", length(correlated_genes),
                             " mRNAs correlated with signature piRNAs"),
           x = "Gene Count", y = "") +
      pub_theme + theme(axis.text.y = element_text(size = 9))

    ggsave("results/functional/kegg_barplot.png",
           p_kegg, width = 10, height = max(5, n_show * 0.35 + 2), dpi = 300)
    cat("    KEGG:", nrow(kegg_sig), "significant pathways, plot saved.\n")
  } else {
    cat("    No significant KEGG pathways.\n")
  }

  # --- GO Enrichment (BP, CC, MF) ---
  cat("\n  Running GO enrichment...\n")
  go_categories <- c("BP", "CC", "MF")
  go_names <- c(BP = "Biological Process", CC = "Cellular Component",
                MF = "Molecular Function")

  for (ont in go_categories) {
    go_res <- tryCatch({
      enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db,
               ont = ont, pvalueCutoff = 0.05, qvalueCutoff = 0.2,
               readable = TRUE, minGSSize = 5, maxGSSize = 500)
    }, error = function(e) { NULL })

    if (!is.null(go_res) && nrow(go_res@result[go_res@result$p.adjust < 0.05, ]) > 0) {
      go_sig <- go_res@result[go_res@result$p.adjust < 0.05, ]
      write.csv(go_sig, paste0("results/functional/go_", tolower(ont), "_enrichment.csv"),
                row.names = FALSE)

      n_show <- min(15, nrow(go_sig))
      go_plot_df <- head(go_sig, n_show) %>%
        mutate(
          Count = as.numeric(Count),
          GeneRatio_num = sapply(GeneRatio, function(x) {
            parts <- as.numeric(strsplit(x, "/")[[1]])
            parts[1] / parts[2]
          }),
          Description = factor(Description,
                               levels = rev(Description[order(Count)]))
        )

      p_go <- ggplot(go_plot_df,
                     aes(x = GeneRatio_num, y = Description,
                         size = Count, color = p.adjust)) +
        geom_point(alpha = 0.85) +
        scale_size_continuous(range = c(3, 10), name = "Gene Count") +
        scale_color_gradient(low = "#C62828", high = "#FFCDD2",
                             name = "FDR") +
        labs(title = paste0("GO Enrichment: ", go_names[ont]),
             x = "Gene Ratio", y = "") +
        pub_theme + theme(axis.text.y = element_text(size = 9))

      ggsave(paste0("results/functional/go_", tolower(ont), "_bubble.png"),
             p_go, width = 10, height = max(5, n_show * 0.35 + 2), dpi = 300)
      cat(sprintf("    GO %s: %d significant terms.\n", ont, nrow(go_sig)))
    } else {
      cat(sprintf("    GO %s: no significant terms.\n", ont))
    }
  }

  # --- Reactome Enrichment ---
  cat("\n  Running Reactome enrichment...\n")
  reactome_res <- tryCatch({
    enrichPathway(gene = entrez_ids, organism = "human",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  readable = TRUE, minGSSize = 5, maxGSSize = 500)
  }, error = function(e) {
    cat("    Reactome failed:", conditionMessage(e), "\n"); NULL
  })

  reactome_sig <- NULL
  if (!is.null(reactome_res) &&
      nrow(reactome_res@result[reactome_res@result$p.adjust < 0.05, ]) > 0) {
    reactome_sig <- reactome_res@result[reactome_res@result$p.adjust < 0.05, ]
    write.csv(reactome_sig, "results/functional/reactome_enrichment.csv",
              row.names = FALSE)

    n_show <- min(15, nrow(reactome_sig))
    react_plot_df <- head(reactome_sig, n_show) %>%
      mutate(Count = as.numeric(Count),
             Description = factor(Description, levels = rev(Description)))

    p_react <- ggplot(react_plot_df,
                      aes(x = Count, y = Description, fill = p.adjust)) +
      geom_col(width = 0.7) +
      scale_fill_gradient(low = "#C62828", high = "#FFCDD2", name = "FDR") +
      labs(title = "Reactome Pathway Enrichment",
           subtitle = "Pathways enriched in piRNA-correlated genes",
           x = "Gene Count", y = "") +
      pub_theme + theme(axis.text.y = element_text(size = 9))

    ggsave("results/functional/reactome_barplot.png",
           p_react, width = 11, height = max(5, n_show * 0.35 + 2), dpi = 300)
    cat("    Reactome:", nrow(reactome_sig), "significant pathways.\n")
  } else {
    cat("    No significant Reactome pathways.\n")
  }


  # ============================================================================
  # B4. FUNCTIONAL NETWORK: piRNA → Correlated Genes
  #     Circular layout: piRNAs as large colored hubs in center,
  #     correlated genes on outer ring, edges colored by direction
  # ============================================================================
  cat("\n========== B4: Functional Interaction Network ==========\n")

  dir.create("results/network", showWarnings = FALSE, recursive = TRUE)

  if (nrow(sig_cors) > 0) {
    # Select top correlated genes per piRNA for the network
    top_genes_net <- 20
    net_cors <- sig_cors %>%
      group_by(piRNA) %>%
      arrange(desc(abs(Pearson_r))) %>%
      slice_head(n = top_genes_net) %>%
      ungroup()

    gene_symbols_in_net <- unique(net_cors$Gene)

    # Build edges: piRNA → gene
    edges_df <- net_cors %>%
      transmute(from = piRNA, to = Gene,
                weight = abs(Pearson_r),
                direction = ifelse(Pearson_r > 0, "positive", "negative"))

    # Build node metadata
    all_node_names <- unique(c(edges_df$from, edges_df$to))
    node_df <- data.frame(
      name = all_node_names,
      type = ifelse(all_node_names %in% top_feats, "piRNA", "gene"),
      stringsAsFactors = FALSE
    )

    # Assign piRNA-specific colors
    pirna_palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#FF7F00", "#984EA3")
    pirna_names <- top_feats[top_feats %in% unique(edges_df$from)]
    pirna_color_map <- setNames(pirna_palette[seq_along(pirna_names)], pirna_names)

    # For each gene, find the piRNA with strongest |r| for grouping
    gene_primary <- edges_df %>%
      group_by(to) %>%
      slice_max(weight, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(name = to, primary_pirna = from)

    node_df <- node_df %>%
      left_join(gene_primary, by = "name") %>%
      mutate(node_size = ifelse(type == "piRNA", 14, 5))

    # Build igraph
    g <- graph_from_data_frame(edges_df, directed = FALSE, vertices = node_df)

    # Custom circular layout: piRNAs in center, genes on outer ring
    layout_mat <- matrix(0, nrow = nrow(node_df), ncol = 2)
    rownames(layout_mat) <- node_df$name

    # Place piRNAs in center
    pirna_nodes <- node_df$name[node_df$type == "piRNA"]
    if (length(pirna_nodes) == 1) {
      layout_mat[pirna_nodes, ] <- c(0, 0)
    } else {
      pirna_angles <- seq(0, 2 * pi, length.out = length(pirna_nodes) + 1)
      pirna_angles <- pirna_angles[-(length(pirna_nodes) + 1)]
      pirna_r <- 0.3
      layout_mat[pirna_nodes, 1] <- cos(pirna_angles) * pirna_r
      layout_mat[pirna_nodes, 2] <- sin(pirna_angles) * pirna_r
    }

    # Place genes on outer ring, grouped by primary piRNA
    gene_order_df <- node_df %>%
      filter(type == "gene") %>%
      mutate(pirna_rank = match(primary_pirna, pirna_names)) %>%
      arrange(pirna_rank, name)
    gene_nodes_ordered <- gene_order_df$name
    gene_angles <- seq(0, 2 * pi, length.out = length(gene_nodes_ordered) + 1)
    gene_angles <- gene_angles[-(length(gene_nodes_ordered) + 1)]
    gene_r <- 1.0
    layout_mat[gene_nodes_ordered, 1] <- cos(gene_angles) * gene_r
    layout_mat[gene_nodes_ordered, 2] <- sin(gene_angles) * gene_r

    # Reorder layout to match igraph vertex order
    layout_mat <- layout_mat[V(g)$name, ]

    set.seed(SEED)

    # Edge colors: blue for positive, red for negative correlation
    edge_colors <- ifelse(E(g)$direction == "positive", "#4A90D9", "#D94A4A")

    # Compute gene label angles for radial orientation
    gene_mask <- V(g)$type == "gene"
    gene_x <- layout_mat[gene_mask, 1]
    gene_y <- layout_mat[gene_mask, 2]
    label_angles <- atan2(gene_y, gene_x) * 180 / pi
    # Flip labels on left half so text reads correctly
    label_angles <- ifelse(abs(label_angles) > 90,
                           label_angles + 180, label_angles)

    p_net <- ggraph(g, layout = "manual",
                    x = layout_mat[, 1], y = layout_mat[, 2]) +
      # Edges: colored by correlation direction, alpha by strength
      geom_edge_link(
        aes(edge_alpha = weight, edge_width = weight),
        edge_colour = edge_colors,
        show.legend = FALSE
      ) +
      scale_edge_alpha_continuous(range = c(0.12, 0.5)) +
      scale_edge_width_continuous(range = c(0.2, 1.3)) +
      # Gene nodes: pink circles on outer ring
      geom_node_point(
        data = . %>% filter(type == "gene"),
        color = "grey50", fill = "#FFB6C1", shape = 21,
        size = 5, stroke = 0.4, alpha = 0.9
      ) +
      # Gene labels: radially oriented around the circle
      geom_node_text(
        data = . %>% filter(type == "gene"),
        aes(label = name, angle = label_angles),
        size = 2.8, color = "grey20",
        hjust = ifelse(gene_x >= 0, 0, 1),
        nudge_x = ifelse(gene_x >= 0, 0.06, -0.06),
        nudge_y = 0
      ) +
      # piRNA nodes: large colored circles in center
      geom_node_point(
        data = . %>% filter(type == "piRNA"),
        aes(fill = name), shape = 21, size = 12,
        color = "grey30", stroke = 1.2
      ) +
      scale_fill_manual(values = pirna_color_map, name = "piRNA") +
      # piRNA labels
      geom_node_text(
        data = . %>% filter(type == "piRNA"),
        aes(label = name), size = 3, fontface = "bold",
        color = "white"
      ) +
      coord_fixed() +
      labs(
        title = "Significant Correlations Between Signature piRNAs and Target Genes",
        subtitle = "Blue edges = positive correlation, Red edges = negative correlation (|r| > 0.3, FDR < 0.05)"
      ) +
      theme_void() +
      theme(
        plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
        plot.margin   = ggplot2::margin(10, 40, 10, 40),
        legend.position = "bottom",
        legend.title  = element_text(face = "bold"),
        legend.text   = element_text(size = 9)
      ) +
      guides(fill = guide_legend(override.aes = list(size = 5)))

    ggsave("results/network/functional_network.png",
           p_net, width = 14, height = 14, dpi = 300)
    cat("  Functional network saved: results/network/functional_network.png\n")

    write.csv(as.data.frame(edges_df),
              "results/network/network_edges.csv", row.names = FALSE)
    write.csv(node_df, "results/network/network_nodes.csv", row.names = FALSE)
  } else {
    cat("  No significant correlations for network construction.\n")
  }


  # ============================================================================
  # B4b. SCATTER PLOTS WITH MARGINAL HISTOGRAMS
  #      For each piRNA, show top 3 correlated genes as scatter + marginal hist
  #      Style inspired by: correlation panels with R, p annotation
  # ============================================================================
  cat("\n========== B4b: piRNA-mRNA Correlation Scatter Plots ==========\n")

  dir.create("results/functional/scatter_plots", showWarnings = FALSE, recursive = TRUE)

  n_top_genes <- 3  # top genes per piRNA

  for (feat in top_feats) {
    feat_cors <- sig_cors[sig_cors$piRNA == feat, ]
    if (nrow(feat_cors) == 0) {
      cat(sprintf("  %s: No significant correlations, skipping scatter plots.\n", feat))
      next
    }

    # Top N genes by |r|
    feat_cors <- feat_cors %>% arrange(desc(abs(Pearson_r)))
    top_gene_names <- head(feat_cors$Gene, n_top_genes)
    top_gene_cors  <- head(feat_cors, n_top_genes)

    scatter_list <- list()

    for (j in seq_along(top_gene_names)) {
      gene_name <- top_gene_names[j]
      r_val  <- top_gene_cors$Pearson_r[j]
      p_val  <- top_gene_cors$P_value[j]

      pirna_vec <- expr_pirna[, feat]
      mrna_vec  <- expr_mrna[, gene_name]
      complete  <- complete.cases(pirna_vec, mrna_vec)

      scat_df <- data.frame(
        piRNA_expr = pirna_vec[complete],
        mRNA_expr  = mrna_vec[complete]
      )

      # Format annotation
      p_label <- if (p_val < 2.2e-16) {
        sprintf("italic(R) == %.2f ~~ italic(p) < 2.2e-16", r_val)
      } else {
        sprintf("italic(R) == %.2f ~~ italic(p) == %.2e", r_val, p_val)
      }

      # Base scatter with colored points and regression line
      p_base <- ggplot(scat_df, aes(x = piRNA_expr, y = mRNA_expr)) +
        geom_point(aes(color = piRNA_expr), size = 1.8, alpha = 0.6,
                   show.legend = FALSE) +
        scale_color_gradient2(low = "#E41A1C", mid = "#8B7BB8", high = "#377EB8",
                              midpoint = median(scat_df$piRNA_expr, na.rm = TRUE)) +
        geom_smooth(method = "lm", color = "#2166AC", fill = "#92C5DE",
                    linewidth = 1, alpha = 0.3, se = TRUE) +
        annotate("text", x = -Inf, y = Inf, label = p_label, parse = TRUE,
                 hjust = -0.05, vjust = 1.5, size = 3.8, fontface = "italic") +
        labs(x = feat, y = gene_name) +
        theme_bw() +
        theme(
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "grey92"),
          axis.text  = element_text(size = 10, color = "black"),
          axis.title = element_text(size = 12, face = "bold"),
          plot.margin = ggplot2::margin(5, 5, 5, 5)
        )

      # Add marginal histograms (top = piRNA, right = mRNA)
      p_marginal <- ggExtra::ggMarginal(
        p_base,
        type = "histogram",
        xFill = "#F4A460",   # sandy brown for piRNA (x-axis)
        yFill = "#8FBC8F",   # dark sea green for mRNA (y-axis)
        margins = "both",
        size = 5
      )

      scatter_list[[j]] <- p_marginal

      # Save individual scatter plot
      fname <- gsub("[^a-zA-Z0-9_.-]", "_", feat)
      gname <- gsub("[^a-zA-Z0-9_.-]", "_", gene_name)
      ggsave(
        paste0("results/functional/scatter_plots/", fname, "_vs_", gname, ".png"),
        p_marginal, width = 5.5, height = 5, dpi = 300
      )
    }

    # Combine top scatter plots into a panel for this piRNA
    if (length(scatter_list) >= 2) {
      fname <- gsub("[^a-zA-Z0-9_.-]", "_", feat)

      combined <- cowplot::plot_grid(
        plotlist = scatter_list,
        ncol = min(3, length(scatter_list)),
        labels = LETTERS[seq_along(scatter_list)],
        label_size = 16, label_fontface = "bold"
      )

      title_grob <- cowplot::ggdraw() +
        cowplot::draw_label(
          paste0("Top Correlated Genes: ", feat),
          fontface = "bold", size = 14, x = 0.5, hjust = 0.5
        )

      final_panel <- cowplot::plot_grid(
        title_grob, combined,
        ncol = 1, rel_heights = c(0.06, 1)
      )

      ggsave(
        paste0("results/functional/scatter_plots/panel_", fname, ".png"),
        final_panel,
        width = min(16, 5.5 * length(scatter_list)),
        height = 5.5, dpi = 300
      )
      cat(sprintf("  %s: %d scatter plots saved (panel + individual).\n",
                  feat, length(scatter_list)))
    }
  }

} else {
  cat("  Gene mapping failed. Cannot run enrichment.\n")
  cat("  Ensure correlated genes use standard HGNC symbols.\n")
}


# ==============================================================================
# B5. PER-piRNA CORRELATED GENE TABLES
# ==============================================================================
cat("\n========== B5: Per-piRNA Correlated Gene Tables ==========\n")

for (feat in top_feats) {
  feat_cors <- sig_cors[sig_cors$piRNA == feat, ]
  if (nrow(feat_cors) == 0) {
    cat(sprintf("  %s: No significant correlated mRNAs.\n", feat))
    next
  }
  feat_cors <- feat_cors %>% arrange(desc(abs(Pearson_r)))
  fname <- gsub("[^a-zA-Z0-9_.-]", "_", feat)
  write.csv(feat_cors,
            paste0("results/functional/correlated_mRNAs_", fname, ".csv"),
            row.names = FALSE)
  cat(sprintf("  %s: %d correlated mRNAs (top: %s, r=%.3f)\n",
              feat, nrow(feat_cors),
              feat_cors$Gene[1], feat_cors$Pearson_r[1]))
}

}  # end if (run_partB)


# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  FUNCTIONAL & META-ANALYSIS SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("piRNA Signature:", paste(top_feats, collapse = ", "), "\n\n")

cat("PART A \u2014 Meta-Analysis:\n")
cat("  A1. Expression heatmap across", length(datasets), "datasets\n")
cat("  A2. SMD forest plots for", length(meta_results_all), "piRNAs\n")
cat("  A3. Spearman correlation heatmaps (all + averaged)\n\n")

cat("PART B \u2014 Functional Prediction (piRNA\u2013mRNA Correlation):\n")
if (run_partB) {
  cat("  B1. mRNA data loaded for piRNA\u2013mRNA correlation\n")
  cat("  B2. Pearson correlations:", nrow(sig_cors), "significant piRNA\u2013mRNA pairs\n")
  cat("      Dot plot: results/functional/piRNA_mRNA_correlation_dotplot.png\n")
  if (exists("kegg_res") && !is.null(kegg_res)) {
    n_kegg <- nrow(kegg_res@result[kegg_res@result$p.adjust < 0.05, ])
    cat("  B3. KEGG pathways:", n_kegg, "significant\n")
  }
  cat("  B3. GO enrichment: BP, CC, MF bubble plots\n")
  if (exists("reactome_sig") && !is.null(reactome_sig)) {
    cat("  B3. Reactome pathways:", nrow(reactome_sig), "significant\n")
  }
  cat("  B4. Functional network: piRNA \u2192 target genes (circular layout)\n")
  cat("  B4b. piRNA\u2013mRNA scatter plots with marginal histograms (top 3 genes/piRNA)\n")
  cat("  B5. Per-piRNA correlated mRNA gene tables\n")
} else {
  cat("  (Skipped \u2014 mRNA expression data not available)\n")
}

cat("\nOutput directories:\n")
cat("  results/meta_analysis/  \u2014 Heatmaps, SMD forest plots, correlation\n")
cat("  results/functional/     \u2014 piRNA\u2013mRNA correlations, KEGG, GO, Reactome\n")
cat("  results/network/        \u2014 Functional interaction network\n")
cat("  results/functional/scatter_plots/ \u2014 piRNA\u2013mRNA scatter + marginal histograms\n")

end_time_fm <- Sys.time()
cat("\nRuntime:", round(difftime(end_time_fm, start_time_fm, units = "mins"), 1), "min\n")
cat("\n*** FUNCTIONAL & META-ANALYSIS COMPLETE ***\n")

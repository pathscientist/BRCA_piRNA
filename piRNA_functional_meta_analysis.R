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
  "gridExtra", "cowplot", "reshape2", "corrplot"
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
            comb.fixed = FALSE, comb.random = TRUE)
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
      plot.margin = margin(10, 100, 10, 10),
      axis.text.y = element_text(
        face = ifelse(plot_df$Dataset[order(plot_df$y_pos)] == "Overall",
                      "bold", "plain"),
        size = 10
      )
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
# ══════════════════════════════════════════════════════════════════════════════

cat("\n")
cat(paste(rep("═", 70), collapse = ""), "\n")
cat("  PART B: FUNCTIONAL PREDICTION\n")
cat(paste(rep("═", 70), collapse = ""), "\n")


# ==============================================================================
# B1. PEARSON'S CORRELATION: piRNAs vs ALL GENES
#     Identify correlated genes for each signature piRNA
# ==============================================================================
cat("\n========== B1: Pearson's Correlation (piRNAs vs Genes) ==========\n")

# Strategy: correlate each signature piRNA against all other piRNAs/genes
# in the expression matrix. If mRNA data is available, it can be loaded
# separately and joined.

# --- Option A: Correlate within piRNA expression data ---
# (Use all non-signature piRNAs as potential target correlates)
non_sig_genes <- setdiff(gene_cols, top_feats)
cat("  Correlating", length(top_feats), "signature piRNAs against",
    length(non_sig_genes), "other features...\n")

# Use only tumor samples for biologically relevant correlations
tumor_idx <- combat_df_all$Group == "Tumor"
expr_sig    <- as.matrix(combat_df_all[tumor_idx, top_feats])
expr_others <- as.matrix(combat_df_all[tumor_idx, non_sig_genes])

# Compute Pearson correlations: each signature piRNA vs all other genes
cor_results <- list()
for (feat in top_feats) {
  cors <- apply(expr_others, 2, function(col) {
    ct <- cor.test(expr_sig[, feat], col, method = "pearson")
    c(r = ct$estimate, p = ct$p.value)
  })
  cor_df <- data.frame(
    piRNA      = feat,
    Gene       = colnames(cors),
    Pearson_r  = cors["r.cor", ],
    P_value    = cors["p", ],
    stringsAsFactors = FALSE
  )
  # Adjust p-values
  cor_df$P_adj <- p.adjust(cor_df$P_value, method = "BH")
  # Filter significant correlations (|r| > 0.3 and adj.p < 0.05)
  cor_df$Significant <- abs(cor_df$Pearson_r) > 0.3 & cor_df$P_adj < 0.05
  cor_results[[feat]] <- cor_df
}

all_cors <- do.call(rbind, cor_results)
sig_cors <- all_cors[all_cors$Significant, ]

cat(sprintf("  Total significant correlations: %d (|r|>0.3, adj.p<0.05)\n",
            nrow(sig_cors)))
for (feat in top_feats) {
  n_sig <- sum(sig_cors$piRNA == feat)
  cat(sprintf("    %s: %d correlated genes\n", feat, n_sig))
}

# Save full correlation table
write.csv(all_cors, "results/functional/pearson_correlations_all.csv", row.names = FALSE)
write.csv(sig_cors, "results/functional/pearson_correlations_significant.csv", row.names = FALSE)

# --- Option B: Load mRNA expression data if available ---
# >>> EDIT THIS: If you have paired mRNA data <<<
# mRNA_expr <- read.table("mRNA_expression.txt", header=TRUE, row.names=1, sep="\t")
# Then replace expr_others above with t(mRNA_expr) to correlate piRNAs vs mRNAs

# --- Correlation plot: top correlated genes per piRNA ---
cat("\n  Generating correlation summary plot...\n")

# Top N correlated genes per piRNA
top_n_corr <- 10
top_cors <- sig_cors %>%
  group_by(piRNA) %>%
  arrange(desc(abs(Pearson_r))) %>%
  slice_head(n = top_n_corr) %>%
  ungroup()

if (nrow(top_cors) > 0) {
  p_cor_bar <- ggplot(top_cors,
                      aes(x = reorder(Gene, abs(Pearson_r)),
                          y = Pearson_r,
                          fill = ifelse(Pearson_r > 0, "Positive", "Negative"))) +
    geom_col(width = 0.7, alpha = 0.85) +
    coord_flip() +
    facet_wrap(~ piRNA, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = c("Positive" = "#2166AC", "Negative" = "#B2182B"),
                      name = "Correlation") +
    labs(
      title = "Top Correlated Genes per Signature piRNA",
      subtitle = "Pearson's correlation (tumor samples, adj.p < 0.05)",
      x = "", y = "Pearson's r"
    ) +
    pub_theme +
    theme(
      strip.background = element_rect(fill = "grey95"),
      legend.position = "top"
    )

  ggsave("results/functional/correlation_top_genes.png",
         p_cor_bar, width = 12,
         height = max(6, ceiling(length(top_feats) / 2) * 4),
         dpi = 300)
  cat("  Correlation bar plot saved.\n")
}


# ==============================================================================
# B2. GENE SET ENRICHMENT: KEGG PATHWAYS
# ==============================================================================
cat("\n========== B2: KEGG Pathway Enrichment ==========\n")

# Collect all significant correlated genes across piRNAs
correlated_genes <- unique(sig_cors$Gene)
cat("  Unique correlated genes for enrichment:", length(correlated_genes), "\n")

# Attempt to map gene names to Entrez IDs
# piRNA names may not map to standard gene names; we try anyway
# If you have mRNA gene names, this will work well

# Try mapping gene names (removing common piRNA prefixes)
gene_names_clean <- gsub("^hsa[-_]", "", correlated_genes)
gene_names_clean <- gsub("^piR[-_]", "", gene_names_clean)

# Try converting to Entrez IDs
entrez_map <- tryCatch({
  bitr(gene_names_clean, fromType = "SYMBOL", toType = "ENTREZID",
       OrgDb = org.Hs.eg.db)
}, error = function(e) {
  cat("  Note: Gene symbol mapping failed. Trying ALIAS...\n")
  tryCatch(
    bitr(gene_names_clean, fromType = "ALIAS", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db),
    error = function(e2) NULL
  )
})

# If piRNA names cannot be mapped, use all genes as background
# and create per-piRNA gene lists from correlation analysis
run_enrichment <- !is.null(entrez_map) && nrow(entrez_map) >= 5

if (run_enrichment) {
  entrez_ids <- unique(entrez_map$ENTREZID)
  cat("  Mapped", length(entrez_ids), "genes to Entrez IDs.\n")

  # --- KEGG Enrichment ---
  kegg_res <- tryCatch({
    enrichKEGG(gene = entrez_ids, organism = "hsa",
               pvalueCutoff = 0.05, qvalueCutoff = 0.1,
               minGSSize = 5, maxGSSize = 500)
  }, error = function(e) {
    cat("  KEGG enrichment failed:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(kegg_res) && nrow(kegg_res@result[kegg_res@result$p.adjust < 0.05, ]) > 0) {
    kegg_sig <- kegg_res@result[kegg_res@result$p.adjust < 0.05, ]
    write.csv(kegg_sig, "results/functional/kegg_enrichment.csv", row.names = FALSE)

    # KEGG Bar Plot
    n_show <- min(20, nrow(kegg_sig))
    kegg_plot_df <- head(kegg_sig, n_show) %>%
      mutate(
        NegLogP = -log10(pvalue),
        Count = as.numeric(Count),
        Description = factor(Description, levels = rev(Description))
      )

    p_kegg <- ggplot(kegg_plot_df, aes(x = Count, y = Description, fill = pvalue)) +
      geom_col(width = 0.7) +
      scale_fill_gradient(low = "grey40", high = "#8FBC8F",
                          name = "p-value",
                          limits = c(0.01, 0.05),
                          oob = squish,
                          breaks = c(0.01, 0.03, 0.05)) +
      labs(
        title = "KEGG Pathway Enrichment Analysis",
        subtitle = paste0(length(entrez_ids), " correlated genes"),
        x = "Number of Enriched Genes",
        y = ""
      ) +
      pub_theme +
      theme(axis.text.y = element_text(size = 9))

    ggsave("results/functional/kegg_barplot.png",
           p_kegg, width = 10, height = max(5, n_show * 0.35 + 2), dpi = 300)
    cat("  KEGG bar plot saved (", nrow(kegg_sig), " significant pathways).\n")
  } else {
    cat("  No significant KEGG pathways found.\n")
  }

  # ==============================================================================
  # B3. GO ENRICHMENT (BP, CC, MF) — BUBBLE PLOTS
  # ==============================================================================
  cat("\n========== B3: GO Enrichment (BP, CC, MF) ==========\n")

  go_categories <- c("BP", "CC", "MF")
  go_full_names <- c(BP = "Biological Process", CC = "Cellular Component",
                     MF = "Molecular Function")

  for (ont in go_categories) {
    cat(sprintf("  Running GO %s enrichment...\n", ont))

    go_res <- tryCatch({
      enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db,
               ont = ont, pvalueCutoff = 0.05, qvalueCutoff = 0.1,
               readable = TRUE, minGSSize = 5, maxGSSize = 500)
    }, error = function(e) {
      cat(sprintf("    GO %s failed: %s\n", ont, conditionMessage(e)))
      NULL
    })

    if (!is.null(go_res) && nrow(go_res@result[go_res@result$p.adjust < 0.05, ]) > 0) {
      go_sig <- go_res@result[go_res@result$p.adjust < 0.05, ]
      write.csv(go_sig, paste0("results/functional/go_", tolower(ont), "_enrichment.csv"),
                row.names = FALSE)

      # Bubble plot
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
        scale_color_gradient(low = "grey40", high = "#8FBC8F",
                             name = "FDR",
                             limits = c(0.01, 0.05), oob = squish,
                             breaks = c(0.01, 0.03, 0.05)) +
        labs(
          title = paste0("GO Enrichment: ", go_full_names[ont]),
          x = "Gene Ratio",
          y = ""
        ) +
        pub_theme +
        theme(axis.text.y = element_text(size = 9))

      ggsave(paste0("results/functional/go_", tolower(ont), "_bubble.png"),
             p_go, width = 10, height = max(5, n_show * 0.35 + 2), dpi = 300)
      cat(sprintf("    GO %s: %d significant terms, bubble plot saved.\n",
                  ont, nrow(go_sig)))
    } else {
      cat(sprintf("    No significant GO %s terms found.\n", ont))
    }
  }

  # ==============================================================================
  # B4. REACTOME PATHWAY ENRICHMENT
  # ==============================================================================
  cat("\n========== B4: Reactome Pathway Enrichment ==========\n")

  reactome_res <- tryCatch({
    enrichPathway(gene = entrez_ids, organism = "human",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.1,
                  readable = TRUE, minGSSize = 5, maxGSSize = 500)
  }, error = function(e) {
    cat("  Reactome enrichment failed:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(reactome_res) &&
      nrow(reactome_res@result[reactome_res@result$p.adjust < 0.05, ]) > 0) {
    react_sig <- reactome_res@result[reactome_res@result$p.adjust < 0.05, ]
    write.csv(react_sig, "results/functional/reactome_enrichment.csv", row.names = FALSE)

    n_show <- min(15, nrow(react_sig))
    react_plot_df <- head(react_sig, n_show) %>%
      mutate(
        Count = as.numeric(Count),
        Description = factor(Description, levels = rev(Description))
      )

    p_react <- ggplot(react_plot_df,
                      aes(x = Count, y = Description, fill = p.adjust)) +
      geom_col(width = 0.7) +
      scale_fill_gradient(low = "grey40", high = "#8FBC8F",
                          name = "FDR", limits = c(0.01, 0.05), oob = squish) +
      labs(
        title = "Reactome Pathway Enrichment",
        x = "Number of Enriched Genes", y = ""
      ) +
      pub_theme +
      theme(axis.text.y = element_text(size = 9))

    ggsave("results/functional/reactome_barplot.png",
           p_react, width = 11, height = max(5, n_show * 0.35 + 2), dpi = 300)
    cat("  Reactome: ", nrow(react_sig), " significant pathways, plot saved.\n")
  } else {
    cat("  No significant Reactome pathways found.\n")
  }

} else {
  cat("  Gene name mapping was insufficient for enrichment analysis.\n")
  cat("  To enable enrichment, provide mRNA expression data with standard gene symbols.\n")
  cat("  Skipping KEGG, GO, and Reactome enrichment.\n\n")

  cat("  >>> ALTERNATIVE APPROACH FOR piRNA FUNCTIONAL PREDICTION <<<\n")
  cat("  Since piRNA names may not map directly to gene ontologies,\n")
  cat("  consider these strategies:\n")
  cat("  1. Use piRBase (http://www.regulatoryrna.org/database/piRNA/) to find\n")
  cat("     piRNA target genes, then run enrichment on those targets.\n")
  cat("  2. Use miRanda or piRNAPredictor for piRNA-mRNA interaction prediction.\n")
  cat("  3. Provide paired mRNA expression data for correlation-based prediction.\n")
  cat("  4. Load pre-computed piRNA target genes from databases.\n")
}


# ==============================================================================
# B5. FUNCTIONAL INTERACTION NETWORK
#     (piRNA → correlated genes → KEGG pathways)
# ==============================================================================
cat("\n========== B5: Functional Interaction Network ==========\n")

# Build network from correlation results and enrichment results
build_pirna_network <- function(sig_cors, kegg_sig = NULL, top_genes = 5) {
  edges <- data.frame()
  nodes <- data.frame()

  # Node: piRNAs (diamond shape)
  pirna_nodes <- data.frame(
    name  = unique(sig_cors$piRNA),
    type  = "piRNA",
    group = "piRNA",
    stringsAsFactors = FALSE
  )

  # Top correlated genes per piRNA
  top_per_pirna <- sig_cors %>%
    group_by(piRNA) %>%
    arrange(desc(abs(Pearson_r))) %>%
    slice_head(n = top_genes) %>%
    ungroup()

  # Node: genes (ellipse shape)
  gene_nodes <- data.frame(
    name  = unique(top_per_pirna$Gene),
    type  = "mRNA",
    group = "gene",
    stringsAsFactors = FALSE
  )

  # Edges: piRNA → gene (correlation)
  pirna_gene_edges <- top_per_pirna %>%
    transmute(
      from   = piRNA,
      to     = Gene,
      weight = abs(Pearson_r),
      type   = ifelse(Pearson_r > 0, "positive", "negative")
    )

  # If KEGG results available, add pathway nodes
  pathway_nodes <- data.frame()
  gene_pathway_edges <- data.frame()

  if (!is.null(kegg_sig) && nrow(kegg_sig) > 0) {
    top_pathways <- head(kegg_sig, 8)  # Top 8 pathways

    pathway_nodes <- data.frame(
      name  = top_pathways$Description,
      type  = "pathway",
      group = top_pathways$Description,
      stringsAsFactors = FALSE
    )

    # Parse gene lists in KEGG results to build gene → pathway edges
    for (i in seq_len(nrow(top_pathways))) {
      pathway_genes <- unlist(strsplit(top_pathways$geneID[i], "/"))
      # Match with our correlated genes
      matched <- intersect(pathway_genes, gene_nodes$name)
      if (length(matched) > 0) {
        gene_pathway_edges <- rbind(gene_pathway_edges, data.frame(
          from   = matched,
          to     = top_pathways$Description[i],
          weight = 1,
          type   = "pathway",
          stringsAsFactors = FALSE
        ))
        # Color genes by pathway
        gene_nodes$group[gene_nodes$name %in% matched] <- top_pathways$Description[i]
      }
    }
  }

  # Combine
  all_nodes <- rbind(pirna_nodes, gene_nodes, pathway_nodes)
  all_edges <- rbind(pirna_gene_edges, gene_pathway_edges)

  list(nodes = all_nodes, edges = all_edges)
}

if (nrow(sig_cors) > 0) {
  kegg_sig_for_net <- NULL
  if (exists("kegg_res") && !is.null(kegg_res)) {
    kegg_sig_for_net <- kegg_res@result[kegg_res@result$p.adjust < 0.05, ]
  }

  net_data <- build_pirna_network(sig_cors, kegg_sig_for_net, top_genes = 5)

  if (nrow(net_data$edges) > 0) {
    # Build igraph object
    g <- graph_from_data_frame(net_data$edges, directed = FALSE,
                               vertices = net_data$nodes)

    # Set visual attributes
    V(g)$shape <- ifelse(V(g)$type == "piRNA", "diamond",
                  ifelse(V(g)$type == "pathway", "square", "circle"))

    # Assign colors by group
    unique_groups <- unique(V(g)$group)
    group_colors <- c(
      "piRNA" = "#E41A1C",
      "gene"  = "#377EB8",
      setNames(
        colorRampPalette(c("#4DAF4A", "#FF7F00", "#984EA3", "#A65628",
                           "#F781BF", "#999999", "#66C2A5", "#FC8D62"))(
          max(1, length(unique_groups) - 2)),
        setdiff(unique_groups, c("piRNA", "gene"))
      )
    )

    # Plot with ggraph
    set.seed(SEED)
    p_net <- ggraph(g, layout = "fr") +
      geom_edge_link(aes(edge_alpha = weight,
                         edge_linetype = ifelse(type == "negative", "dashed", "solid")),
                     edge_colour = "grey60", edge_width = 0.5,
                     show.legend = FALSE) +
      geom_node_point(aes(color = group, shape = type), size = 5, alpha = 0.9) +
      geom_node_text(aes(label = name), repel = TRUE, size = 2.8,
                     max.overlaps = 20) +
      scale_shape_manual(values = c("piRNA" = 18, "mRNA" = 16, "pathway" = 15),
                         name = "Node Type") +
      scale_color_manual(values = group_colors, name = "Group") +
      labs(
        title = "Functional Interaction Network",
        subtitle = "piRNAs — Correlated Genes — KEGG Pathways"
      ) +
      theme_void() +
      theme(
        plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
        legend.position = "right",
        legend.text = element_text(size = 8)
      )

    ggsave("results/network/functional_network.png",
           p_net, width = 14, height = 10, dpi = 300)
    cat("  Functional interaction network saved.\n")

    # Also save as edge list
    write.csv(net_data$edges, "results/network/network_edges.csv", row.names = FALSE)
    write.csv(net_data$nodes, "results/network/network_nodes.csv", row.names = FALSE)
  } else {
    cat("  Not enough data to build network.\n")
  }
} else {
  cat("  No significant correlations found. Network skipped.\n")
}


# ==============================================================================
# B6. PER-piRNA SUMMARY: CORRELATED GENES TABLE
# ==============================================================================
cat("\n========== B6: Per-piRNA Summary ==========\n")

for (feat in top_feats) {
  feat_cors <- sig_cors[sig_cors$piRNA == feat, ]
  if (nrow(feat_cors) == 0) {
    cat(sprintf("  %s: No significant correlated genes.\n", feat))
    next
  }

  feat_cors <- feat_cors %>% arrange(desc(abs(Pearson_r)))
  fname <- gsub("[^a-zA-Z0-9_.-]", "_", feat)
  write.csv(feat_cors,
            paste0("results/functional/correlated_genes_", fname, ".csv"),
            row.names = FALSE)

  cat(sprintf("  %s: %d correlated genes (top: %s, r=%.3f)\n",
              feat, nrow(feat_cors),
              feat_cors$Gene[1], feat_cors$Pearson_r[1]))
}


# ==============================================================================
# SUMMARY
# ==============================================================================
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("  FUNCTIONAL & META-ANALYSIS SUMMARY\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("piRNA Signature:", paste(top_feats, collapse = ", "), "\n\n")

cat("PART A — Meta-Analysis:\n")
cat("  A1. Expression heatmap across", length(datasets), "datasets\n")
cat("  A2. SMD forest plots for", length(meta_results_all), "piRNAs\n")
cat("  A3. Spearman correlation heatmaps (all + averaged)\n\n")

cat("PART B — Functional Prediction:\n")
cat("  B1. Pearson correlations:", nrow(sig_cors), "significant gene-piRNA pairs\n")
if (exists("kegg_res") && !is.null(kegg_res)) {
  n_kegg <- nrow(kegg_res@result[kegg_res@result$p.adjust < 0.05, ])
  cat("  B2. KEGG pathways:", n_kegg, "significant\n")
}
cat("  B3. GO enrichment: BP, CC, MF bubble plots\n")
if (exists("reactome_res") && !is.null(reactome_res)) {
  n_react <- nrow(reactome_res@result[reactome_res@result$p.adjust < 0.05, ])
  cat("  B4. Reactome pathways:", n_react, "significant\n")
}
cat("  B5. Functional interaction network\n")
cat("  B6. Per-piRNA correlated gene tables\n\n")

cat("Output directories:\n")
cat("  results/meta_analysis/  — Heatmaps, SMD forest plots, correlation\n")
cat("  results/functional/     — KEGG, GO, Reactome, correlation tables\n")
cat("  results/network/        — Interaction network + edge/node lists\n")

end_time_fm <- Sys.time()
cat("\nRuntime:", round(difftime(end_time_fm, start_time_fm, units = "mins"), 1), "min\n")
cat("\n*** FUNCTIONAL & META-ANALYSIS COMPLETE ***\n")

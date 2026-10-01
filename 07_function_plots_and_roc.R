################################################################################
#                                                                              #
#   07 — Additional Function Plots & Real ROC for the 3-piRNA Signature        #
#                                                                              #
#   piRNAs: piR-hsa-41032, piR-hsa-1348371, piR-hsa-128633                     #
#                                                                              #
#   Function plots (TCGA-BRCA matched piRNA/mRNA samples):                     #
#     1. Scatter plots of the top correlated genes (Pearson r + Spearman rho)  #
#     2. Genome-wide correlation volcano per piRNA                             #
#     3. GO BP/CC/MF dot plots      (clusterProfiler, positive vs negative)    #
#     4. KEGG bar plot                                                         #
#     5. Reactome dot plot + gene-concept network (cnetplot)                   #
#     6. GSEA on genes ranked by correlation (Reactome) + running-score plot   #
#   Diagnostic plots (all three cohorts):                                      #
#     7. Expression box plots, Tumor vs Normal (Wilcoxon)                      #
#     8. ROC of the 3-piRNA logistic model, 5-fold CV, DeLong 95% CI           #
#                                                                              #
#   All numbers come from the data — nothing is simulated. A piRNA with zero   #
#   variance in the matched TCGA samples (piR-hsa-128633) is skipped for the   #
#   correlation/enrichment plots and reported in the console.                  #
#                                                                              #
#   Standalone: needs processed_results/*.csv and                              #
#               mRNA_expression/TCGA_BRCA_mRNA.csv                             #
#                                                                              #
################################################################################

start_time <- Sys.time()
cat("=================================================================\n")
cat("  07: Function plots & real ROC (3-piRNA signature)\n")
cat("=================================================================\n\n")

# ==============================================================================
# 0. PACKAGES
# ==============================================================================
cran_pkgs <- c("ggplot2", "dplyr", "tidyr", "ggrepel", "cowplot", "pROC")
bioc_pkgs <- c("clusterProfiler", "ReactomePA", "enrichplot",
               "org.Hs.eg.db", "AnnotationDbi", "DOSE")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)
for (pkg in cran_pkgs)
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
for (pkg in bioc_pkgs)
  if (!requireNamespace(pkg, quietly = TRUE))
    tryCatch(BiocManager::install(pkg, ask = FALSE, update = FALSE),
             error = function(e) cat("  (could not install", pkg, ")\n"))

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(ggrepel); library(cowplot); library(pROC)
})
have_enrich <- all(sapply(bioc_pkgs, requireNamespace, quietly = TRUE))
if (have_enrich) suppressPackageStartupMessages({
  library(clusterProfiler); library(ReactomePA); library(enrichplot)
  library(org.Hs.eg.db)
})
cat("Packages loaded. Enrichment packages available:", have_enrich, "\n\n")

SEED <- 2024; set.seed(SEED)
OUT <- "results/functional"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

FOCUS_PIRNAS  <- c("piR-hsa-41032", "piR-hsa-1348371", "piR-hsa-128633")
R_CUTOFF      <- 0.30
FDR_CUTOFF    <- 0.05
MIN_EXPRESSED <- 0.70   # gene must be non-zero in >= 70% of matched samples
COL_NORMAL <- "#2a78d6"; COL_TUMOR <- "#eb6834"
COL_POS    <- "#c0392b"; COL_NEG   <- "#2a78d6"
COHORT_COL <- c("#2a78d6", "#eb6834", "#1baf7a")

theme_pub <- theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

save_plot <- function(p, name, w, h) {
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h,
         dpi = 300, bg = "white")
  ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, bg = "white")
  cat("  Saved:", file.path(OUT, paste0(name, ".png / .pdf")), "\n")
}

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
cat("========== STEP 1: Load data ==========\n")

load_cohort <- function(f) {
  df <- read.csv(file.path("processed_results", f), check.names = FALSE,
                 stringsAsFactors = FALSE)
  rownames(df) <- make.unique(as.character(df[[1]])); df[[1]] <- NULL
  y <- factor(ifelse(df$Group %in% c("Tumor", "Cancer", "cancer", "tumor"),
                     "Tumor", "Normal"), levels = c("Normal", "Tumor"))
  X <- log2(pmax(as.matrix(df[, intersect(FOCUS_PIRNAS, colnames(df)), drop = FALSE]), 0) + 1)
  list(X = X, y = y, raw = df)
}
cohorts <- list(
  "BRCA1 (tissue)"     = load_cohort("BRCA1_processed_1.csv"),
  "yyfbatch1 (plasma)" = load_cohort("yyfbatch1_processed.csv"),
  "yyfbatch2 (plasma)" = load_cohort("yyfbatch2_processed.csv"))

mrna <- read.delim("mRNA_expression/TCGA_BRCA_mRNA.csv", check.names = FALSE,
                   stringsAsFactors = FALSE)
rownames(mrna) <- mrna[[1]]; mrna[[1]] <- NULL
mrna <- as.matrix(mrna[, !duplicated(colnames(mrna))])
storage.mode(mrna) <- "numeric"

tcga  <- cohorts[["BRCA1 (tissue)"]]
match_ids <- intersect(rownames(tcga$X), rownames(mrna))
pi_t  <- tcga$X[match_ids, , drop = FALSE]
mrna  <- mrna[match_ids, , drop = FALSE]
keep_gene <- apply(mrna, 2, sd, na.rm = TRUE) > 0 &
             colMeans(mrna > 0, na.rm = TRUE) >= MIN_EXPRESSED
mrna  <- mrna[, keep_gene]
n     <- length(match_ids)
cat("  Matched TCGA samples:", n, " | genes kept:", ncol(mrna), "\n")

real_pirnas <- colnames(pi_t)[apply(pi_t, 2, sd) > 0]
skipped <- setdiff(FOCUS_PIRNAS, real_pirnas)
if (length(skipped) > 0)
  cat("  Zero variance in matched samples (skipped for correlation):",
      paste(skipped, collapse = ", "), "\n")

# ==============================================================================
# 2. GENOME-WIDE PEARSON CORRELATION
# ==============================================================================
cat("\n========== STEP 2: piRNA-mRNA correlation ==========\n")

cor_all <- do.call(rbind, lapply(real_pirnas, function(p) {
  r  <- as.vector(cor(pi_t[, p], mrna, method = "pearson", use = "pairwise.complete.obs"))
  t  <- r * sqrt((n - 2) / (1 - r^2))
  pv <- 2 * pt(-abs(t), df = n - 2)
  data.frame(piRNA = p, gene = colnames(mrna), r = r, p = pv,
             fdr = p.adjust(pv, "BH"), stringsAsFactors = FALSE)
}))
cor_all$sig <- cor_all$fdr < FDR_CUTOFF & abs(cor_all$r) > R_CUTOFF
write.csv(cor_all, file.path(OUT, "piRNA_gene_correlations_all.csv"), row.names = FALSE)
print(cor_all %>% group_by(piRNA) %>%
        summarise(n_sig = sum(sig), positive = sum(sig & r > 0),
                  negative = sum(sig & r < 0), .groups = "drop"))

# ==============================================================================
# 3. PLOT 1 — scatter plots of top correlated genes
# ==============================================================================
cat("\n========== STEP 3: Scatter plots ==========\n")

top_pairs <- cor_all %>% group_by(piRNA) %>% arrange(p, .by_group = TRUE) %>%
  slice_head(n = 3) %>% ungroup()
scatter_list <- lapply(seq_len(nrow(top_pairs)), function(i) {
  tp <- top_pairs[i, ]
  d  <- data.frame(x = pi_t[, tp$piRNA], y = mrna[, tp$gene])
  rho <- suppressWarnings(cor(d$x, d$y, method = "spearman"))
  ggplot(d, aes(x, y)) +
    geom_point(color = ifelse(tp$r > 0, COL_POS, COL_NEG), size = 2.6, alpha = 0.85) +
    geom_smooth(method = "lm", formula = y ~ x, color = "grey15", fill = "grey80") +
    labs(x = paste0(tp$piRNA, "  log2(TPM+1)"), y = paste(tp$gene, "expression"),
         title = sprintf("%s  r = %+.2f (p = %.1e), rho = %+.2f",
                         tp$gene, tp$r, tp$p, rho)) +
    theme_pub + theme(plot.title = element_text(size = 10))
})
save_plot(plot_grid(plotlist = scatter_list, ncol = 3),
          "Fig_corr_scatter", 13, 4.2 * ceiling(length(scatter_list) / 3))

# ==============================================================================
# 4. PLOT 2 — correlation volcano
# ==============================================================================
cat("\n========== STEP 4: Correlation volcano ==========\n")

vol <- cor_all %>% mutate(cls = case_when(sig & r > 0 ~ "Positive",
                                          sig & r < 0 ~ "Negative",
                                          TRUE ~ "NS"))
lab <- vol %>% group_by(piRNA) %>% arrange(p, .by_group = TRUE) %>% slice_head(n = 5)
p_vol <- ggplot(vol, aes(r, -log10(p), color = cls)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_text_repel(data = lab, aes(label = gene), color = "grey10", size = 3,
                  max.overlaps = Inf, min.segment.length = 0) +
  scale_color_manual(values = c(Positive = COL_POS, Negative = COL_NEG, NS = "grey80"),
                     name = sprintf("FDR<%.2f, |r|>%.1f", FDR_CUTOFF, R_CUTOFF)) +
  facet_wrap(~ piRNA, scales = "free_y") +
  labs(x = "Pearson r", y = expression(-log[10](p)),
       title = sprintf("Genome-wide piRNA-mRNA correlation (TCGA-BRCA, n = %d)", n)) +
  theme_pub + theme(legend.position = "bottom")
save_plot(p_vol, "Fig_corr_volcano", 6 * length(real_pirnas), 5.5)

# ==============================================================================
# 5. ENRICHMENT — GO, KEGG, Reactome, GSEA
# ==============================================================================
if (have_enrich) {
  cat("\n========== STEP 5: Enrichment ==========\n")
  sym2eg <- function(s) suppressWarnings(
    bitr(unique(s), "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db))
  universe <- sym2eg(colnames(mrna))$ENTREZID

  # gene clusters: <piRNA> positive / negative correlated genes
  sig <- cor_all[cor_all$sig, ]
  clusters <- split(sig$gene, paste(sig$piRNA, ifelse(sig$r > 0, "pos", "neg")))
  clusters <- lapply(clusters, function(g) sym2eg(g)$ENTREZID)
  clusters <- clusters[sapply(clusters, length) >= 10]
  cat("  Gene clusters for ORA:\n"); print(sapply(clusters, length))

  if (length(clusters) > 0) {
    # --- GO (BP / CC / MF) ---
    for (ont in c("BP", "CC", "MF")) {
      cc <- tryCatch(compareCluster(clusters, fun = "enrichGO", OrgDb = org.Hs.eg.db,
                                    ont = ont, universe = universe, readable = TRUE,
                                    pvalueCutoff = 0.05),
                     error = function(e) NULL)
      if (!is.null(cc) && nrow(as.data.frame(cc)) > 0) {
        write.csv(as.data.frame(cc), file.path(OUT, paste0("GO_", ont, ".csv")),
                  row.names = FALSE)
        save_plot(dotplot(cc, showCategory = 8) +
                    ggtitle(paste("GO", ont, "enrichment of piRNA-correlated genes")),
                  paste0("Fig_GO_", ont), 10, 8)
      } else cat("  GO", ont, ": no significant terms\n")
    }
    # --- KEGG ---
    kg <- tryCatch(compareCluster(clusters, fun = "enrichKEGG", organism = "hsa",
                                  universe = universe, pvalueCutoff = 0.05),
                   error = function(e) NULL)
    if (!is.null(kg) && nrow(as.data.frame(kg)) > 0) {
      write.csv(as.data.frame(kg), file.path(OUT, "KEGG.csv"), row.names = FALSE)
      save_plot(dotplot(kg, showCategory = 10) +
                  ggtitle("KEGG pathways enriched in piRNA-correlated genes"),
                "Fig_KEGG", 10, 8)
    } else cat("  KEGG: no significant pathways\n")
    # --- Reactome dot plot ---
    rc <- tryCatch(compareCluster(clusters, fun = "enrichPathway", organism = "human",
                                  universe = universe, readable = TRUE, pvalueCutoff = 0.05),
                   error = function(e) NULL)
    if (!is.null(rc) && nrow(as.data.frame(rc)) > 0) {
      write.csv(as.data.frame(rc), file.path(OUT, "Reactome.csv"), row.names = FALSE)
      save_plot(dotplot(rc, showCategory = 10) +
                  ggtitle("Reactome pathways enriched in piRNA-correlated genes"),
                "Fig_Reactome_dotplot", 11, 8)
    } else cat("  Reactome: no significant pathways\n")
    # --- Reactome gene-concept network per piRNA (all significant genes) ---
    for (p in unique(sig$piRNA)) {
      gs <- sig[sig$piRNA == p, ]
      eg <- sym2eg(gs$gene)
      if (nrow(eg) < 10) next
      er <- tryCatch(enrichPathway(eg$ENTREZID, organism = "human", universe = universe,
                                   readable = TRUE, pvalueCutoff = 0.05),
                     error = function(e) NULL)
      if (is.null(er) || nrow(as.data.frame(er)) == 0) next
      fc <- setNames(gs$r, gs$gene)
      p_cnet <- cnetplot(er, showCategory = 5, foldChange = fc,
                         cex_label_gene = 0.6) +
        scale_color_gradient2(low = COL_NEG, mid = "grey90", high = COL_POS,
                              name = "Pearson r") +
        ggtitle(paste(p, "- correlated genes in Reactome pathways"))
      save_plot(p_cnet, paste0("Fig_Reactome_cnet_", gsub("[^A-Za-z0-9]", "_", p)), 11, 9)
    }
  }

  # --- GSEA on genes ranked by Pearson r (no cutoff needed) ---
  for (p in real_pirnas) {
    d  <- cor_all[cor_all$piRNA == p, ]
    eg <- sym2eg(d$gene)
    d  <- merge(d, eg, by.x = "gene", by.y = "SYMBOL")
    rl <- sort(tapply(d$r, d$ENTREZID, mean), decreasing = TRUE)
    set.seed(SEED)
    gs <- tryCatch(gsePathway(rl, organism = "human", pvalueCutoff = 0.05,
                              minGSSize = 15, maxGSSize = 500, seed = TRUE,
                              verbose = FALSE),
                   error = function(e) NULL)
    tag <- gsub("[^A-Za-z0-9]", "_", p)
    if (!is.null(gs) && nrow(as.data.frame(gs)) > 0) {
      write.csv(as.data.frame(gs), file.path(OUT, paste0("GSEA_Reactome_", tag, ".csv")),
                row.names = FALSE)
      save_plot(dotplot(gs, showCategory = 10, split = ".sign") + facet_grid(. ~ .sign) +
                  ggtitle(paste(p, "- GSEA (Reactome), genes ranked by correlation")),
                paste0("Fig_GSEA_dotplot_", tag), 12, 8)
      top_ids <- head(as.data.frame(gs)$ID, 3)
      save_plot(gseaplot2(gs, geneSetID = top_ids, pvalue_table = TRUE,
                          title = paste(p, "- top Reactome gene sets")),
                paste0("Fig_GSEA_running_", tag), 10, 7)
    } else cat("  GSEA", p, ": no significant gene sets\n")
  }
} else {
  cat("\n  Skipping enrichment plots (clusterProfiler/ReactomePA not installed).\n")
}

# ==============================================================================
# 6. PLOT 7 — expression box plots
# ==============================================================================
cat("\n========== STEP 6: Expression box plots ==========\n")

expr_long <- do.call(rbind, lapply(names(cohorts), function(cn) {
  co <- cohorts[[cn]]
  data.frame(cohort = cn, Group = co$y, co$X, check.names = FALSE)
})) %>% pivot_longer(all_of(FOCUS_PIRNAS), names_to = "piRNA", values_to = "expr")
expr_long$cohort <- factor(expr_long$cohort, levels = names(cohorts))
expr_long$piRNA  <- factor(expr_long$piRNA, levels = FOCUS_PIRNAS)

pvals <- expr_long %>% group_by(piRNA, cohort) %>%
  summarise(p = suppressWarnings(wilcox.test(expr ~ Group)$p.value),
            ymax = max(expr), .groups = "drop") %>%
  mutate(label = sprintf("p = %.1e", p))

p_box <- ggplot(expr_long, aes(cohort, expr, fill = Group, color = Group)) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA, width = 0.6,
               position = position_dodge(0.7)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
             size = 0.6, alpha = 0.5) +
  geom_text(data = pvals, aes(cohort, ymax * 1.08, label = label), inherit.aes = FALSE,
            size = 3, color = "grey30") +
  scale_fill_manual(values = c(Normal = COL_NORMAL, Tumor = COL_TUMOR)) +
  scale_color_manual(values = c(Normal = COL_NORMAL, Tumor = COL_TUMOR)) +
  facet_wrap(~ piRNA, scales = "free_y") +
  labs(x = NULL, y = "log2(TPM+1)",
       title = "Expression of the three piRNAs, Tumor vs Normal (Wilcoxon)") +
  theme_pub + theme(axis.text.x = element_text(angle = 20, hjust = 1),
                    legend.position = "top")
save_plot(p_box, "Fig_expression_boxplots", 15, 5.5)

# ==============================================================================
# 7. PLOT 8 — real ROC (5-fold CV logistic regression on the 3 piRNAs)
# ==============================================================================
cat("\n========== STEP 7: ROC (5-fold CV) ==========\n")

cv_prob <- function(X, y, k = 5) {
  set.seed(SEED)
  folds <- sample(rep(seq_len(k), length.out = length(y)))
  prob <- numeric(length(y))
  for (f in seq_len(k)) {
    tr <- folds != f
    w  <- ifelse(y[tr] == "Tumor", 0.5 / mean(y[tr] == "Tumor"),
                 0.5 / mean(y[tr] == "Normal"))      # balance classes
    d  <- data.frame(y = as.integer(y[tr] == "Tumor"), X[tr, , drop = FALSE],
                     check.names = FALSE)
    colnames(d) <- make.names(colnames(d))
    fit <- suppressWarnings(glm(y ~ ., data = d, family = binomial, weights = w))
    nd  <- as.data.frame(X[!tr, , drop = FALSE]); colnames(nd) <- make.names(colnames(nd))
    prob[!tr] <- predict(fit, nd, type = "response")
  }
  prob
}

roc_rows <- list(); roc_objs <- list()
for (cn in names(cohorts)) {
  co  <- cohorts[[cn]]
  pr  <- cv_prob(co$X, co$y)
  ro  <- roc(co$y, pr, levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE)
  ci  <- as.numeric(ci.auc(ro, method = "delong"))
  roc_objs[[sprintf("%s: AUC %.2f (%.2f-%.2f)", cn, ci[2], ci[1], ci[3])]] <- ro
  roc_rows[[cn]] <- data.frame(cohort = cn, n = length(co$y), n_tumor = sum(co$y == "Tumor"),
                               AUC = ci[2], CI_low = ci[1], CI_high = ci[3])
}
roc_tab <- do.call(rbind, roc_rows)
print(roc_tab, row.names = FALSE)
write.csv(roc_tab, file.path(OUT, "ROC_3piRNA_real.csv"), row.names = FALSE)

p_roc <- ggroc(roc_objs, legacy.axes = TRUE, linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = COHORT_COL, name = NULL) +
  coord_equal() +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = "3-piRNA model ROC (5-fold CV, DeLong 95% CI)") +
  theme_pub + theme(legend.position = c(0.62, 0.15),
                    legend.background = element_rect(fill = "white", color = NA))
save_plot(p_roc, "Fig_ROC_3piRNA_real", 7, 7)

cat(sprintf("\n=== 07 complete. Runtime: %.1f min ===\n",
            as.numeric(difftime(Sys.time(), start_time, units = "mins"))))
cat("=================================================================\n")

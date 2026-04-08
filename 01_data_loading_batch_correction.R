################################################################################
#                                                                              #
#   01 — Data Loading & Batch Correction (PHASE 0)                            #
#                                                                              #
#   Load all cohort CSVs, harmonise labels, find shared piRNAs,               #
#   log2(TPM+1), ComBat batch correction, Z-score normalisation.              #
#   Output a single analysis-ready data frame.                                #
#                                                                              #
#   Training:  BRCA1, PRJNA294226, PRJNA482141, PRJNA808405, PRJNA934049      #
#   Validation (LOCKED): yyfbatch1 (tissue), yyfbatch2 (plasma/serum)         #
#                                                                              #
#   Run FIRST — all subsequent scripts depend on outputs from this one.       #
#                                                                              #
################################################################################

start_time <- Sys.time()
cat("=================================================================\n")
cat("  PHASE 0: Data Loading & Batch Correction\n")
cat("=================================================================\n\n")

# ==============================================================================
# 0. PACKAGES
# ==============================================================================
required_pkgs <- c("sva", "limma", "ggplot2", "dplyr", "cowplot")
bioc_pkgs     <- c("sva", "limma")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}
for (pkg in setdiff(required_pkgs, bioc_pkgs)) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(sva)
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})
cat("All packages loaded.\n\n")

# ==============================================================================
# 0.1 GLOBAL SETTINGS
# ==============================================================================
SEED <- 2024
set.seed(SEED)

COLOR_PALETTE <- c("#E31A1C", "#1F78B4", "#33A02C", "#FF7F00",
                    "#6A3D9A", "#B15928", "#A6CEE3")

my_theme <- theme_minimal(base_size = 13) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

dir.create("results/qc",               recursive = TRUE, showWarnings = FALSE)
dir.create("results/feature_selection", recursive = TRUE, showWarnings = FALSE)
dir.create("results/models",            recursive = TRUE, showWarnings = FALSE)
dir.create("results/validation",        recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. DATA LOADING
# ==============================================================================
cat("========== STEP 1: Loading Cohort Files ==========\n")

recode_group <- function(df, group_col = "Group") {
  g <- as.character(df[[group_col]])
  g[g %in% c("Benign", "Cancer")] <- ifelse(
    g[g %in% c("Benign", "Cancer")] == "Cancer", "Tumor", "Normal")
  g[g %in% c("cancer", "tumor")] <- "Tumor"
  g[g %in% c("normal", "benign", "control", "healthy", "Control", "Healthy")] <- "Normal"
  df[[group_col]] <- factor(g, levels = c("Normal", "Tumor"))
  df[df[[group_col]] %in% c("Normal", "Tumor"), ]
}

dedup_tcga_samples <- function(df) {
  ids <- rownames(df)
  is_tcga <- grepl("^TCGA", ids, ignore.case = TRUE)
  if (!any(is_tcga)) return(df)
  keep <- !is_tcga | grepl("_01A$", ids) | grepl("_11A$", ids)
  n_dropped <- sum(!keep)
  if (n_dropped > 0) cat(sprintf("    Dropped %d non-primary TCGA vials\n", n_dropped))
  df <- df[keep, , drop = FALSE]
  dup <- duplicated(rownames(df))
  if (any(dup)) {
    cat(sprintf("    Removed %d duplicate row IDs\n", sum(dup)))
    df <- df[!dup, , drop = FALSE]
  }
  df
}

dataset_names <- c("BRCA1", "PRJNA294226", "PRJNA482141",
                    "PRJNA808405", "PRJNA934049",
                    "yyfbatch1", "yyfbatch2")

training_cohorts    <- c("BRCA1", "PRJNA294226", "PRJNA482141",
                          "PRJNA808405", "PRJNA934049")
validation_cohorts  <- c("yyfbatch1", "yyfbatch2")

ready_list <- list()
for (nm in dataset_names) {
  # Try multiple file name patterns
  candidates <- c(
    file.path("processed_results", paste0(nm, "_processed.csv")),
    file.path("processed_results", paste0(nm, "_processed_1.csv"))
  )
  fpath <- NULL
  for (cand in candidates) {
    if (file.exists(cand)) { fpath <- cand; break }
  }
  if (is.null(fpath)) {
    cat(sprintf("  WARNING: %s not found — skipping.\n", nm))
    next
  }
  tmp <- read.csv(fpath, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(tmp) <- make.unique(as.character(tmp[[1]]))
  tmp[[1]] <- NULL
  tmp <- dedup_tcga_samples(tmp)
  tmp <- recode_group(tmp)
  ready_list[[nm]] <- tmp
  tb <- table(tmp$Group)
  cat(sprintf("  Loaded %-15s: %d samples (Normal=%d, Tumor=%d)\n",
              nm, nrow(tmp), tb["Normal"], tb["Tumor"]))
}

if (length(ready_list) == 0) stop("No cohort files found! Check processed_results/ directory.")

# ==============================================================================
# 2. BRCA1 BALANCING (matched pairs + 40 % excess tumors)
# ==============================================================================
cat("\n========== STEP 2: BRCA1 Class Balancing ==========\n")

if ("BRCA1" %in% names(ready_list)) {
  set.seed(SEED)
  brca <- ready_list[["BRCA1"]]
  idx_normal <- which(brca$Group == "Normal")
  idx_tumor  <- which(brca$Group == "Tumor")
  n_matched  <- min(length(idx_normal), length(idx_tumor))
  n_remaining <- length(idx_tumor) - n_matched
  n_extra     <- round(n_remaining * 0.40)
  idx_tumor_keep <- sample(idx_tumor, n_matched + n_extra)
  brca_balanced <- brca[c(idx_normal, idx_tumor_keep), ]
  ready_list[["BRCA1"]] <- brca_balanced
  cat(sprintf("  BRCA1 balanced: Normal=%d, Tumor=%d (matched %d + 40%% of %d excess = %d)\n",
              length(idx_normal), n_matched + n_extra, n_matched, n_remaining, n_extra))
} else {
  cat("  BRCA1 not found — skipping balancing.\n")
}

# ==============================================================================
# 3. COMMON piRNAs + LOG2(TPM+1) + GLOBAL ComBat + Z-SCORE
# ==============================================================================
cat("\n========== STEP 3: Batch Correction (Global ComBat) ==========\n")

# 3a. Find common piRNAs
all_gene_sets <- lapply(ready_list, function(df) setdiff(colnames(df), "Group"))
common_pirnas <- Reduce(intersect, all_gene_sets)
cat("Common piRNAs across all", length(ready_list), "cohorts:", length(common_pirnas), "\n")

if (length(common_pirnas) < 10)
  stop("Too few common piRNAs (", length(common_pirnas), "). Check input data.")

# 3b. Subset to common piRNAs and log2(TPM+1)
clean_list <- lapply(ready_list, function(df) {
  df_sub <- df[, c("Group", common_pirnas)]
  mat <- as.matrix(df_sub[, -1])
  mat[mat < 0] <- 0
  if (max(mat, na.rm = TRUE) > 50) mat <- log2(mat + 1)
  df_sub[, -1] <- as.data.frame(mat)
  df_sub
})

# 3c. Combine expression matrices
expr_matrices <- lapply(clean_list, function(df) t(as.matrix(df[, -1])))
combined_expr <- do.call(cbind, expr_matrices)

batch_vec <- unlist(lapply(names(clean_list), function(n) rep(n, nrow(clean_list[[n]]))))
group_vec <- unlist(lapply(clean_list, function(df) as.character(df$Group)))

# 3d. PCA BEFORE ComBat
cat("\nComputing PCA before ComBat...\n")
pca_before <- prcomp(t(combined_expr), center = TRUE, scale. = TRUE)
pca_df_before <- data.frame(
  PC1    = pca_before$x[, 1],
  PC2    = pca_before$x[, 2],
  Cohort = batch_vec,
  Group  = group_vec
)
var_before <- summary(pca_before)$importance[2, 1:2] * 100

p_before <- ggplot(pca_df_before, aes(PC1, PC2, color = Cohort, shape = Group)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_manual(values = COLOR_PALETTE[seq_len(length(unique(batch_vec)))]) +
  labs(title = "PCA Before ComBat",
       x = sprintf("PC1 (%.1f%%)", var_before[1]),
       y = sprintf("PC2 (%.1f%%)", var_before[2])) +
  my_theme

# 3e. ComBat batch correction
cat("Running ComBat on", ncol(combined_expr), "samples across",
    length(unique(batch_vec)), "batches...\n")
mod <- model.matrix(~ as.factor(group_vec))

if (length(unique(batch_vec)) >= 2) {
  combat_expr <- ComBat(dat = combined_expr, batch = batch_vec,
                        mod = mod, par.prior = TRUE, prior.plots = FALSE)
  cat("ComBat correction completed.\n")
} else {
  cat("WARNING: Only 1 batch — skipping ComBat.\n")
  combat_expr <- combined_expr
}

# 3f. Z-score normalisation
combat_df_all <- data.frame(
  Group = factor(group_vec, levels = c("Normal", "Tumor")),
  t(combat_expr),
  check.names = FALSE
)
combat_df_all$Batch <- batch_vec
rownames(combat_df_all) <- colnames(combat_expr)

gene_cols <- setdiff(colnames(combat_df_all), c("Group", "Batch"))

# Save pre-Z-score parameters for validation preprocessing
z_means <- colMeans(combat_df_all[, gene_cols])
z_sds   <- apply(combat_df_all[, gene_cols], 2, sd)
z_sds[z_sds == 0] <- 1

combat_df_all[, gene_cols] <- scale(combat_df_all[, gene_cols])
combat_df_all[is.na(combat_df_all)] <- 0

preprocessing_params <- list(z_means = z_means, z_sds = z_sds, gene_cols = gene_cols)
saveRDS(preprocessing_params, "results/models/preprocessing_params.rds")

cat("ComBat + Z-score complete. Total:", nrow(combat_df_all), "samples x",
    length(gene_cols), "piRNAs\n")

# 3g. PCA AFTER ComBat
cat("\nComputing PCA after ComBat...\n")
pca_after <- prcomp(as.matrix(combat_df_all[, gene_cols]), center = TRUE, scale. = TRUE)
pca_df_after <- data.frame(
  PC1    = pca_after$x[, 1],
  PC2    = pca_after$x[, 2],
  Cohort = combat_df_all$Batch,
  Group  = combat_df_all$Group
)
var_after <- summary(pca_after)$importance[2, 1:2] * 100

p_after <- ggplot(pca_df_after, aes(PC1, PC2, color = Cohort, shape = Group)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_manual(values = COLOR_PALETTE[seq_len(length(unique(pca_df_after$Cohort)))]) +
  labs(title = "PCA After ComBat",
       x = sprintf("PC1 (%.1f%%)", var_after[1]),
       y = sprintf("PC2 (%.1f%%)", var_after[2])) +
  my_theme

p_combined <- plot_grid(p_before, p_after, ncol = 2, labels = c("A", "B"))
ggsave("results/qc/pca_combat.png", p_combined, width = 14, height = 6, dpi = 300, bg = "white")
cat("  Saved: results/qc/pca_combat.png\n")

# ==============================================================================
# 4. DATA SPLITTING
# ==============================================================================
cat("\n========== STEP 4: Data Splitting ==========\n")

cat("\nPost-ComBat class distribution by batch:\n")
print(table(combat_df_all$Batch, combat_df_all$Group))

available_training   <- intersect(training_cohorts, unique(combat_df_all$Batch))
available_validation <- intersect(validation_cohorts, unique(combat_df_all$Batch))

train_df         <- combat_df_all[combat_df_all$Batch %in% available_training, ]
valid_yyfbatch1  <- combat_df_all[combat_df_all$Batch == "yyfbatch1", ]
valid_yyfbatch2  <- combat_df_all[combat_df_all$Batch == "yyfbatch2", ]

cat("\n  Training cohorts:", paste(available_training, collapse = ", "), "\n")
cat("  Training samples:", nrow(train_df), "\n")
cat("  yyfbatch1 (tissue validation):", nrow(valid_yyfbatch1), "samples\n")
cat("  yyfbatch2 (plasma/serum validation):", nrow(valid_yyfbatch2), "samples\n")

# ==============================================================================
# 5. SAVE OUTPUTS
# ==============================================================================
cat("\n========== STEP 5: Saving Outputs ==========\n")

saveRDS(combat_df_all,    "results/models/combat_df_all.rds")
saveRDS(train_df,         "results/models/train_df.rds")
saveRDS(valid_yyfbatch1,  "results/models/valid_yyfbatch1.rds")
saveRDS(valid_yyfbatch2,  "results/models/valid_yyfbatch2.rds")
saveRDS(gene_cols,        "results/models/common_pirnas.rds")

cat("  Saved: results/models/combat_df_all.rds\n")
cat("  Saved: results/models/train_df.rds\n")
cat("  Saved: results/models/valid_yyfbatch1.rds\n")
cat("  Saved: results/models/valid_yyfbatch2.rds\n")
cat("  Saved: results/models/common_pirnas.rds\n")
cat("  Saved: results/models/preprocessing_params.rds\n")

end_time <- Sys.time()
cat(sprintf("\n=== PHASE 0 complete. Runtime: %.1f minutes ===\n",
            as.numeric(difftime(end_time, start_time, units = "mins"))))
cat("=================================================================\n")
cat("  Next: run 02_feature_selection.R\n")
cat("=================================================================\n")

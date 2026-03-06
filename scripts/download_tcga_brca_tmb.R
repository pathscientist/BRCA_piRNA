################################################################################
#                                                                              #
#   Download & Compute TCGA-BRCA TMB (Tumor Mutation Burden) Data             #
#                                                                              #
#   This script downloads somatic mutation data (MAF) from TCGA-BRCA via      #
#   TCGAbiolinks, computes TMB (nonsynonymous mutations per megabase),        #
#   and outputs tmb_data.csv for use in piRNA_advanced_analysis.R.            #
#                                                                              #
#   Output: tmb_data.csv with columns:                                         #
#     - SampleID        : TCGA barcode (first 15 chars, e.g. TCGA-XX-XXXX-01) #
#     - Tumor_Sample_Barcode : Full barcode from MAF                          #
#     - TMB             : Mutations per megabase (nonsynonymous / 38 Mb)      #
#     - Total_Mutations : Total nonsynonymous mutation count                   #
#     - TMB_Category    : "High" (>=10 mut/Mb) or "Low" (<10 mut/Mb)         #
#                                                                              #
#   Usage: Rscript scripts/download_tcga_brca_tmb.R                           #
#                                                                              #
################################################################################

cat("=== TCGA-BRCA TMB Data Download & Computation ===\n\n")

# ==============================================================================
# 1. INSTALL & LOAD PACKAGES
# ==============================================================================

cran_pkgs <- c("dplyr", "data.table")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, dependencies = TRUE, quiet = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", quiet = TRUE)

bioc_pkgs <- c("TCGAbiolinks", "maftools")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(maftools)
  library(dplyr)
  library(data.table)
})

cat("Packages loaded.\n\n")

# ==============================================================================
# 2. DOWNLOAD TCGA-BRCA SOMATIC MUTATIONS (MAF)
# ==============================================================================

cat("Step 1: Querying GDC for TCGA-BRCA somatic mutations...\n")

# Query masked somatic mutations (open-access, MuTect2 caller)
query_maf <- GDCquery(
  project    = "TCGA-BRCA",
  data.category = "Simple Nucleotide Variation",
  data.type     = "Masked Somatic Mutation",
  access        = "open"
)

cat("Step 2: Downloading MAF files from GDC...\n")
GDCdownload(query_maf)

cat("Step 3: Preparing MAF data...\n")
maf_data <- GDCprepare(query_maf)

cat(sprintf("  Downloaded %d mutation records.\n", nrow(maf_data)))
cat(sprintf("  Unique tumor samples: %d\n",
            length(unique(maf_data$Tumor_Sample_Barcode))))

# ==============================================================================
# 3. COMPUTE TMB (Tumor Mutation Burden)
# ==============================================================================

cat("\nStep 4: Computing TMB...\n")

# Define nonsynonymous variant classifications
# These are the protein-altering mutations used for standard TMB calculation
nonsynonymous_classes <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "In_Frame_Del",
  "In_Frame_Ins",
  "Splice_Site",
  "Nonstop_Mutation",
  "Translation_Start_Site"
)

# Filter to nonsynonymous mutations only
maf_nonsyn <- maf_data[maf_data$Variant_Classification %in% nonsynonymous_classes, ]

cat(sprintf("  Total mutations: %d\n", nrow(maf_data)))
cat(sprintf("  Nonsynonymous mutations: %d\n", nrow(maf_nonsyn)))

# Count nonsynonymous mutations per sample
mut_counts <- maf_nonsyn %>%
  group_by(Tumor_Sample_Barcode) %>%
  summarise(Total_Mutations = n(), .groups = "drop")

# WES exome capture size: ~38 Mb (standard for TCGA WES)
EXOME_SIZE_MB <- 38

# Compute TMB = nonsynonymous mutations / exome size
mut_counts$TMB <- mut_counts$Total_Mutations / EXOME_SIZE_MB

# Create short SampleID (first 15 characters = patient-level + sample type)
# TCGA barcode: TCGA-XX-XXXX-01A-11D-A123-09
#   Chars 1-12: Patient barcode (TCGA-XX-XXXX)
#   Chars 13-15: Sample type (01A = primary tumor)
mut_counts$SampleID <- substr(mut_counts$Tumor_Sample_Barcode, 1, 15)

# Categorize TMB
# Standard threshold: 10 mutations/Mb (used in clinical practice, e.g. FDA)
mut_counts$TMB_Category <- ifelse(mut_counts$TMB >= 10, "High", "Low")

# Also add a research-oriented tertile classification
tmb_tertiles <- quantile(mut_counts$TMB, probs = c(1/3, 2/3), na.rm = TRUE)
mut_counts$TMB_Tertile <- cut(
  mut_counts$TMB,
  breaks = c(-Inf, tmb_tertiles[1], tmb_tertiles[2], Inf),
  labels = c("Low", "Medium", "High")
)

# Sort by TMB descending
mut_counts <- mut_counts %>% arrange(desc(TMB))

# ==============================================================================
# 4. SUMMARY STATISTICS
# ==============================================================================

cat("\n--- TMB Summary Statistics ---\n")
cat(sprintf("  Samples with TMB data: %d\n", nrow(mut_counts)))
cat(sprintf("  Median TMB: %.2f mutations/Mb\n", median(mut_counts$TMB)))
cat(sprintf("  Mean TMB:   %.2f mutations/Mb\n", mean(mut_counts$TMB)))
cat(sprintf("  Range:      %.2f - %.2f mutations/Mb\n",
            min(mut_counts$TMB), max(mut_counts$TMB)))
cat(sprintf("  TMB-High (>=10 mut/Mb): %d (%.1f%%)\n",
            sum(mut_counts$TMB_Category == "High"),
            100 * mean(mut_counts$TMB_Category == "High")))
cat(sprintf("  TMB-Low  (<10 mut/Mb):  %d (%.1f%%)\n",
            sum(mut_counts$TMB_Category == "Low"),
            100 * mean(mut_counts$TMB_Category == "Low")))

cat("\n  TMB Tertile distribution:\n")
print(table(mut_counts$TMB_Tertile))

# Mutation type breakdown
cat("\n  Mutation type breakdown (top 10):\n")
var_class_counts <- sort(table(maf_nonsyn$Variant_Classification), decreasing = TRUE)
print(head(var_class_counts, 10))

# ==============================================================================
# 5. SAVE OUTPUT
# ==============================================================================

# Select and order output columns
tmb_output <- mut_counts %>%
  select(SampleID, Tumor_Sample_Barcode, TMB, Total_Mutations,
         TMB_Category, TMB_Tertile)

output_path <- "tmb_data.csv"
write.csv(tmb_output, output_path, row.names = FALSE)
cat(sprintf("\nTMB data saved to: %s\n", output_path))
cat(sprintf("  Rows: %d, Columns: %d\n", nrow(tmb_output), ncol(tmb_output)))

# Also save a quick-reference version with just SampleID and TMB
# using the 12-char patient barcode for easier matching
tmb_simple <- tmb_output %>%
  mutate(PatientID = substr(SampleID, 1, 12)) %>%
  select(PatientID, SampleID, TMB, TMB_Category)

write.csv(tmb_simple, "tmb_data_simple.csv", row.names = FALSE)
cat("Simple version saved to: tmb_data_simple.csv\n")

# ==============================================================================
# 6. OPTIONAL: READ MAF WITH maftools FOR ADDITIONAL SUMMARIES
# ==============================================================================

cat("\nStep 5: Generating maftools summary...\n")

maf_obj <- tryCatch({
  read.maf(maf = maf_data, isTCGA = TRUE)
}, error = function(e) {
  cat("  maftools summary skipped:", conditionMessage(e), "\n")
  NULL
})

if (!is.null(maf_obj)) {
  # Save oncoplot for top mutated genes (informational)
  dir.create("results/tmb", recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    png("results/tmb/tcga_brca_oncoplot_top20.png",
        width = 14, height = 8, units = "in", res = 300)
    oncoplot(maf = maf_obj, top = 20,
             titleText = "TCGA-BRCA: Top 20 Mutated Genes")
    dev.off()
    cat("  Oncoplot saved to results/tmb/tcga_brca_oncoplot_top20.png\n")
  }, error = function(e) {
    cat("  Oncoplot generation skipped:", conditionMessage(e), "\n")
    try(dev.off(), silent = TRUE)
  })

  # Save TMB plot from maftools
  tryCatch({
    png("results/tmb/tcga_brca_tmb_distribution.png",
        width = 10, height = 6, units = "in", res = 300)
    tcgaCompare(maf = maf_obj, cohortName = "TCGA-BRCA")
    dev.off()
    cat("  TMB comparison plot saved to results/tmb/tcga_brca_tmb_distribution.png\n")
  }, error = function(e) {
    cat("  TMB comparison plot skipped:", conditionMessage(e), "\n")
    try(dev.off(), silent = TRUE)
  })

  # Gene summary
  gene_summary <- getSampleSummary(maf_obj)
  write.csv(gene_summary, "results/tmb/tcga_brca_sample_mutation_summary.csv",
            row.names = FALSE)
  cat("  Sample mutation summary saved.\n")
}

cat("\n=== TMB DATA DOWNLOAD COMPLETE ===\n")
cat("\nNext steps:\n")
cat("  1. Run piRNA_advanced_analysis.R — it will auto-detect tmb_data.csv\n")
cat("  2. TMB will be matched to your TCGA-BRCA samples by barcode\n")
cat("  3. Violin plots, stage analysis, and T-Score correlation will use real data\n")

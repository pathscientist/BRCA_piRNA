setwd("~/GEO_BRCA/meta/brca_dataset_1101tpm")

rm(list = ls())

install.packages("GGally")
install.packages("ggfortify")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")


# Load required libraries
library(limma)
library(sva)

library(GGally)  # Make sure 'GGally' package is installed and loaded
library(ggbiplot)

BRCA1 <- read.table("BRCA1.txt", header = TRUE, row.names = 1, sep = "\t")
#plasma <- read.table("plasma.txt", header = TRUE, row.names = 1, sep = "\t")
PRJNA294226 <- read.table("PRJNA294226.txt", header = TRUE, row.names = 1, sep = "\t")
PRJNA482141 <- read.table("PRJNA482141.txt", header = TRUE, row.names = 1, sep = "\t")
PRJNA808405 <- read.table("PRJNA808405.txt", header = TRUE, row.names = 1, sep = "\t")
PRJNA934049 <- read.table("PRJNA934049.txt", header = TRUE, row.names = 1, sep = "\t")
#tissue <- read.table("tissue.txt", header = TRUE, row.names = 1, sep = "\t")
yyfbatch1 <- read.table("yyfbatch1.txt", header = TRUE, row.names = 1, sep = "\t")
yyfbatch2 <- read.table("yyfbatch2.txt", header = TRUE, row.names = 1, sep = "\t")

# Combine datasets and create batch indicator
combined_data <- cbind(BRCA1, PRJNA294226, PRJNA482141
                       , PRJNA808405, PRJNA934049, yyfbatch1, yyfbatch2)

tpm_data <- combined_data
# Assuming `tpm_data` is your DataFrame with TPM values
log_tpm_data <- log2(tpm_data + 1)

batch_indicator <- factor(c(rep("BRCA1", ncol(BRCA1)), rep("PRJNA294226", ncol(PRJNA294226))
                            , rep("PRJNA482141", ncol(PRJNA482141)), rep("PRJNA808405", ncol(PRJNA808405)), rep("PRJNA934049", ncol(PRJNA934049))
                            , rep("yyfbatch1", ncol(yyfbatch1)), rep("yyfbatch2", ncol(yyfbatch2))))

exp_m_t <- as.data.frame(t(log_tpm_data))
design <- read.table("group.txt", header = TRUE, sep = "\t")

adjusted_data <- ComBat(dat = log_tpm_data, batch = batch_indicator, mod = NULL, par.prior = TRUE)



modtype <- design$Type
batchtype <- design$Group
mod <- model.matrix(~as.factor(modtype))

adjusted_data2 <- ComBat(dat = log_tpm_data, batch = batch_indicator, mod = mod, par.prior = TRUE)


# Simple PCA Plot before and after adjustment
library(ggplot2)

# PCA on original data
pca_original <- prcomp(t(log_tpm_data))
original_df <- data.frame(PC1 = pca_original$x[,1], PC2 = pca_original$x[,2], Batch = batch_indicator)

# PCA on adjusted data
pca_adjusted <- prcomp(t(adjusted_data))
adjusted_df <- data.frame(PC1 = pca_adjusted$x[,1], PC2 = pca_adjusted$x[,2], Batch = batch_indicator)

# Plotting
ggplot(original_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point() +
  ggtitle("Original Data") +
  theme_classic()

ggplot(adjusted_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point() +
  ggtitle("Adjusted Data") +
  theme_classic()



# PCA on original data
pca_original <- prcomp(t(log_tpm_data))
original_df <- data.frame(PC1 = pca_original$x[,1], PC2 = pca_original$x[,2], Batch = modtype)

# PCA on adjusted data
pca_adjusted <- prcomp(t(adjusted_data2))
adjusted_df <- data.frame(PC1 = pca_adjusted$x[,1], PC2 = pca_adjusted$x[,2], Batch = modtype)

# Plotting
ggplot(original_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point() +
  stat_ellipse() +
  ggtitle("Original Data") +
  theme_classic()

ggplot(adjusted_df, aes(x = PC1, y = PC2, color = Batch)) +
  geom_point() +
  stat_ellipse() +
  ggtitle("Adjusted Data") +
  theme_classic()

# Write the normalized data to a text file
write.table(adjusted_data2, file = "normalized_combat_20240506.txt", sep = "\t", quote = FALSE, col.names = FALSE)
gene_exp_combat = adjusted_data2
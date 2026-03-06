#!/bin/bash
#SBATCH --job-name=piRNA_mRNA_corr
#SBATCH --output=piRNA_mRNA_corr_%j.out
#SBATCH --error=piRNA_mRNA_corr_%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=standard

# ===========================================================
# piRNA-mRNA Spearman correlation analysis - SLURM submission
#
# Usage:
#   sbatch piRNA_mRNA_correlation.sh
#
# Adjust --cpus-per-task, --mem, --time, --partition as needed
# for your HPC cluster.
# ===========================================================

echo "Job started: $(date)"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"

# Load modules (uncomment/edit for your HPC)
# module load python/3.10
# module load anaconda3

# Activate conda env if needed
# conda activate brca_pirna

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"

python3 "${SCRIPT_DIR}/piRNA_mRNA_correlation.py" \
    --mrna "${PROJECT_DIR}/processed_results/TCGA_BRCA_tpm_all_samples.csv" \
    --pirna "${PROJECT_DIR}/processed_results/BRCA1_processed.csv" \
    --outdir "${PROJECT_DIR}/results" \
    --workers "${SLURM_CPUS_PER_TASK:-16}" \
    --fdr-cutoff 0.05

echo "Job finished: $(date)"

#!/usr/bin/env python3
"""
Calculate Spearman correlations between piRNA and mRNA expression
using TCGA BRCA data.

Usage:
  # Single node (uses all cores):
  python3 piRNA_mRNA_correlation.py

  # HPC with SLURM (submit via the companion .sh script):
  sbatch piRNA_mRNA_correlation.sh

Input files (place in same directory or edit paths below):
  - TCGA_BRCA_tpm_all_samples.csv   (mRNA TPM, tab-separated)
  - BRCA1_processed.csv             (piRNA expression, tab-separated)

Output:
  - piRNA_mRNA_correlations.csv       (all piRNA-mRNA pairs with rho & p-value)
  - piRNA_mRNA_correlations_sig.csv   (FDR < 0.05 pairs only)
"""

import os
import sys
import time
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from multiprocessing import Pool, cpu_count

# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------
INPUT_DIR = os.path.dirname(os.path.abspath(__file__))
MRNA_FILE = os.path.join(INPUT_DIR, "..", "processed_results", "TCGA_BRCA_tpm_all_samples.csv")
PIRNA_FILE = os.path.join(INPUT_DIR, "..", "processed_results", "BRCA1_processed.csv")
OUTPUT_DIR = os.path.join(INPUT_DIR, "..", "results")


def load_data(mrna_path, pirna_path):
    """Load and align mRNA and piRNA expression matrices on common samples."""
    print(f"Loading mRNA data from: {mrna_path}")
    mrna = pd.read_csv(mrna_path, sep="\t", index_col=0)
    print(f"  mRNA matrix: {mrna.shape[0]} samples x {mrna.shape[1]} genes")

    print(f"Loading piRNA data from: {pirna_path}")
    pirna = pd.read_csv(pirna_path, sep="\t", index_col=0)
    # Drop the Group column if present
    if "Group" in pirna.columns:
        pirna = pirna.drop(columns=["Group"])
    print(f"  piRNA matrix: {pirna.shape[0]} samples x {pirna.shape[1]} piRNAs")

    # Find common samples
    common = mrna.index.intersection(pirna.index)
    print(f"  Common samples: {len(common)}")
    if len(common) == 0:
        sys.exit("ERROR: No common samples found. Check sample ID formats.")

    mrna = mrna.loc[common].astype(np.float64)
    pirna = pirna.loc[common].astype(np.float64)

    # Filter: remove features with zero variance (constant across all samples)
    mrna_var = mrna.var()
    pirna_var = pirna.var()
    mrna = mrna.loc[:, mrna_var > 0]
    pirna = pirna.loc[:, pirna_var > 0]
    print(f"  After variance filter: {mrna.shape[1]} genes, {pirna.shape[1]} piRNAs")

    return mrna, pirna


def _correlate_chunk(args):
    """Worker function: compute Spearman correlation for a chunk of piRNAs."""
    pirna_names, pirna_mat, mrna_mat, mrna_names = args
    results = []
    for i, pi_name in enumerate(pirna_names):
        pi_vals = pirna_mat[i, :]
        for j, gene_name in enumerate(mrna_names):
            rho, pval = stats.spearmanr(pi_vals, mrna_mat[j, :])
            results.append((pi_name, gene_name, rho, pval))
    return results


def run_correlations(pirna, mrna, n_workers=None):
    """Compute all piRNA-mRNA Spearman correlations in parallel."""
    if n_workers is None:
        n_workers = min(cpu_count(), 32)

    pirna_names = pirna.columns.tolist()
    mrna_names = mrna.columns.tolist()
    n_pirna = len(pirna_names)
    total_pairs = n_pirna * len(mrna_names)
    print(f"\nComputing {total_pairs:,} piRNA-mRNA correlations "
          f"({n_pirna} piRNAs x {len(mrna_names)} genes) "
          f"using {n_workers} workers...")

    # Transpose to (features x samples) numpy arrays for fast access
    pirna_mat = pirna.values.T  # (n_pirna, n_samples)
    mrna_mat = mrna.values.T    # (n_genes, n_samples)

    # Split piRNAs into chunks for parallel processing
    chunk_size = max(1, n_pirna // n_workers)
    chunks = []
    for start in range(0, n_pirna, chunk_size):
        end = min(start + chunk_size, n_pirna)
        chunks.append((
            pirna_names[start:end],
            pirna_mat[start:end, :],
            mrna_mat,
            mrna_names,
        ))

    t0 = time.time()
    all_results = []
    with Pool(n_workers) as pool:
        for i, chunk_results in enumerate(pool.imap_unordered(_correlate_chunk, chunks)):
            all_results.extend(chunk_results)
            elapsed = time.time() - t0
            done = len(all_results)
            print(f"  Progress: {done:,}/{total_pairs:,} pairs "
                  f"({100*done/total_pairs:.1f}%) [{elapsed:.0f}s]",
                  flush=True)

    elapsed = time.time() - t0
    print(f"  Completed in {elapsed:.1f}s")

    df = pd.DataFrame(all_results,
                      columns=["piRNA", "Gene", "Spearman_rho", "p_value"])
    return df


def adjust_pvalues(df):
    """Apply Benjamini-Hochberg FDR correction."""
    from statsmodels.stats.multitest import multipletests
    print("Applying Benjamini-Hochberg FDR correction...")
    _, fdr, _, _ = multipletests(df["p_value"].values, method="fdr_bh")
    df["FDR"] = fdr
    return df


def main():
    parser = argparse.ArgumentParser(
        description="piRNA-mRNA Spearman correlation analysis")
    parser.add_argument("--mrna", default=MRNA_FILE,
                        help="Path to mRNA TPM file (tab-separated)")
    parser.add_argument("--pirna", default=PIRNA_FILE,
                        help="Path to piRNA expression file (tab-separated)")
    parser.add_argument("--outdir", default=OUTPUT_DIR,
                        help="Output directory")
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel workers (default: all CPUs)")
    parser.add_argument("--fdr-cutoff", type=float, default=0.05,
                        help="FDR threshold for significant pairs")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load and align data
    mrna, pirna = load_data(args.mrna, args.pirna)

    # Run correlations
    results = run_correlations(pirna, mrna, n_workers=args.workers)

    # FDR correction
    results = adjust_pvalues(results)

    # Sort by absolute correlation
    results = results.sort_values("Spearman_rho", key=abs, ascending=False)

    # Save all results
    out_all = os.path.join(args.outdir, "piRNA_mRNA_correlations.csv")
    results.to_csv(out_all, index=False)
    print(f"\nAll results saved: {out_all} ({len(results):,} pairs)")

    # Save significant results
    sig = results[results["FDR"] < args.fdr_cutoff]
    out_sig = os.path.join(args.outdir, "piRNA_mRNA_correlations_sig.csv")
    sig.to_csv(out_sig, index=False)
    print(f"Significant results (FDR < {args.fdr_cutoff}): "
          f"{out_sig} ({len(sig):,} pairs)")

    # Summary
    print(f"\n{'='*60}")
    print(f"Summary:")
    print(f"  Total pairs tested: {len(results):,}")
    print(f"  Significant (FDR < {args.fdr_cutoff}): {len(sig):,}")
    if len(sig) > 0:
        pos = sig[sig["Spearman_rho"] > 0]
        neg = sig[sig["Spearman_rho"] < 0]
        print(f"    Positive correlations: {len(pos):,}")
        print(f"    Negative correlations: {len(neg):,}")
        print(f"  Top 10 by |rho|:")
        for _, row in sig.head(10).iterrows():
            print(f"    {row['piRNA']} ~ {row['Gene']}: "
                  f"rho={row['Spearman_rho']:.4f}, FDR={row['FDR']:.2e}")


if __name__ == "__main__":
    main()

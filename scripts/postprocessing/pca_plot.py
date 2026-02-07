#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import os


def perform_pca_plot(df_path, output_fig="PCA_Tanimoto_Publication.png",
                     output_csv="pca_results.csv"):
    """
    Perform PCA on similarity matrix and create scatter plot.

    Args:
        df_path (str): Path to CSV file containing similarity matrix
        output_fig (str): Path to save output figure
        output_csv (str): Path to save PCA results CSV
    """
    # Load similarity matrix
    if not os.path.exists(df_path):
        raise FileNotFoundError(f"Similarity matrix file not found: {df_path}")

    df = pd.read_csv(df_path, index_col=0)
    print(f"Loaded similarity matrix: {df.shape[0]} x {df.shape[1]}")

    # Standardize and perform PCA
    scaled_data = StandardScaler().fit_transform(df.values)
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(scaled_data)
    pca_df = pd.DataFrame(principal_components, columns=["PC1", "PC2"], index=df.index)

    # Plot (no labels, single style)
    sns.set(style="white", context="talk")
    plt.figure(figsize=(10, 8))
    plt.scatter(
        pca_df["PC1"],
        pca_df["PC2"],
        s=80,
        edgecolor="k",
        linewidth=0.5,
        alpha=0.9,
        c="steelblue",
    )

    for axis in ['x', 'y']:
        plt.gca().tick_params(axis=axis, which='both', direction='out', length=6, width=1,
                              colors='black', bottom=True, left=True)

    explained_var = pca.explained_variance_ratio_ * 100
    plt.xlabel(f"PC1 ({explained_var[0]:.1f}%)", fontsize=22, weight="bold")
    plt.ylabel(f"PC2 ({explained_var[1]:.1f}%)", fontsize=22, weight="bold")
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)
    sns.despine()
    plt.tight_layout()

    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_fig) if os.path.dirname(output_fig) else '.'
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created output directory: {output_dir}")

    plt.savefig(output_fig, dpi=600)
    print(f"Saved figure to: {output_fig}")
    plt.show()

    # Save PCA results
    pca_df.to_csv(output_csv, index=True)
    print(f"Saved PCA results to: {output_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="Perform PCA on similarity matrix and plot results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python pca_plot.py --df similarity_matrix.csv
  python pca_plot.py -d matrix.csv -o output.png
        """
    )

    parser.add_argument(
        "--df", "-d",
        type=str,
        required=True,
        help="Path to CSV file containing similarity matrix (with index column)"
    )

    parser.add_argument(
        "--output-fig", "-o",
        type=str,
        default="PCA_Tanimoto_Publication.png",
        help="Path to save output figure (default: PCA_Tanimoto_Publication.png)"
    )

    parser.add_argument(
        "--output-csv",
        type=str,
        default="pca_results.csv",
        help="Path to save PCA results CSV (default: pca_results.csv)"
    )

    args = parser.parse_args()

    try:
        perform_pca_plot(
            df_path=args.df,
            output_fig=args.output_fig,
            output_csv=args.output_csv
        )
    except Exception as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())

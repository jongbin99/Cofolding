#!/usr/bin/env python3
"""Calculate and plot rolling hit rate curves from Excel data.

This script calculates rolling hit rates with Wilson 95% confidence intervals
from Excel files containing score and active/inactive labels.

Usage
-----

    python hitrate.py \\
        --file data.xlsx \\
        --sheet-name "D4" \\
        --score-col "L-pLDDT score" \\
        --window 50
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def rolling_hit_rate(
    df,
    score_col="L-pLDDT score",
    label_col="Active",
    window=50,
    higher_is_better=True,
    min_periods=None
):
    
    dd = df[[score_col, label_col]].dropna().copy()
    dd[label_col] = dd[label_col].astype(int)

    # sort so that "best" scores are on the left of the plot
    dd = dd.sort_values(score_col, ascending=not higher_is_better).reset_index(drop=True)

    if min_periods is None:
        min_periods = max(25, window // 5)

    # rolling counts and means over the ranked list
    rcount = dd[label_col].rolling(window, min_periods=min_periods).count()
    rmean  = dd[label_col].rolling(window, min_periods=min_periods).mean()

    # Wilson 95% confidence interval for a proportion
    z = 1.96
    n = rcount.astype(float)
    p = rmean.fillna(0.0)
    denom  = 1.0 + (z**2)/n.replace(0, np.nan)
    center = (p + (z**2)/(2*n.replace(0, np.nan))) / denom
    halfw  = (z*np.sqrt((p*(1.0-p)/n.replace(0, np.nan)) + (z**2)/(4*n.replace(0, np.nan)**2))) / denom
    lo95 = (center - halfw).clip(lower=0.0)
    hi95 = (center + halfw).clip(upper=1.0)

    out = pd.DataFrame({
        "score": dd[score_col].to_numpy(),
        "hit_rate": rmean.to_numpy(),
        "lo95": lo95.to_numpy(),
        "hi95": hi95.to_numpy(),
        "n": n.to_numpy()
    })

    # Drop positions before the first complete-ish window
    out = out[~np.isnan(out["hit_rate"])].reset_index(drop=True)
    return out


def plot_rolling_hit_rate(
    roll_df,
    baseline=None,
    method_name="AF3",
    color="tab:red",
    score_label="Score",
    title="Rolling hit-rate vs. score",
    output_path=None
):
    """
    Plot rolling hit-rate with Wilson 95% CI and optional baseline line.
    
    Parameters
    ----------
    roll_df : pd.DataFrame
        DataFrame with columns: score, hit_rate, lo95, hi95, n
    baseline : float, optional
        Baseline hit rate to plot as horizontal line
    method_name : str, default="AF3"
        Method name for legend
    color : str, default="tab:red"
        Color for plot
    score_label : str, default="Score"
        X-axis label
    title : str, default="Rolling hit-rate vs. score"
        Plot title
    output_path : str, optional
        Path to save plot (if None, displays plot)
    """
    plt.rcParams.update({
        "font.size": 14,
        "axes.labelsize": 16,
        "axes.titlesize": 18,
        "legend.fontsize": 14,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12
    })

    fig = plt.figure(figsize=(8,6))
    x = roll_df["score"].to_numpy()
    y = roll_df["hit_rate"].to_numpy()
    lo = roll_df["lo95"].to_numpy()
    hi = roll_df["hi95"].to_numpy()

    plt.plot(x, y, color=color, linewidth=2, label=f"{method_name} (rolling)")
    plt.fill_between(x, lo, hi, color=color, alpha=0.2, linewidth=0)

    if baseline is not None:
        plt.axhline(baseline, linestyle="--", color="gray", linewidth=2,
                    label=f"Random baseline = {baseline:.3f}")

    plt.xlabel(score_label)
    plt.ylabel("Hit rate (actives fraction)")
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"[INFO] Plot saved to: {output_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Calculate and plot rolling hit rate curves from Excel data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--file', required=True,
        help='Path to Excel file (must have score and Active columns)'
    )
    parser.add_argument(
        '--sheet-name', default=None,
        help='Sheet name in Excel file (default: first sheet)'
    )
    parser.add_argument(
        '--score-col', default="L-pLDDT score",
        help='Column name for scores'
    )
    parser.add_argument(
        '--label-col', default="Active",
        help='Column name for active/inactive labels'
    )
    parser.add_argument(
        '--window', type=int, default=50,
        help='Rolling window size for hit rate calculation'
    )
    parser.add_argument(
        '--higher-is-better', action='store_true', default=True,
        help='Higher scores are better (default: True)'
    )
    parser.add_argument(
        '--lower-is-better', action='store_true',
        help='Lower scores are better (overrides --higher-is-better)'
    )
    parser.add_argument(
        '--method-name', default="AF3",
        help='Method name for plot legend'
    )
    parser.add_argument(
        '--color', default="tab:red",
        help='Plot color'
    )
    parser.add_argument(
        '--score-label', default=None,
        help='X-axis label (default: same as score column name)'
    )
    parser.add_argument(
        '--title', default=None,
        help='Plot title (default: auto-generated)'
    )
    parser.add_argument(
        '--output', default=None,
        help='Output file path for plot (if not provided, displays plot)'
    )
    parser.add_argument(
        '--show-baseline', action='store_true',
        help='Show random baseline line'
    )
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.file).exists():
        print(f"[ERROR] File not found: {args.file}", file=sys.stderr)
        return 1
    
    # Determine if higher is better
    higher_is_better = args.higher_is_better and not args.lower_is_better
    
    # Load Excel file
    try:
        print(f"[INFO] Loading data from: {args.file}")
        df = pd.read_excel(args.file, sheet_name=args.sheet_name)
        print(f"[INFO] Loaded {len(df)} rows")
    except Exception as e:
        print(f"[ERROR] Failed to read Excel file: {e}", file=sys.stderr)
        return 1
    
    # Validate columns
    if args.score_col not in df.columns:
        print(f"[ERROR] Score column '{args.score_col}' not found. Available columns: {list(df.columns)}", file=sys.stderr)
        return 1
    
    if args.label_col not in df.columns:
        print(f"[ERROR] Label column '{args.label_col}' not found. Available columns: {list(df.columns)}", file=sys.stderr)
        return 1
    
    # Calculate rolling hit rate
    print(f"[INFO] Calculating rolling hit rate (window={args.window})...")
    roll_df = rolling_hit_rate(
        df=df,
        score_col=args.score_col,
        label_col=args.label_col,
        window=args.window,
        higher_is_better=higher_is_better
    )
    
    print(f"[INFO] Calculated hit rates for {len(roll_df)} positions")
    print(f"[INFO] Hit rate range: {roll_df['hit_rate'].min():.3f} - {roll_df['hit_rate'].max():.3f}")
    
    # Calculate baseline if requested
    baseline = None
    if args.show_baseline:
        active_data = df[[args.label_col]].dropna().copy()
        active_data[args.label_col] = active_data[args.label_col].astype(int)
        baseline = active_data[args.label_col].mean()
        print(f"[INFO] Baseline hit rate: {baseline:.3f}")
    
    # Generate labels
    score_label = args.score_label if args.score_label else args.score_col
    title = args.title if args.title else f"Rolling hit-rate vs. {score_label}"
    
    # Plot
    plot_rolling_hit_rate(
        roll_df=roll_df,
        baseline=baseline,
        method_name=args.method_name,
        color=args.color,
        score_label=score_label,
        title=title,
        output_path=args.output
    )
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
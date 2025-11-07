#!/usr/bin/env python3
"""Plot correlation scatter plots with color coding from Excel data.

This script creates scatter plots with correlation analysis, allowing flexible
selection of x and y columns. Points are color-coded based on conditions:
- Black: No/zero affinity data
- Green: With affinity data
- Red: COM distance > threshold (among rows with affinity data)

Usage
-----

    # Basic usage with column names
    python plot_correlation.py \\
        --input data.xlsx \\
        --x-col "Boltz_LRMSD" \\
        --y-col "DOCK_LRMSD"

    # Using column indices (0-based)
    python plot_correlation.py \\
        --input data.xlsx \\
        --x-col 5 \\
        --y-col 2 \\
        --condition-col 1 \\
        --com-cols 6,7,8,9

    # Save to file
    python plot_correlation.py \\
        --input data.xlsx \\
        --x-col "Boltz_LRMSD" \\
        --y-col "DOCK_LRMSD" \\
        --output correlation_plot.png
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import linregress, pearsonr


def parse_column_spec(col_spec, df):
    """
    Parse column specification (either name or index).
    
    Parameters
    ----------
    col_spec : str or int
        Column name or index (0-based)
    df : pd.DataFrame
        DataFrame to extract column from
        
    Returns
    -------
    pd.Series
        Column data
    """
    if isinstance(col_spec, int):
        return df.iloc[:, col_spec]
    elif col_spec.isdigit():
        return df.iloc[:, int(col_spec)]
    else:
        if col_spec in df.columns:
            return df[col_spec]
        else:
            raise ValueError(f"Column '{col_spec}' not found in dataframe. Available columns: {list(df.columns)}")


def plot_correlation(
    input_path,
    x_col,
    y_col,
    condition_col=None,
    com_cols=None,
    com_threshold=2.5,
    output_path=None,
    x_label=None,
    y_label=None,
    title=None,
    figsize=(8, 6),
    dpi=300
):
    """
    Plot correlation scatter plot with color coding.
    
    Parameters
    ----------
    input_path : str
        Path to Excel file
    x_col : str or int
        Column name or index for x-axis
    y_col : str or int
        Column name or index for y-axis
    condition_col : str or int, optional
        Column name or index for condition/affinity data
    com_cols : list of str or int, optional
        Column names or indices for COM distance columns
    com_threshold : float, default=2.5
        Threshold for COM distance (points exceeding this will be red)
    output_path : str, optional
        Path to save plot (if None, displays plot)
    x_label : str, optional
        X-axis label (auto-generated if None)
    y_label : str, optional
        Y-axis label (auto-generated if None)
    title : str, optional
        Plot title (auto-generated with correlation if None)
    figsize : tuple, default=(8, 6)
        Figure size
    dpi : int, default=300
        Resolution for saved plots
    """
    # Load file
    print(f"[INFO] Loading data from: {input_path}")
    try:
        data = pd.read_excel(input_path)
    except Exception as e:
        print(f"[ERROR] Failed to read Excel file: {e}", file=sys.stderr)
        return 1
    
    print(f"[INFO] Loaded {len(data)} rows, {len(data.columns)} columns")
    
    # Extract columns
    try:
        x = parse_column_spec(x_col, data)
        y = parse_column_spec(y_col, data)
    except ValueError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1
    
    # Get column names for labels
    if isinstance(x_col, int) or (isinstance(x_col, str) and x_col.isdigit()):
        x_name = f"Column_{x_col}"
    else:
        x_name = x_col
    
    if isinstance(y_col, int) or (isinstance(y_col, str) and y_col.isdigit()):
        y_name = f"Column_{y_col}"
    else:
        y_name = y_col
    
    # Build dataframe
    df = pd.DataFrame({
        "x": x,
        "y": y
    })
    
    # Add condition column if specified
    if condition_col is not None:
        try:
            df["Condition_raw"] = parse_column_spec(condition_col, data)
            df["Condition"] = pd.to_numeric(df["Condition_raw"], errors="coerce")
        except ValueError as e:
            print(f"[WARN] Could not parse condition column: {e}", file=sys.stderr)
            df["Condition"] = None
    else:
        df["Condition"] = None
    
    # Add COM columns if specified
    if com_cols is not None:
        com_list = [c.strip() for c in com_cols.split(",")]
        for i, com_col in enumerate(com_list):
            try:
                df[f"COM_{i}"] = pd.to_numeric(
                    parse_column_spec(com_col, data),
                    errors="coerce"
                )
            except ValueError as e:
                print(f"[WARN] Could not parse COM column {com_col}: {e}", file=sys.stderr)
                df[f"COM_{i}"] = None
    else:
        # No COM columns specified
        pass
    
    # Drop rows lacking x or y
    initial_len = len(df)
    df.dropna(subset=["x", "y"], inplace=True)
    print(f"[INFO] Dropped {initial_len - len(df)} rows with missing x or y values")
    
    # Create masks
    if df["Condition"].notna().any():
        green_mask = (df["Condition"].notna()) & (df["Condition"] != 0)
    else:
        green_mask = pd.Series([False] * len(df), index=df.index)
    
    # Red mask: COM distance > threshold (only among rows with affinity data)
    if com_cols is not None and any(f"COM_{i}" in df.columns for i in range(len(com_list))):
        com_conditions = []
        for i in range(len(com_list)):
            if f"COM_{i}" in df.columns:
                com_conditions.append(df[f"COM_{i}"] > com_threshold)
        if com_conditions:
            com_mask = pd.concat(com_conditions, axis=1).any(axis=1)
            red_mask = com_mask & (df["Condition"].notna())
        else:
            red_mask = pd.Series([False] * len(df), index=df.index)
    else:
        red_mask = pd.Series([False] * len(df), index=df.index)
    
    print(f"[INFO] {green_mask.sum()} points will be colored GREEN (with affinity data)")
    print(f"[INFO] {red_mask.sum()} points will be colored RED (COM distance > {com_threshold} Å)")
    
    # Correlation & regression (all points)
    valid_mask = df["x"].notna() & df["y"].notna()
    if valid_mask.sum() < 2:
        print(f"[ERROR] Not enough valid data points for correlation (need at least 2, got {valid_mask.sum()})", file=sys.stderr)
        return 1
    
    corr_coeff, p_value = pearsonr(df.loc[valid_mask, "x"], df.loc[valid_mask, "y"])
    slope, intercept, r_value, p_val_reg, std_err = linregress(
        df.loc[valid_mask, "x"],
        df.loc[valid_mask, "y"]
    )
    
    x_min, x_max = df["x"].min(), df["x"].max()
    x_vals = np.linspace(x_min, x_max, 500)
    y_vals = slope * x_vals + intercept
    
    # Generate labels
    if x_label is None:
        x_label = x_name
    if y_label is None:
        y_label = y_name
    if title is None:
        title = f"Pearson r = {corr_coeff:.3f}, p = {p_value:.2e}"
    
    # Plot
    plt.figure(figsize=figsize, dpi=dpi)
    plt.scatter(
        df.loc[~green_mask, "x"],
        df.loc[~green_mask, "y"],
        color="black",
        s=30,
        alpha=0.5,
        label="No/zero affinity data",
        zorder=1
    )
    plt.scatter(
        df.loc[green_mask, "x"],
        df.loc[green_mask, "y"],
        color="green",
        s=30,
        alpha=0.8,
        label="With affinity data",
        zorder=2
    )
    plt.scatter(
        df.loc[red_mask, "x"],
        df.loc[red_mask, "y"],
        color="red",
        s=40,
        alpha=0.9,
        label=f"COM distance > {com_threshold} Å",
        zorder=3
    )
    
    plt.plot(x_vals, y_vals, color="blue", linewidth=2, label="Regression line", zorder=2)
    
    plt.xlabel(x_label, fontsize=14, fontweight="bold", fontname="Arial")
    plt.ylabel(y_label, fontsize=14, fontweight="bold", fontname="Arial")
    plt.title(title, fontsize=14, fontname="Arial")
    plt.tick_params(axis='both', labelsize=14)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        print(f"[INFO] Plot saved to: {output_path}")
    else:
        plt.show()
    
    print(f"[INFO] Correlation: r = {corr_coeff:.3f}, p = {p_value:.2e}")
    print(f"[INFO] Regression: y = {slope:.3f}x + {intercept:.3f}")
    
    return 0


def main():
    parser = argparse.ArgumentParser(
        description='Plot correlation scatter plots with color coding from Excel data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--input', required=True,
                       help='Path to Excel file')
    parser.add_argument('--x-col', required=True,
                       help='X-axis column (name or 0-based index)')
    parser.add_argument('--y-col', required=True,
                       help='Y-axis column (name or 0-based index)')
    parser.add_argument('--condition-col', default=None,
                       help='Condition/affinity column (name or index, optional)')
    parser.add_argument('--com-cols', default=None,
                       help='COM distance columns, comma-separated (e.g., "6,7,8,9", optional)')
    parser.add_argument('--com-threshold', type=float, default=2.5,
                       help='COM distance threshold for red coloring (Å)')
    parser.add_argument('--output', default=None,
                       help='Output file path (if not provided, displays plot)')
    parser.add_argument('--x-label', default=None,
                       help='X-axis label (auto-generated if not provided)')
    parser.add_argument('--y-label', default=None,
                       help='Y-axis label (auto-generated if not provided)')
    parser.add_argument('--title', default=None,
                       help='Plot title (auto-generated with correlation if not provided)')
    parser.add_argument('--figsize', type=float, nargs=2, default=[8, 6],
                       metavar=('WIDTH', 'HEIGHT'),
                       help='Figure size (width, height)')
    parser.add_argument('--dpi', type=int, default=300,
                       help='Resolution for saved plots')
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.input).exists():
        print(f"[ERROR] Input file not found: {args.input}", file=sys.stderr)
        return 1
    
    # Convert column specs (try to parse as int if possible)
    x_col = args.x_col
    y_col = args.y_col
    condition_col = args.condition_col
    
    # Run plotting
    return plot_correlation(
        input_path=args.input,
        x_col=x_col,
        y_col=y_col,
        condition_col=condition_col,
        com_cols=args.com_cols,
        com_threshold=args.com_threshold,
        output_path=args.output,
        x_label=args.x_label,
        y_label=args.y_label,
        title=args.title,
        figsize=tuple(args.figsize),
        dpi=args.dpi
    )


if __name__ == "__main__":
    sys.exit(main())


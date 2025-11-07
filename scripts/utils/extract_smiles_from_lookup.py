#!/usr/bin/env python3
"""Extract SMILES and IDs from lookup table based on filters.

This script extracts SMILES and IDs from a lookup CSV file (e.g., Mols_labelled_hits.csv)
and writes them to a .smi file based on specified filters.

Usage
-----

    # Extract AmpC Hits
    python extract_smiles_from_lookup.py \\
        --lookup /path/to/Mols_labelled_hits.csv \\
        --output ampc_hits.smi \\
        --target AmpC \\
        --hits-status Hits

    # Extract all hits (any target)
    python extract_smiles_from_lookup.py \\
        --lookup /path/to/Mols_labelled_hits.csv \\
        --output all_hits.smi \\
        --hits-status Hits
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


def extract_smiles(lookup_path, output_path, target=None, hits_status=None):
    """
    Extract SMILES and IDs from lookup table based on filters.
    
    Parameters
    ----------
    lookup_path : str
        Path to lookup CSV file
    output_path : str
        Path to output .smi file
    target : str, optional
        Filter by Target column (e.g., 'AmpC')
    hits_status : str, optional
        Filter by Hits_or_nonhits column (e.g., 'Hits', 'Non-hits')
    
    Returns
    -------
    int
        Number of entries extracted
    """
    # Load lookup table
    try:
        df = pd.read_csv(lookup_path)
    except Exception as e:
        print(f"[ERROR] Failed to read lookup file: {e}", file=sys.stderr)
        return 0
    
    print(f"[INFO] Loaded {len(df)} entries from lookup table")
    
    # Find column names (case-insensitive)
    id_col = None
    smiles_col = None
    target_col = None
    hits_col = None
    
    for col in df.columns:
        col_upper = col.upper()
        if col_upper in ['ID', 'ZINC_ID', 'MOLECULE_ID', 'ZINC', 'ZINCID']:
            id_col = col
        elif col_upper in ['SMILES', 'SMILES_STRING', 'SMILES_STR']:
            smiles_col = col
        elif col_upper == 'TARGET':
            target_col = col
        elif col_upper in ['HITS_OR_NONHITS', 'HITS_OR_NON_HITS', 'HIT_STATUS']:
            hits_col = col
    
    # Check required columns
    if id_col is None:
        print(f"[ERROR] Could not find ID column. Available columns: {list(df.columns)}", file=sys.stderr)
        return 0
    
    if smiles_col is None:
        print(f"[ERROR] Could not find SMILES column. Available columns: {list(df.columns)}", file=sys.stderr)
        return 0
    
    # Apply filters
    filtered_df = df.copy()
    filters_applied = []
    
    if target is not None:
        if target_col is None:
            print(f"[ERROR] Target filter specified but no 'Target' column found. Available columns: {list(df.columns)}", file=sys.stderr)
            return 0
        filtered_df = filtered_df[filtered_df[target_col].astype(str).str.strip().str.lower() == target.strip().lower()]
        filters_applied.append(f"Target='{target}'")
    
    if hits_status is not None:
        if hits_col is None:
            print(f"[ERROR] Hits status filter specified but no 'Hits_or_nonhits' column found. Available columns: {list(df.columns)}", file=sys.stderr)
            return 0
        filtered_df = filtered_df[filtered_df[hits_col].astype(str).str.strip().str.lower() == hits_status.strip().lower()]
        filters_applied.append(f"Hits_or_nonhits='{hits_status}'")
    
    # Drop rows with missing SMILES or ID
    filtered_df = filtered_df.dropna(subset=[smiles_col, id_col])
    
    # Remove duplicates (keep first)
    filtered_df = filtered_df.drop_duplicates(subset=[smiles_col, id_col], keep='first')
    
    if len(filtered_df) == 0:
        print(f"[WARN] No entries found matching filters: {', '.join(filters_applied)}", file=sys.stderr)
        return 0
    
    # Write to .smi file (SMILES ID format, space-separated)
    try:
        with open(output_path, 'w', encoding='utf-8') as f:
            for _, row in filtered_df.iterrows():
                smiles = str(row[smiles_col]).strip()
                mol_id = str(row[id_col]).strip()
                f.write(f"{smiles} {mol_id}\n")
    except Exception as e:
        print(f"[ERROR] Failed to write output file: {e}", file=sys.stderr)
        return 0
    
    print(f"[INFO] Filters applied: {', '.join(filters_applied) if filters_applied else 'None'}")
    print(f"[INFO] Extracted {len(filtered_df)} entries to {output_path}")
    
    return len(filtered_df)


def main():
    parser = argparse.ArgumentParser(
        description='Extract SMILES and IDs from lookup table based on filters',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--lookup', required=True,
                       help='Path to lookup CSV file (e.g., Mols_labelled_hits.csv)')
    parser.add_argument('--output', required=True,
                       help='Output .smi file path')
    parser.add_argument('--target', default=None,
                       help='Filter by Target column value (e.g., "AmpC")')
    parser.add_argument('--hits-status', default=None,
                       choices=['Hits', 'Non-hits'],
                       help='Filter by Hits_or_nonhits column value')
    
    args = parser.parse_args()
    
    # Validate input file
    if not Path(args.lookup).exists():
        print(f"[ERROR] Lookup file not found: {args.lookup}", file=sys.stderr)
        return 1
    
    # Extract SMILES
    count = extract_smiles(
        lookup_path=args.lookup,
        output_path=args.output,
        target=args.target,
        hits_status=args.hits_status
    )
    
    if count == 0:
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())


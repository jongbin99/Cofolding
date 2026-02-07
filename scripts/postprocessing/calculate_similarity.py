#!/usr/bin/env python3
"""Calculate molecular similarity matrices from SMILES files.

This script computes pairwise similarity between molecules using either:
1. Tanimoto similarity (ECFP4 fingerprints)
2. Maximum Common Substructure (MCS) percentage overlap

Input files should be whitespace-delimited .smi files with columns: SMILES ID

Usage
-----

    # Tanimoto similarity
    python calculate_similarity.py \\
        --input /path/to/molecules.smi \\
        --output similarity_matrix.csv \\
        --method tanimoto

    # MCS percentage overlap
    python calculate_similarity.py \\
        --input /path/to/molecules.smi \\
        --output mcs_matrix.csv \\
        --method mcs
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem import rdFMCS
from rdkit.Chem import rdFingerprintGenerator

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_molecules_from_smi(file_path: str) -> Tuple[List[Chem.Mol], List[str]]:
    """Load molecules and IDs from a .smi file.
    
    Parameters
    ----------
    file_path : str
        Path to whitespace-delimited .smi file (col1=SMILES, col2=ID)
        
    Returns
    -------
    Tuple[List[Chem.Mol], List[str]]
        List of RDKit molecule objects and corresponding IDs
    """
    try:
        raw = pd.read_csv(
            file_path,
            sep=r"\s+",
            header=None,
            usecols=[0, 1],
            names=["smiles", "id"],
            engine="python"
        )
    except Exception as e:
        logger.error(f"Error reading SMILES file: {e}")
        sys.exit(1)
    
    # Drop missing and deduplicate by SMILES (keep first)
    df_in = (raw.dropna(subset=["smiles", "id"])
             .drop_duplicates(subset=["smiles"], keep="first")
             .reset_index(drop=True))
    
    # Build molecules; filter invalid
    mols = []
    ids = []
    skipped_count = 0
    
    for smi, _id in zip(df_in["smiles"], df_in["id"]):
        mol = Chem.MolFromSmiles(str(smi))
        if mol is None:
            skipped_count += 1
            logger.debug(f"Skipping invalid SMILES for ID {_id}: {smi[:50]}...")
            continue
        mols.append(mol)
        ids.append(str(_id))
    
    if skipped_count > 0:
        logger.warning(f"Skipped {skipped_count} invalid SMILES entries")
    
    logger.info(f"Successfully loaded {len(mols)} valid molecules from {len(df_in)} entries")
    return mols, ids


def calculate_tanimoto_similarity(mols: List[Chem.Mol], ids: List[str], 
                                  radius: int = 2, n_bits: int = 2048) -> pd.DataFrame:
    """Calculate pairwise Tanimoto similarity matrix.
    
    Parameters
    ----------
    mols : List[Chem.Mol]
        List of RDKit molecule objects
    ids : List[str]
        List of molecule IDs
    radius : int, default=2
        Morgan fingerprint radius (ECFP4 = radius 2)
    n_bits : int, default=2048
        Number of bits in the fingerprint
        
    Returns
    -------
    pd.DataFrame
        Symmetric similarity matrix with IDs as row/column labels
    """
    logger.info(f"Generating Morgan fingerprints (radius={radius}, nBits={n_bits})...")
    # Use the new MorganGenerator API (replaces deprecated GetMorganFingerprintAsBitVect)
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    fps = [morgan_gen.GetFingerprint(m) for m in mols]
    
    logger.info("Computing pairwise Tanimoto similarity...")
    n = len(fps)
    sim_matrix = np.zeros((n, n), dtype=float)
    
    for i in range(n):
        sim_matrix[i, i] = 1.0  # Self-similarity = 1.0
        for j in range(i + 1, n):
            s = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            sim_matrix[i, j] = s
            sim_matrix[j, i] = s  # Symmetric matrix
    
    # Create DataFrame labeled by IDs
    df = pd.DataFrame(sim_matrix, columns=ids, index=ids)
    return df


def calculate_mcs_overlap(mols: List[Chem.Mol], ids: List[str]) -> pd.DataFrame:
    """Calculate pairwise MCS percentage overlap.
    
    Parameters
    ----------
    mols : List[Chem.Mol]
        List of RDKit molecule objects
    ids : List[str]
        List of molecule IDs
        
    Returns
    -------
    pd.DataFrame
        Symmetric similarity matrix with IDs as row/column labels
    """
    import itertools
    
    logger.info("Computing pairwise MCS percentage overlap...")
    n = len(mols)
    mcs_matrix = pd.DataFrame(index=ids, columns=ids)
    
    # Calculate MCS for each pair
    total_pairs = n * (n - 1) // 2
    processed = 0
    
    for i, j in itertools.combinations(range(n), 2):
        m1, m2 = mols[i], mols[j]
        if m1 and m2:
            try:
                mcs_result = rdFMCS.FindMCS([m1, m2], completeRingsOnly=True)
                mcs_mol = Chem.MolFromSmarts(mcs_result.smartsString)
                mcs_atoms = mcs_mol.GetNumAtoms()
                
                # Percentage overlap: MCS atoms / average of the two molecule sizes
                avg_size = (m1.GetNumAtoms() + m2.GetNumAtoms()) / 2.0
                pct_overlap = (mcs_atoms / avg_size * 100.0) if avg_size > 0 else 0.0
            except Exception as e:
                logger.debug(f"Error computing MCS for {ids[i]} vs {ids[j]}: {e}")
                pct_overlap = 0.0
        else:
            pct_overlap = 0.0
        
        mcs_matrix.iloc[i, j] = pct_overlap
        mcs_matrix.iloc[j, i] = pct_overlap  # Symmetric matrix
        
        processed += 1
        if processed % 100 == 0:
            logger.info(f"Processed {processed}/{total_pairs} pairs...")
    
    # Set diagonal to 100% for self-comparison
    for i in range(n):
        mcs_matrix.iloc[i, i] = 100.0
    
    logger.info(f"Finished computing MCS for {total_pairs} pairs")
    return mcs_matrix


def calculate_average_excluding_diagonal(df: pd.DataFrame) -> float:
    """Calculate average similarity excluding diagonal (self-comparison) entries.
    
    Parameters
    ----------
    df : pd.DataFrame
        Symmetric similarity matrix
        
    Returns
    -------
    float
        Average similarity excluding diagonal entries
    """
    # Convert to numpy array and mask diagonal
    matrix = df.values
    n = len(matrix)
    
    # Create mask to exclude diagonal
    mask = ~np.eye(n, dtype=bool)
    
    # Calculate mean excluding diagonal
    avg = matrix[mask].mean()
    return avg


def main():
    parser = argparse.ArgumentParser(
        description='Calculate molecular similarity matrices from SMILES files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--input', required=True,
                       help='Input .smi file (columns: SMILES ID)')
    parser.add_argument('--output', required=True,
                       help='Output CSV file path for similarity matrix')
    parser.add_argument('--method', required=True,
                       choices=['tanimoto', 'mcs'],
                       help='Similarity metric: tanimoto (ECFP4) or mcs (MCS % overlap)')
    
    # Optional arguments
    parser.add_argument('--radius', type=int, default=2,
                       help='Morgan fingerprint radius for Tanimoto (ECFP4=2)')
    parser.add_argument('--n-bits', type=int, default=2048,
                       help='Number of fingerprint bits for Tanimoto')
    parser.add_argument('--decimals', type=int, default=2,
                       help='Number of decimal places in output')
    parser.add_argument('--log-level', type=str, default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                       help='Set the logging level')
    
    args = parser.parse_args()
    
    # Set logging level
    logger.setLevel(args.log_level)
    
    # Validate input file
    if not Path(args.input).exists():
        logger.error(f"Input file not found: {args.input}")
        return 1
    
    # Load molecules
    mols, ids = load_molecules_from_smi(args.input)
    
    if len(mols) == 0:
        logger.error("No valid molecules found in input file")
        return 1
    
    # Calculate similarity matrix
    logger.info(f"Computing {args.method} similarity matrix...")
    if args.method == 'tanimoto':
        df = calculate_tanimoto_similarity(mols, ids, radius=args.radius, n_bits=args.n_bits)
    else:  # mcs
        df = calculate_mcs_overlap(mols, ids)
    
    # Calculate and report average excluding diagonal (self-comparisons)
    avg_excluding_diagonal = calculate_average_excluding_diagonal(df)
    metric_name = "Tanimoto similarity" if args.method == 'tanimoto' else "MCS%"
    logger.info(f"Average {metric_name} (excluding self-comparisons): {avg_excluding_diagonal:.{args.decimals}f}")
    
    # Round values
    df_rounded = df.round(args.decimals)
    
    # Save to CSV
    try:
        df_rounded.to_csv(args.output)
        logger.info(f"Successfully wrote similarity matrix to {args.output}")
        logger.info(f"Matrix size: {len(df)}x{len(df)}")
    except Exception as e:
        logger.error(f"Error writing output file: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())


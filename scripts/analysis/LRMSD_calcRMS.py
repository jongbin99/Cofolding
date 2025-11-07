import os
import csv
import argparse
import logging
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolAlign import CalcRMS
from rdkit.Chem.rdmolops import RemoveHs, RemoveStereochemistry
from rdkit.Geometry import Point3D

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def get_com(mol):
    """Calculate center of mass for a molecule"""
    conf = mol.GetConformer()
    coords = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
    com_array = np.mean(coords, axis=0)
    return Point3D(float(com_array[0]), float(com_array[1]), float(com_array[2]))

def calculate_com_distance(mol1, mol2):
    """Calculate distance between centers of mass of two molecules"""
    com1 = get_com(mol1)
    com2 = get_com(mol2)
    return np.sqrt((com1.x - com2.x)**2 + (com1.y - com2.y)**2 + (com1.z - com2.z)**2)

def process_ligand_pair(ref_file, ref_dir, pred_dir, df, args):
    """Process a single ligand pair and return RMSD results"""
    
    if not ref_file.endswith(".pdb"):
        return None

    ref_prefix = os.path.splitext(ref_file)[0]
    ref_path = os.path.join(ref_dir, ref_file)
    pred_path = os.path.join(pred_dir, f"{ref_prefix}.pdb")
    
    # Check if predicted file exists
    if not os.path.exists(pred_path):
        logger.warning(f"Predicted file not found: {pred_path}")
        return None

    # Get SMILES from dataframe
    smiles_row = df.loc[df['Dataset_ID'] == ref_prefix, 'SMILES']
    if smiles_row.empty:
        logger.warning(f"No SMILES found for {ref_prefix}")
        return None
    
    smiles = smiles_row.values[0]
    tmplt_mol = Chem.MolFromSmiles(smiles)
    if tmplt_mol is None:
        logger.error(f"Invalid SMILES for {ref_prefix}: {smiles}")
        return None

    # Load molecules
    try:
        reference = Chem.MolFromPDBFile(ref_path, sanitize=False, removeHs=False)
        if reference is None:
            logger.error(f"Failed to load reference PDB: {ref_path}")
            return None
    except Exception as e:
        logger.error(f"Error loading reference {ref_path}: {e}")
        return None

    try:
        predicted = Chem.MolFromPDBFile(pred_path, sanitize=False, removeHs=False)
        if predicted is None:
            logger.error(f"Failed to load predicted PDB: {pred_path}")
            return None
    except Exception as e:
        logger.error(f"Error loading predicted {pred_path}: {e}")
        return None

    # Process molecules
    try:
        RemoveStereochemistry(reference)
        reference = RemoveHs(reference, sanitize=False)
        RemoveStereochemistry(predicted)
        predicted = RemoveHs(predicted, sanitize=False)
    except Exception as e:
        logger.error(f"Error removing stereochemistry for {ref_prefix}: {e}")
        return None

    # Calculate COM distance
    try:
        com_dist = calculate_com_distance(reference, predicted)   
    except Exception as e:
        logger.error(f"Error calculating COM distance for {ref_prefix}: {e}")
        return None

    # Assign bond orders from template
    try:
        reference = AllChem.AssignBondOrdersFromTemplate(tmplt_mol, reference)
        predicted = AllChem.AssignBondOrdersFromTemplate(tmplt_mol, predicted)
    except Exception as e:
        logger.error(f"Error assigning bond orders for {ref_prefix}: {e}")
        return None

    # Align atom ordering by substructure match
    try:
        match = predicted.GetSubstructMatch(reference)
        if not match:
            logger.warning(f"No substructure match found for {ref_prefix}")
            return None
        predicted = Chem.RenumberAtoms(predicted, list(match))
    except Exception as e:
        logger.error(f"Error aligning atoms for {ref_prefix}: {e}")
        return None

    # Compute RMSD
    try:
        calc_rmsd = CalcRMS(predicted, reference)  # CalcRMS from RDKit
        
        logger.info(f"Success {ref_prefix}: CalcRMSD={calc_rmsd:.3f}, COM_dist={com_dist:.3f}")
        return [ref_prefix, calc_rmsd, com_dist]
        
    except Exception as e:
        logger.error(f"Error calculating RMSD for {ref_prefix}: {e}")
        return None

def main():
    parser = argparse.ArgumentParser(
        description='Calculate ligand RMSD between experimental and predicted structures',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--excel', required=True,
                       help='Path to Excel file with SMILES (must have Dataset_ID and SMILES columns)')
    parser.add_argument('--ref-dir', required=True,
                       help='Directory containing reference/experimental PDB files')
    parser.add_argument('--pred-dir', required=True,
                       help='Directory containing predicted/docked PDB files')
    parser.add_argument('--output', required=True,
                       help='Output CSV file path')
    
    # Optional arguments
    parser.add_argument('--log-level', type=str, default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                       help='Set the logging level')
    
    args = parser.parse_args()
    
    # Set logging level
    logger.setLevel(args.log_level)
    
    # Load Excel file
    try:
        df = pd.read_excel(args.excel)
        df['Dataset_ID'] = df['Dataset_ID'].astype(str)
        df['SMILES'] = df['SMILES'].astype(str)
        logger.info(f"Loaded {len(df)} entries from Excel file")
    except Exception as e:
        logger.error(f"Error reading Excel file: {e}")
        return 1
    
    # Process all reference files
    rmsd_results = []
    ref_files = [f for f in os.listdir(args.ref_dir) if f.endswith('.pdb')]
    logger.info(f"Found {len(ref_files)} reference PDB files")
    
    for i, ref_file in enumerate(ref_files, 1):
        logger.info(f"Processing {i}/{len(ref_files)}: {ref_file}")
        result = process_ligand_pair(ref_file, args.ref_dir, args.pred_dir, df, args)
        if result:
            rmsd_results.append(result)
    
    # Write results
    if rmsd_results:
        try:
            with open(args.output, "w", newline="") as f:
                writer = csv.writer(f)
                writer.writerow(["ID", "CalcRMSD", "COM_Distance"])
                writer.writerows(rmsd_results)
            logger.info(f"Successfully wrote {len(rmsd_results)} results to {args.output}")
        except Exception as e:
            logger.error(f"Error writing output file: {e}")
            return 1
    else:
        logger.warning("No successful RMSD calculations to write")
    
    # Summary
    logger.info(f"Summary: {len(rmsd_results)}/{len(ref_files)} successfully processed")
    return 0

if __name__ == "__main__":
    exit(main())
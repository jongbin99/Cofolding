#!/usr/bin/env python3
"""Align predicted PDB structures to reference complexes using PyMOL.

This script replicates the alignment logic from `aligning_docked_exp.ipynb` and
exposes it as a command-line tool. It uses PyMOL to align predicted structures
to reference complexes and optionally extracts ligand coordinates.

Usage
-----

    python align_pdb_structures.py \\
        --complexes-dir /path/to/reference \\
        --predicted-dir /path/to/predicted \\
        --aligned-dir /path/to/output/aligned \\
        --ligand-dir /path/to/output/ligands

The script will:
1. Find matching PDB files in complexes-dir and predicted-dir
2. Align predicted structures to reference complexes using PyMOL
3. Save aligned structures to aligned-dir
4. Optionally extract ligand entries (lines containing ' LIG ') to ligand-dir

Optional flags control whether to extract ligands and the ligand label to search for.
"""

import argparse
import logging
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def align_with_pymol_api(
    complex_path: str,
    predicted_path: str,
    aligned_path: str
) -> bool:
    """Align structures using PyMOL Python API.
    
    Parameters
    ----------
    complex_path : str
        Path to reference complex PDB file
    predicted_path : str
        Path to predicted PDB file
    aligned_path : str
        Output path for aligned structure
        
    Returns
    -------
    bool
        True if alignment succeeded, False otherwise
    """
    try:
        from pymol import cmd
        
        # Initialize PyMOL
        cmd.reinitialize()
        
        # Load structures
        cmd.load(complex_path, "ref")
        cmd.load(predicted_path, "pred")
        
        # Align
        cmd.align("pred", "ref")
        
        # Save aligned structure
        cmd.save(aligned_path, "pred")
        
        # Clean up
        cmd.delete("ref")
        cmd.delete("pred")
        
        return True
    except ImportError:
        logger.debug("PyMOL Python API not available")
        return False
    except Exception as e:
        logger.debug(f"PyMOL API alignment failed: {e}")
        return False


def align_and_extract(
    complexes_dir: str,
    predicted_dir: str,
    aligned_dir: str,
    ligand_dir: Optional[str] = None,
    ligand_label: str = " LIG ",
    use_pymol_api: bool = False
):
    """Align predicted PDB structures to reference complexes and optionally extract ligands.
    
    Parameters
    ----------
    complexes_dir : str
        Directory containing reference/complex PDB files
    predicted_dir : str
        Directory containing predicted PDB files
    aligned_dir : str
        Output directory for aligned structures
    ligand_dir : str, optional
        If provided, extracts ligand entries to this directory
    ligand_label : str, default=" LIG "
        String pattern to identify ligand lines in PDB files
    use_pymol_api : bool, default=False
        If True, use PyMOL Python API instead of subprocess
        
    Returns
    -------
    int
        Number of successfully processed structures
    """
    # Ensure output directories exist
    os.makedirs(aligned_dir, exist_ok=True)
    if ligand_dir:
        os.makedirs(ligand_dir, exist_ok=True)
    
    # Get all PDB files in both directories
    try:
        complex_files = {f for f in os.listdir(complexes_dir) if f.endswith(".pdb")}
        predicted_files = {f for f in os.listdir(predicted_dir) if f.endswith(".pdb")}
    except OSError as e:
        logger.error(f"Error reading directories: {e}")
        return 0
    
    if not complex_files:
        logger.warning(f"No PDB files found in {complexes_dir}")
        return 0
    
    if not predicted_files:
        logger.warning(f"No PDB files found in {predicted_dir}")
        return 0
    
    # Find matching filenames in both folders
    matching_files = complex_files.intersection(predicted_files)
    logger.info(f"Found {len(matching_files)} matching PDB files")
    
    if not matching_files:
        logger.warning("No matching files found between complexes and predicted directories")
        return 0
    
    # Use a temporary directory for PyMOL scripts to avoid clutter
    temp_dir = tempfile.mkdtemp(prefix="pymol_align_")
    
    success_count = 0
    failed_files = []
    
    for idx, pdb_file in enumerate(sorted(matching_files), 1):
        try:
            complex_path = os.path.join(complexes_dir, pdb_file)
            predicted_path = os.path.join(predicted_dir, pdb_file)
            aligned_pdb_path = os.path.join(aligned_dir, f"aligned_{pdb_file}")
            
            logger.info(f"Processing {idx}/{len(matching_files)}: {pdb_file}")
            
            # Try PyMOL API first if requested, otherwise fall back to subprocess
            if use_pymol_api:
                success = align_with_pymol_api(complex_path, predicted_path, aligned_pdb_path)
                if success:
                    # Ligand extraction happens below
                    pass
                else:
                    # Fall back to subprocess
                    logger.debug("PyMOL API failed, trying subprocess method")
                    success = False
            
            if not use_pymol_api or not success:
                # Use subprocess method
                pymol_script = f"""load {complex_path}, complex
load {predicted_path}, predicted
align predicted, complex
save {aligned_pdb_path}, predicted
quit
"""
                
                # Save script as a temporary file
                script_path = os.path.join(temp_dir, f"{pdb_file}.pml")
                with open(script_path, "w") as f:
                    f.write(pymol_script)
                
                # Run PyMOL
                subprocess.run(
                    ["pymol", "-cq", script_path],
                    check=True,
                    capture_output=True,
                    text=True
                )
            
            # Extract ligand entries if requested
            if ligand_dir:
                ligand_output_path = os.path.join(ligand_dir, f"ligand_{pdb_file}")
                with open(aligned_pdb_path, "r") as f:
                    ligand_entries = [line for line in f if ligand_label in line]
                
                # Save extracted ligand data
                with open(ligand_output_path, "w") as f:
                    f.writelines(ligand_entries)
                
                if ligand_entries:
                    logger.debug(f"Extracted {len(ligand_entries)} ligand lines")
                else:
                    logger.warning(f"No ligand entries found (pattern: '{ligand_label}')")
            
            success_count += 1
            
        except subprocess.CalledProcessError as e:
            logger.error(f"PyMOL error for {pdb_file}: {e}")
            if e.stderr:
                logger.debug(f"PyMOL stderr: {e.stderr}")
            failed_files.append(pdb_file)
        except FileNotFoundError:
            logger.error("PyMOL not found. Please install PyMOL or add it to your PATH.")
            return success_count
        except Exception as e:
            logger.error(f"Unexpected error for {pdb_file}: {e}")
            failed_files.append(pdb_file)
    
    # Clean up temporary directory
    try:
        import shutil
        shutil.rmtree(temp_dir)
    except Exception as e:
        logger.debug(f"Could not remove temp directory {temp_dir}: {e}")
    
    # Summary
    logger.info(f"\n{'='*60}")
    logger.info(f"Successfully processed: {success_count}/{len(matching_files)}")
    if failed_files:
        logger.warning(f"Failed files: {len(failed_files)}")
        for failed in failed_files[:10]:  # Show first 10
            logger.warning(f"  - {failed}")
        if len(failed_files) > 10:
            logger.warning(f"  ... and {len(failed_files) - 10} more")
    
    return success_count


def main():
    parser = argparse.ArgumentParser(
        description='Align predicted PDB structures to reference complexes using PyMOL',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--complexes-dir', required=True,
                       help='Directory containing reference/complex PDB files')
    parser.add_argument('--predicted-dir', required=True,
                       help='Directory containing predicted PDB files')
    parser.add_argument('--aligned-dir', required=True,
                       help='Output directory for aligned structures')
    
    # Optional arguments
    parser.add_argument('--ligand-dir', default=None,
                       help='Output directory for extracted ligand coordinates (optional)')
    parser.add_argument('--ligand-label', default=' LIG ',
                       help='String pattern to identify ligand lines in PDB files')
    parser.add_argument('--use-pymol-api', action='store_true',
                       help='Use PyMOL Python API instead of subprocess (requires pymol module)')
    parser.add_argument('--log-level', type=str, default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                       help='Set the logging level')
    
    args = parser.parse_args()
    
    # Set logging level
    logger.setLevel(args.log_level)
    
    # Validate directories
    if not os.path.isdir(args.complexes_dir):
        logger.error(f"Directory not found: {args.complexes_dir}")
        return 1
    
    if not os.path.isdir(args.predicted_dir):
        logger.error(f"Directory not found: {args.predicted_dir}")
        return 1
    
    # Run alignment
    success_count = align_and_extract(
        complexes_dir=args.complexes_dir,
        predicted_dir=args.predicted_dir,
        aligned_dir=args.aligned_dir,
        ligand_dir=args.ligand_dir,
        ligand_label=args.ligand_label,
        use_pymol_api=args.use_pymol_api
    )
    
    return 0 if success_count > 0 else 1


if __name__ == "__main__":
    sys.exit(main())


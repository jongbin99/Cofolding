#!/usr/bin/env python3
"""
Master workflow script for co-folding analysis pipeline.

This script orchestrates the complete analysis workflow:
1. PDB file processing (clean and renumber)
2. Protein RMSD calculation
3. Structure alignment and ligand extraction
4. Ligand RMSD calculation

Usage:
    python scripts/run_analysis_workflow.py --config workflow_config.json
    python scripts/run_analysis_workflow.py --step 1 --base-dir /path/to/pdbs
"""

import argparse
import json
import logging
import os
import subprocess
import sys
from pathlib import Path

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Script paths (relative to repository root)
SCRIPTS_DIR = Path(__file__).parent
REPO_ROOT = SCRIPTS_DIR.parent
PROCESS_PDB = REPO_ROOT / "scripts" / "postprocessing" / "process_pdb_residues.py"
PROTEIN_RMSD = REPO_ROOT / "scripts" / "postprocessing" / "All_protein_RMSD.py"
ALIGN_STRUCTURES = REPO_ROOT / "scripts" / "postprocessing" / "align_pdb_structures.py"
LIGAND_RMSD = REPO_ROOT / "scripts" / "postprocessing" / "LRMSD_calcRMS.py"


def run_step1_process_pdb(base_dir, new_start=3, recursive=False):
    """Step 1: Process PDB files (clean and renumber)."""
    logger.info("=" * 60)
    logger.info("Step 1: Processing PDB files")
    logger.info("=" * 60)
    
    if not os.path.isdir(base_dir):
        logger.error(f"Base directory not found: {base_dir}")
        return False
    
    cmd = [
        sys.executable,
        str(PROCESS_PDB),
        "--base-dir", base_dir,
        "--new-start", str(new_start)
    ]
    
    if recursive:
        cmd.append("--recursive")
    
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Step 1 failed: {result.stderr}")
        return False
    
    logger.info("Step 1 completed successfully")
    return True


def run_step2_protein_rmsd(input_dir, output_csv):
    """Step 2: Calculate protein RMSD."""
    logger.info("=" * 60)
    logger.info("Step 2: Calculating protein RMSD")
    logger.info("=" * 60)
    
    if not os.path.isdir(input_dir):
        logger.error(f"Input directory not found: {input_dir}")
        return False
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    cmd = [
        sys.executable,
        str(PROTEIN_RMSD),
        "--input_dir", input_dir,
        "--output_csv", output_csv
    ]
    
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Step 2 failed: {result.stderr}")
        return False
    
    logger.info(f"Step 2 completed successfully. Output: {output_csv}")
    return True


def run_step3_align_structures(complexes_dir, predicted_dir, aligned_dir, 
                                ligand_dir=None, ligand_label=" LIG "):
    """Step 3: Align structures and extract ligands."""
    logger.info("=" * 60)
    logger.info("Step 3: Aligning structures and extracting ligands")
    logger.info("=" * 60)
    
    if not os.path.isdir(complexes_dir):
        logger.error(f"Complexes directory not found: {complexes_dir}")
        return False
    
    if not os.path.isdir(predicted_dir):
        logger.error(f"Predicted directory not found: {predicted_dir}")
        return False
    
    # Create output directories
    os.makedirs(aligned_dir, exist_ok=True)
    if ligand_dir:
        os.makedirs(ligand_dir, exist_ok=True)
    
    cmd = [
        sys.executable,
        str(ALIGN_STRUCTURES),
        "--complexes-dir", complexes_dir,
        "--predicted-dir", predicted_dir,
        "--aligned-dir", aligned_dir,
        "--ligand-label", ligand_label
    ]
    
    if ligand_dir:
        cmd.extend(["--ligand-dir", ligand_dir])
    
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Step 3 failed: {result.stderr}")
        return False
    
    logger.info("Step 3 completed successfully")
    return True


def run_step4_ligand_rmsd(excel_file, ref_dir, pred_dir, output_csv):
    """Step 4: Calculate ligand RMSD."""
    logger.info("=" * 60)
    logger.info("Step 4: Calculating ligand RMSD")
    logger.info("=" * 60)
    
    if not os.path.isfile(excel_file):
        logger.error(f"Excel file not found: {excel_file}")
        return False
    
    if not os.path.isdir(ref_dir):
        logger.error(f"Reference directory not found: {ref_dir}")
        return False
    
    if not os.path.isdir(pred_dir):
        logger.error(f"Predicted directory not found: {pred_dir}")
        return False
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    
    cmd = [
        sys.executable,
        str(LIGAND_RMSD),
        "--excel", excel_file,
        "--ref-dir", ref_dir,
        "--pred-dir", pred_dir,
        "--output", output_csv
    ]
    
    logger.info(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        logger.error(f"Step 4 failed: {result.stderr}")
        return False
    
    logger.info(f"Step 4 completed successfully. Output: {output_csv}")
    return True


def run_all_steps(config):
    """Run all workflow steps in sequence."""
    logger.info("Starting complete workflow")
    
    success = True
    
    # Step 1: Process PDB files
    if config.get("step1", {}).get("enabled", False):
        step1_config = config["step1"]
        if not run_step1_process_pdb(
            step1_config["base_dir"],
            step1_config.get("new_start", 3),
            step1_config.get("recursive", False)
        ):
            success = False
            if config.get("stop_on_error", True):
                return False
    
    # Step 2: Protein RMSD
    if config.get("step2", {}).get("enabled", False):
        step2_config = config["step2"]
        if not run_step2_protein_rmsd(
            step2_config["input_dir"],
            step2_config["output_csv"]
        ):
            success = False
            if config.get("stop_on_error", True):
                return False
    
    # Step 3: Align structures
    if config.get("step3", {}).get("enabled", False):
        step3_config = config["step3"]
        if not run_step3_align_structures(
            step3_config["complexes_dir"],
            step3_config["predicted_dir"],
            step3_config["aligned_dir"],
            step3_config.get("ligand_dir"),
            step3_config.get("ligand_label", " LIG ")
        ):
            success = False
            if config.get("stop_on_error", True):
                return False
    
    # Step 4: Ligand RMSD
    if config.get("step4", {}).get("enabled", False):
        step4_config = config["step4"]
        if not run_step4_ligand_rmsd(
            step4_config["excel"],
            step4_config["ref_dir"],
            step4_config["pred_dir"],
            step4_config["output"]
        ):
            success = False
    
    return success


def main():
    parser = argparse.ArgumentParser(
        description="Run co-folding analysis workflow",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all steps from config file
  python scripts/run_analysis_workflow.py --config workflow_config.json
  
  # Run a specific step
  python scripts/run_analysis_workflow.py --step 1 --base-dir /path/to/pdbs
  
  # Run step 2
  python scripts/run_analysis_workflow.py --step 2 \\
      --input-dir /path/to/input \\
      --output-csv protein_rmsd.csv
        """
    )
    
    parser.add_argument(
        "--config",
        type=str,
        help="Path to JSON configuration file"
    )
    
    parser.add_argument(
        "--step",
        type=int,
        choices=[1, 2, 3, 4],
        help="Run a specific step (1-4)"
    )
    
    # Step 1 arguments
    parser.add_argument("--base-dir", help="Base directory for PDB processing (Step 1)")
    parser.add_argument("--new-start", type=int, default=3, help="Starting residue number (Step 1)")
    parser.add_argument("--recursive", action="store_true", help="Process recursively (Step 1)")
    
    # Step 2 arguments
    parser.add_argument("--input-dir", help="Input directory for protein RMSD (Step 2)")
    parser.add_argument("--output-csv", help="Output CSV file (Step 2)")
    
    # Step 3 arguments
    parser.add_argument("--complexes-dir", help="Reference/complexes directory (Step 3)")
    parser.add_argument("--predicted-dir", help="Predicted directory (Step 3)")
    parser.add_argument("--aligned-dir", help="Aligned output directory (Step 3)")
    parser.add_argument("--ligand-dir", help="Ligand output directory (Step 3)")
    parser.add_argument("--ligand-label", default=" LIG ", help="Ligand label pattern (Step 3)")
    
    # Step 4 arguments
    parser.add_argument("--excel", help="Excel file with SMILES (Step 4)")
    parser.add_argument("--ref-dir", help="Reference directory (Step 4)")
    parser.add_argument("--pred-dir", help="Predicted directory (Step 4)")
    parser.add_argument("--output", help="Output CSV file (Step 4)")
    
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set logging level"
    )
    
    args = parser.parse_args()
    
    # Set logging level
    logger.setLevel(args.log_level)
    
    # If config file is provided, load it and run
    if args.config:
        if not os.path.isfile(args.config):
            logger.error(f"Config file not found: {args.config}")
            return 1
        
        with open(args.config, 'r') as f:
            config = json.load(f)
        
        success = run_all_steps(config)
        return 0 if success else 1
    
    # Otherwise, run a specific step
    if args.step == 1:
        if not args.base_dir:
            logger.error("--base-dir is required for Step 1")
            return 1
        success = run_step1_process_pdb(args.base_dir, args.new_start, args.recursive)
    
    elif args.step == 2:
        if not args.input_dir or not args.output_csv:
            logger.error("--input-dir and --output-csv are required for Step 2")
            return 1
        success = run_step2_protein_rmsd(args.input_dir, args.output_csv)
    
    elif args.step == 3:
        if not args.complexes_dir or not args.predicted_dir or not args.aligned_dir:
            logger.error("--complexes-dir, --predicted-dir, and --aligned-dir are required for Step 3")
            return 1
        success = run_step3_align_structures(
            args.complexes_dir,
            args.predicted_dir,
            args.aligned_dir,
            args.ligand_dir,
            args.ligand_label
        )
    
    elif args.step == 4:
        if not args.excel or not args.ref_dir or not args.pred_dir or not args.output:
            logger.error("--excel, --ref-dir, --pred-dir, and --output are required for Step 4")
            return 1
        success = run_step4_ligand_rmsd(args.excel, args.ref_dir, args.pred_dir, args.output)
    
    else:
        logger.error("Either --config or --step must be specified")
        parser.print_help()
        return 1
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())


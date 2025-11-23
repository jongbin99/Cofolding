# Cofolding Analysis Workflow

This document describes the complete workflow for analyzing co-folded protein-ligand complexes.

## Workflow Overview

The analysis pipeline consists of four main steps:

```
┌─────────────────────────────────────────────────────────────────┐
│ Step 1: PDB File Processing                                      │
│                                                                   │
│  base_directory ──┐                                              │
│  starting_number ─┘──> process_pdb_residues.py                  │
│  (default: 3)                                                    │
│                                                                   │
│  Output: Cleaned and renumbered PDB files                        │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ Step 2: Protein RMSD Calculation                                │
│                                                                   │
│  input_directory ──┐                                            │
│                    └──> All_protein_RMSD.py                     │
│  output.csv <───────┘                                            │
│                                                                   │
│  Output: CSV with protein RMSD metrics                          │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ Step 3: Structure Alignment and Ligand Extraction               │
│                                                                   │
│  reference_directory ──┐                                        │
│  predicted_directory ──┼──> align_pdb_structures.py            │
│                        └──┐                                     │
│                           │                                     │
│  aligned_predicted_directory <──┘                             │
│  aligned_ligand_directory <────┘                                │
│                                                                   │
│  Output: Aligned structures and extracted ligands               │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│ Step 4: Ligand RMSD Calculation                                 │
│                                                                   │
│  excel.xlsx ─────────┐                                           │
│  reference_directory ┼──> LRMSD_calcRMS.py                      │
│  predicted_directory ┘                                          │
│                                                                   │
│  output.csv <────────┘                                           │
│                                                                   │
│  Output: CSV with ligand RMSD and COM distance                  │
└─────────────────────────────────────────────────────────────────┘
```

## Detailed Workflow Steps

### Step 1: PDB File Processing
**Script:** `scripts/utils/process_pdb_residues.py`

Clean and renumber PDB files by removing specific residues and converting CIF to PDB format.

**Inputs:**
- `base_directory`: Directory containing PDB files to process
- `starting_number`: Residue number to assign to the first remaining residue (default: 3)

**Command:**
```bash
python scripts/utils/process_pdb_residues.py \
    --base-dir base_directory \
    --new-start starting_number
```

**Output:** Modified PDB files in place (cleaned and renumbered)

---

### Step 2: Protein RMSD Calculation
**Script:** `scripts/analysis/All_protein_RMSD.py`

Calculate RMSD between reference and predicted protein structures using PyMOL.

**Inputs:**
- `input_directory`: Directory containing both reference and predicted PDB files
  - Reference files: named with `mac-x` prefix
  - Predicted files: contain `_pred_chainA` in filename

**Command:**
```bash
python scripts/analysis/All_protein_RMSD.py \
    --input_dir input_directory \
    --output_csv output.csv
```

**Output:** CSV file with columns:
- `Protein`: Name of the protein
- `Full Structure RMSD (Å)`: RMSD between full structures
- `Pocket RMSD (Å)`: RMSD between pocket residues
- `Side-Chain RMSD (Å)`: RMSD between side chains
- `Backbone RMSD (Å)`: RMSD between backbone

---

### Step 3: Structure Alignment and Ligand Extraction
**Script:** `scripts/utils/align_pdb_structures.py`

Align predicted structures to reference complexes and extract ligand coordinates.

**Inputs:**
- `reference_directory`: Directory containing reference/complex PDB files
- `predicted_directory`: Directory containing predicted PDB files

**Command:**
```bash
python scripts/utils/align_pdb_structures.py \
    --complexes-dir reference_directory \
    --predicted-dir predicted_directory \
    --aligned-dir aligned_predicted_directory \
    --ligand-dir aligned_ligand_directory
```

**Outputs:**
- `aligned_predicted_directory`: Aligned predicted structures
- `aligned_ligand_directory`: Extracted ligand coordinates (optional)

---

### Step 4: Ligand RMSD Calculation
**Script:** `scripts/analysis/LRMSD_calcRMS.py`

Calculate ligand RMSD and center of mass distance between experimental and predicted structures.

**Inputs:**
- `excel.xlsx`: Excel file with `Dataset_ID` and `SMILES` columns
- `reference_directory`: Directory containing reference/experimental PDB files
- `predicted_directory`: Directory containing predicted/docked PDB files

**Command:**
```bash
python scripts/analysis/LRMSD_calcRMS.py \
    --excel excel.xlsx \
    --ref-dir reference_directory \
    --pred-dir predicted_directory \
    --output output.csv
```

**Output:** CSV file with columns:
- `ID`: Dataset ID
- `CalcRMSD`: RMSD between predicted and reference ligand structures
- `COM_Distance`: Distance between centers of mass

---

## Complete Workflow Example

Here's an example of running the complete workflow:

```bash
# Step 1: Process PDB files
python scripts/utils/process_pdb_residues.py \
    --base-dir /path/to/pdb_files \
    --new-start 3

# Step 2: Calculate protein RMSD
python scripts/analysis/All_protein_RMSD.py \
    --input_dir /path/to/input_directory \
    --output_csv protein_rmsd.csv

# Step 3: Align structures and extract ligands
python scripts/utils/align_pdb_structures.py \
    --complexes-dir /path/to/reference_directory \
    --predicted-dir /path/to/predicted_directory \
    --aligned-dir /path/to/aligned_predicted_directory \
    --ligand-dir /path/to/aligned_ligand_directory

# Step 4: Calculate ligand RMSD
python scripts/analysis/LRMSD_calcRMS.py \
    --excel ligands.xlsx \
    --ref-dir /path/to/reference_directory \
    --pred-dir /path/to/predicted_directory \
    --output ligand_rmsd.csv
```

## Notes

- Steps can be run independently if you have the required inputs
- Step 3 (alignment) should be run before Step 4 (ligand RMSD) if you want to use aligned structures
- The `starting_number` parameter in Step 1 defaults to 3, but can be adjusted based on your needs
- All scripts support logging for debugging (use `--log-level` where available)


# Cofolding

Co-folding for prospective pose prediction and rescoring (Chai-1, AF3, Boltz-2)

## Repository Structure

```
Cofolding/
├── README.md                          # This file
├── TROUBLESHOOTING.md                 # Troubleshooting guide
├── config/
│   └── fold_input.json               # Template JSON for protein sequences
├── scripts/
│   ├── preprocessing/
│   │   └── input_json_generator.py   # Generate input JSON files for co-folding
│   ├── analysis/
│   │   ├── All_protein_RMSD.py       # Protein RMSD calculation using PyMOL
│   │   ├── LRMSD_calcRMS.py          # Ligand RMSD calculation using RDKit
│   │   ├── AF3_scores.py             # Parse AF3 confidence metrics
│   │   ├── Chai_scores.py            # Parse Chai-1 confidence metrics
│   │   └── hitrate.py                # Generate hit rate curves and plots
│   └── utils/
│       ├── process_pdb_residues.py   # Clean PDB files (remove residues, convert CIF to PDB)
│       ├── align_pdb_structures.py   # Align structures and extract ligands
│       ├── calculate_similarity.py   # Calculate pairwise similarity (Tanimoto or MCS)
│       ├── extract_smiles_from_lookup.py  # Extract SMILES from lookup table
│       ├── plot_correlation.py       # Plot correlation scatter plots
│       └── random_pki_generator.py   # Generate random pKi values
└── jobs/
    ├── af3-job.sh                    # SLURM job script for AlphaFold3
    ├── chai-job.sh                   # SLURM job script for Chai-1
    ├── boltz-job.sh                  # SLURM job script for Boltz-2
    └── interactions_csv.sh           # Combine IFP results
```

## Getting Started

### Setup

1. **Input configuration:**
   - Use `config/fold_input.json` as a template for protein sequence input
   - Replace the sequence with your protein-of-interest

2. **Excel file requirements:**
   - For ligand analysis scripts, Excel files must contain:
     - `Dataset_ID` column: distinguishing ID for each compound
     - `SMILES` column: SMILES strings for ligands

3. **Generate input JSON files:**
   ```bash
   python scripts/preprocessing/input_json_generator.py \
       --input_json config/fold_input.json \
       --smiles_file ligands.xlsx \
       --output_dir output_directory
   ```

4. **Run co-folding using job scripts:**
   - For Chai-1: `sbatch jobs/chai-job.sh`
   - For AF3: `sbatch jobs/af3-job.sh`
   - For Boltz-2: `sbatch jobs/boltz-job.sh`
   
   **Note:** Update paths in job scripts to match your environment before running.

## Analysis Tools

### 1. PDB File Processing

**`process_pdb_residues.py`** - Clean PDB files by removing specific residues and converting CIF to PDB format.

For Mac1, we had to remove specific residues to clean up the co-folded structure for an accurate comparison with the crystal structure. Also, the original CIF file is converted to PDB structure file, after cleaning up redundant lines. If you want to change some filters, modify the code directly.

```bash
python scripts/utils/process_pdb_residues.py \
    --base-dir base_directory \
    --new-start starting_number
```

**Arguments:**
- `--base-dir` (required): Base directory containing PDB files to process
- `--new-start` (optional, default: 3): New starting residue number for renumbering
- `--recursive`: Process PDB files in subdirectories recursively

### 2. Protein RMSD Calculation

**`All_protein_RMSD.py`** - Align the protein-ligand complex by proteins Cα and then conduct RMSD calculation for proteins, pockets, backbone, and side chains.

```bash
python scripts/analysis/All_protein_RMSD.py \
    --input_dir input_directory \
    --output_csv output.csv
```

**Arguments:**
- `--input_dir` (required): Directory containing PDB files (reference and predicted)
- `--output_csv` (required): Output CSV file for RMSD results

### 3. Structure Alignment and Ligand Extraction

**`align_pdb_structures.py`** - Align predicted to ground-truth structures and save the aligned predicted pose. Also retrieve ligand PDB file only.

```bash
python scripts/utils/align_pdb_structures.py \
    --complexes-dir reference_directory \
    --predicted-dir predicted_directory \
    --aligned-dir aligned_predicted_directory \
    --ligand-dir aligned_ligand_directory
```

**Arguments:**
- `--complexes-dir` (required): Directory containing reference/complex PDB files
- `--predicted-dir` (required): Directory containing predicted PDB files
- `--aligned-dir` (required): Output directory for aligned structures
- `--ligand-dir` (optional): Output directory for extracted ligand coordinates
- `--ligand-label` (optional, default: " LIG "): String pattern to identify ligand lines
- `--use-pymol-api` (optional): Use PyMOL Python API instead of subprocess
- `--log-level` (optional): Set logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)

### 4. Ligand RMSD Calculation

**`LRMSD_calcRMS.py`** - Ligand RMSD calculation post-alignment by protein backbone, including RDKit basic sanity checks for ligands.

The Excel file (excel.xlsx) must contain columns with ID labeled "Dataset_ID" and "SMILES".

```bash
python scripts/analysis/LRMSD_calcRMS.py \
    --excel excel.xlsx \
    --ref-dir reference_directory \
    --pred-dir predicted_directory \
    --output output.csv
```

**Arguments:**
- `--excel` (required): Path to Excel file with SMILES (must have `Dataset_ID` and `SMILES` columns)
- `--ref-dir` (required): Directory containing reference/experimental PDB files
- `--pred-dir` (required): Directory containing predicted/docked PDB files
- `--output` (required): Output CSV file path
- `--log-level` (optional): Set logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)

### 5. Chai-1 Confidence Metrics

**`Chai_scores.py`** - Parse output for Chai-1 confidence metrics.

The chaifold_directory will have folders in `*.out` format, within which there will be `scores.model_idx_#.npz` files. This code will extract `_idx_0.npz` files by default and save all scores into an Excel file.

```bash
python scripts/analysis/Chai_scores.py \
    --base-directory chaifold_directory \
    --output-file chai_scores_idx_0.xlsx
```

**Arguments:**
- `--base-directory` (required): Base directory containing .out subdirectories
- `--output-file` (required): Output Excel file path
- `--model-idx` (optional, default: 0): Model index for scores file (scores.model_idx_{model_idx}.npz)

### 6. AF3 Confidence Metrics

**`AF3_scores.py`** - Parse output for AF3 confidence metrics.

The input directory contains JSON output files from AF3.

```bash
python scripts/analysis/AF3_scores.py \
    --input_dir input_directory \
    --csv_file output.csv
```

**Arguments:**
- `--input_dir` (required): Directory containing AF3 JSON output files
- `--csv_file` (required): Output CSV file path

### 7. Pairwise Similarity Calculation

**`calculate_similarity.py`** - Calculate pairwise similarity matrix. Outputs a CSV file that contains a pairwise matrix of Tanimoto Coefficient (TC) or MCS% between molecules.

**For Tanimoto Coefficient (TC):**
```bash
python scripts/utils/calculate_similarity.py \
    --input ligands.smi \
    --output TC_matrix.csv \
    --method tanimoto
```

**For MCS%:**
```bash
python scripts/utils/calculate_similarity.py \
    --input ligands.smi \
    --output MCS_matrix.csv \
    --method mcs
```

**Arguments:**
- `--input` (required): Input .smi file (columns: SMILES ID)
- `--output` (required): Output CSV file path for similarity matrix
- `--method` (required): Similarity metric: `tanimoto` (ECFP4) or `mcs` (MCS % overlap)
- `--radius` (optional, default: 2): Morgan fingerprint radius for Tanimoto (ECFP4=2)
- `--n-bits` (optional, default: 2048): Number of fingerprint bits for Tanimoto
- `--decimals` (optional, default: 2): Number of decimal places in output
- `--log-level` (optional): Set logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)

### 8. Random pKi Value Generator

**`random_pki_generator.py`** - Generate random pKi values.

```bash
python scripts/utils/random_pki_generator.py \
    --upper upper_bound \
    --num number_of_values \
    --out output.csv
```

**Arguments:**
- `--upper` (required): Upper bound for pKi values
- `--num` (required): Number of random values to generate
- `--out` (required): Output CSV file path

### 9. Hit Rate Curves

**`hitrate.py`** - Generate hit rate curves with rolling window analysis and Wilson 95% confidence intervals.

```bash
python scripts/analysis/hitrate.py \
    --file data.xlsx \
    --sheet-name "Sheet1" \
    --score-col "L-pLDDT score" \
    --window 50
```

**Arguments:**
- `--file` (required): Path to Excel file (must have score and Active columns)
- `--sheet-name` (optional): Sheet name in Excel file (default: first sheet)
- `--score-col` (optional, default: "L-pLDDT score"): Column name for scores
- `--label-col` (optional, default: "Active"): Column name for active/inactive labels
- `--window` (optional, default: 50): Rolling window size for hit rate calculation
- `--higher-is-better` (optional, default: True): Higher scores are better
- `--lower-is-better` (optional): Lower scores are better (overrides --higher-is-better)
- `--method-name` (optional, default: "AF3"): Method name for plot legend
- `--color` (optional, default: "tab:red"): Plot color
- `--output` (optional): Output file path for plot (if not provided, displays plot)
- `--show-baseline` (optional): Show random baseline line

## Additional Utilities

### Extract SMILES from Lookup Table

**`extract_smiles_from_lookup.py`** - Extract SMILES and IDs from lookup table based on filters (e.g., extract AmpC Hits).

```bash
python scripts/utils/extract_smiles_from_lookup.py \
    --lookup Mols_labelled_hits.csv \
    --output ampc_hits.smi \
    --target AmpC \
    --hits-status Hits
```

### Plot Correlation

**`plot_correlation.py`** - Plot correlation scatter plots with color coding from Excel data.

```bash
python scripts/utils/plot_correlation.py \
    --input data.xlsx \
    --x-col "Boltz_LRMSD" \
    --y-col "DOCK_LRMSD" \
    --output correlation_plot.png
```

## Additional Tools (Referenced but not included in this repository)

These tools are mentioned but should be obtained from your lab cluster:
- `get_fingerprints_and_calc_tc_freechem.py` - Tanimoto calculation
- `ifp_interactions.py` - IFP calculation
- `get_MCS.py` - MCS calculation

For IFP of multiple protein-ligand complexes (normally we find IFP for docked ligands on the protein of interest, but for co-folding we have different protein structure for EVERY model), use:
```bash
bash jobs/interactions_csv.sh
```

## Dependencies

- Python 3.x
- RDKit (for ligand calculations and similarity)
- PyMOL (for protein RMSD calculations and alignment)
- pandas, numpy (for data processing)
- matplotlib, seaborn (for plotting)
- scipy (for statistical analysis)
- scikit-learn (for PCA and preprocessing)

## Notes

- All SLURM job scripts (`jobs/*.sh`) contain hardcoded paths specific to the original environment. Update these paths before running.
- This repository focuses on the core co-folding workflow; additional analysis tools may be available in your cluster environment.
- For troubleshooting, see `TROUBLESHOOTING.md`.

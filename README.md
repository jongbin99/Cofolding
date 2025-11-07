# Cofolding
Co-folding for prospective pose prediction and rescoring (Chai-1, AF3, Boltz-2)

## Repository Structure

```
Cofolding/
├── README.md                          # This file
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
│       └── random_pki_generator.py   # Generate random pKi values
└── jobs/
    ├── af3-job.sh                    # SLURM job script for AlphaFold3
    ├── chai-job.sh                   # SLURM job script for Chai-1
    ├── boltz-job.sh                  # SLURM job script for Boltz-2
    └── interactions_csv.sh          # Combine IFP results
```

## Getting Started

### To run co-folding:

1. **Setup input configuration:**
   - Use `config/fold_input.json` as a template for how protein sequence will be saved – just replace the sequence with protein-of-interest.
   - The excel file with SMILES must have two columns: "Dataset_ID" as a distinguishing ID, and "SMILES" for ligands.

2. **Generate input JSON files:**
   ```bash
   python scripts/preprocessing/input_json_generator.py --input_json config/fold_input.json --smiles_file **.xlsx --output_dir ../output_directory
   ```

3. **Run co-folding using job scripts:**
   - For Chai-1: `sbatch jobs/chai-job.sh`
   - For AF3: `sbatch jobs/af3-job.sh`
   - For Boltz-2: `sbatch jobs/boltz-job.sh`

   Note: You'll need to update the paths in these job scripts to match your environment before running.

## Analysis Tools

### Protein RMSD calculation
Align the protein-ligand complex by proteins Cα and then conduct RMSD calculation for proteins, pockets, backbone, side chains:
```bash
python scripts/analysis/All_protein_RMSD.py --input_dir ../input_directory --output_csv **.csv
```

### Ligand RMSD calculation
Ligand RMSD calculation post-alignment by protein backbone, including RDKit basic sanity checks for ligands:
```bash
python scripts/analysis/LRMSD_calcRMS.py --excel **.xlsx --ref-dir ../reference_directory --pred-dir ../predicted_directory --output **.csv
```

Optional arguments:
- `--log-level`: Set logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL). Default: INFO
- `--rename-output`: Output aligned ligands with matched atom names

### Parse confidence metrics

**Chai-1:**
```bash
python scripts/analysis/Chai_scores.py
```
(change directories inside the script)

**AF3:**
```bash
python scripts/analysis/AF3_scores.py --input_dir ../input_directory --csv_file **.csv
```

### Hit rate curves
```bash
python scripts/analysis/hitrate.py
```
(change directories inside the script)

### Random pKi generator
```bash
python scripts/utils/random_pki_generator.py --upper # --num ## --out output.csv
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
- RDKit (for ligand calculations)
- PyMOL (for protein RMSD calculations)
- pandas, numpy (for data processing)
- matplotlib (for plotting)

## Notes

- All SLURM job scripts (`jobs/*.sh`) contain hardcoded paths specific to the original environment. Update these paths before running.
- Analysis scripts may contain hardcoded paths that need to be updated for your environment.
- This repository focuses on the core co-folding workflow; additional analysis tools may be available in your cluster environment.

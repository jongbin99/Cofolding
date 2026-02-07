# Cofolding

Pipeline for **co-folding** protein–ligand complexes with **Chai-1**, **AlphaFold 3 (AF3)**, and **Boltz-2**, and for analyzing predicted poses (RMSD, confidence scores, similarity, hit rates).

## Overview

This repository supports:

- **Preparing inputs** (FASTA, JSON, YAML) and running co-folding jobs (SLURM scripts for Chai-1, AF3, Boltz-2).
- **Processing structures** (PDB cleanup, alignment, ligand extraction).
- **Analysis** (protein and ligand RMSD, confidence metrics, Tanimoto/MCS similarity, PCA, hit rate curves).

For step-by-step flow and data dependencies, see [WORKFLOW.md](WORKFLOW.md). For common issues, see [TROUBLESHOOTING.md](TROUBLESHOOTING.md).

## Repository structure

```
Cofolding/
├── README.md
├── WORKFLOW.md
├── TROUBLESHOOTING.md
├── config/
│   ├── fold_input.json              # Template for protein sequence input
│   └── workflow_config.json.example  # Example workflow config
├── jobs/
│   ├── af3-job.sh                   # SLURM: AlphaFold 3
│   ├── boltz-job.sh                 # SLURM: Boltz-2
│   ├── chai-job.sh                  # SLURM: Chai-1
│   └── interactions_csv.sh          # Combine IFP results
└── scripts/
    ├── run_analysis_workflow.py     # Orchestrates full analysis pipeline
    ├── preprocessing/              # Input generation and PDB prep
    │   ├── fold_input.json
    │   ├── input_fasta_generator.py
    │   ├── input_json_generator.py
    │   ├── input_yaml_generator.py
    │   ├── make_of3_input.py
    │   ├── process_pdb_residues.py   # Clean/renumber PDB; CIF→PDB
    │   └── All_protein_RMSD.py       # Protein RMSD (PyMOL)
    └── postprocessing/             # RMSD, scores, alignment, similarity, hit rates
        ├── All_protein_RMSD.py       # Protein RMSD
        ├── align_pdb_structures.py  # Align structures; extract ligands
        ├── LRMSD_calcRMS.py          # Ligand RMSD (RDKit)
        ├── AF3_scores.py             # AF3 confidence (L-pLDDT, L-PAE)
        ├── Chai_scores.py            # Chai-1 confidence
        ├── calc_mpae.py
        ├── process_pdb_residues.py
        ├── calculate_similarity.py   # Tanimoto or MCS matrix
        ├── pca_plot.py               # PCA on similarity matrix
        ├── posebusters.py
        ├── cluster_by_MCS.py
        ├── get_fingerprints_and_calc_tc_freechem.py
        ├── get_fingerprints_bfc.py
        ├── get_MCS_to_trained_set.py
        └── ifp_interactions.py      # IFP; use jobs/interactions_csv.sh to combine
```

## Quick start

### 1. Inputs and co-folding

- Put your protein sequence in `config/fold_input.json` (or use the template there).
- For ligand-based runs, use an Excel file with at least `Dataset_ID` and `SMILES`.

Generate input JSONs for co-folding:

```bash
python scripts/preprocessing/input_json_generator.py \
    --input_json config/fold_input.json \
    --smiles_file ligands.xlsx \
    --output_dir output_directory
```

Run co-folding (edit paths in the job scripts first):

```bash
sbatch jobs/chai-job.sh    # Chai-1
sbatch jobs/af3-job.sh     # AlphaFold 3
sbatch jobs/boltz-job.sh    # Boltz-2
```

### 2. Analysis workflow

The analysis pipeline: (1) process PDBs, (2) protein RMSD, (3) align structures and extract ligands, (4) ligand RMSD.

Using the orchestrator (uses a config file):

```bash
python scripts/run_analysis_workflow.py --config config/workflow_config.json
```

Running steps manually:

```bash
# Step 1: Clean and renumber PDBs
python scripts/postprocessing/process_pdb_residues.py --base-dir /path/to/pdbs --new-start 3

# Step 2: Protein RMSD
python scripts/postprocessing/All_protein_RMSD.py --input_dir /path/to/input --output_csv protein_rmsd.csv

# Step 3: Align and extract ligands
python scripts/postprocessing/align_pdb_structures.py \
    --complexes-dir /path/to/reference \
    --predicted-dir /path/to/predicted \
    --aligned-dir /path/to/aligned \
    --ligand-dir /path/to/ligands

# Step 4: Ligand RMSD
python scripts/postprocessing/LRMSD_calcRMS.py \
    --excel ligands.xlsx --ref-dir /path/to/ref --pred-dir /path/to/pred --output ligand_rmsd.csv
```

See `config/workflow_config.json.example` for a config template.

## Main scripts (by task)

### Preprocessing and PDB handling

- **`preprocessing/input_json_generator.py`** – Build co-folding input JSON from `fold_input.json` and a SMILES/Excel file.
- **`preprocessing/input_fasta_generator.py`**, **`input_yaml_generator.py`**, **`make_of3_input.py`** – Other input prep for folding.
- **`postprocessing/process_pdb_residues.py`** – Clean PDBs (remove residues, renumber), convert CIF to PDB. Options: `--base-dir`, `--new-start`, `--recursive`.

### Structure comparison

- **`postprocessing/All_protein_RMSD.py`** – Protein RMSD (full, pocket, backbone, side chain) via PyMOL. Input dir: reference + predicted PDBs. Output: CSV.
- **`postprocessing/align_pdb_structures.py`** – Align predicted to reference (Cα), write aligned PDBs and optional ligand-only PDBs.
- **`postprocessing/LRMSD_calcRMS.py`** – Ligand RMSD and COM distance (RDKit). Requires Excel with `Dataset_ID` and `SMILES`; ref and pred PDB dirs.

### Confidence and scores

- **`postprocessing/AF3_scores.py`** – Parse AF3 JSON outputs → CSV (e.g. L-pLDDT, L-PAE).
- **`postprocessing/Chai_scores.py`** – Parse Chai-1 output (e.g. `scores.model_idx_*.npz`) → Excel.
- **`postprocessing/calc_mpae.py`** – mPAE calculation (for AF3)

### Similarity and clustering

- **`postprocessing/calculate_similarity.py`** – Pairwise similarity matrix (Tanimoto or MCS%) from a `.smi` file. Options: `--method tanimoto|mcs`, `--radius`, `--n-bits`, etc.
- **`postprocessing/pca_plot.py`** – PCA on a similarity matrix CSV; saves scatter plot and PCA results CSV. No lookup/labels.
- **`postprocessing/cluster_by_MCS.py`**, **`get_MCS_to_trained_set.py`**, **`get_fingerprints_*.py`** – MCS and fingerprint-based analysis.

### Other

- **`postprocessing/posebusters.py`** – PoseBusters checks with following lines under gimel:
source /nfs/home/jkim/miniconda3/bin/activate
conda activate py3.10
pip show posebusters
bust ***_predicted.sdf -l ***_reference.sdf protein.pdb --outfmt long


## Dependencies

- Python 3.x  
- RDKit (ligand handling, similarity)  
- PyMOL (protein RMSD, alignment)  
- pandas, numpy, matplotlib, seaborn, scipy, scikit-learn  

## Notes

- All workflow and analysis scripts live under **`scripts/preprocessing/`** and **`scripts/postprocessing/`** (there is no `scripts/analysis/` or `scripts/utils/`).
- **`jobs/`** (SLURM scripts for Chai-1, AF3, Boltz-2) is optional; if present, update paths inside the scripts for your cluster before running.
- **Excel inputs** for ligand workflows must include `Dataset_ID` and `SMILES` where required.
- For the full pipeline and commands, see [WORKFLOW.md](WORKFLOW.md). For common issues, see [TROUBLESHOOTING.md](TROUBLESHOOTING.md).

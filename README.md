# Cofolding

Pipeline for **co-folding** protein–ligand complexes (Chai-1, AlphaFold 3, Boltz-2) and for **analyzing** predicted poses: RMSD, confidence scores, similarity, and PCA.
The ReadMe is organized using CursorAI

## Folder layout

```
Cofolding/
├── README.md
├── WORKFLOW.md
├── TROUBLESHOOTING.md
├── .gitignore
├── config/
│   ├── fold_input.json
│   └── workflow_config.json.example
└── scripts/
    ├── run_analysis_workflow.py
    ├── preprocessing/
    │   ├── fold_input.json
    │   ├── input_fasta_generator.py
    │   ├── input_json_generator.py
    │   ├── input_yaml_generator.py
    │   ├── make_of3_input.py
    │   ├── process_pdb_residues.py
    │   └── All_protein_RMSD.py
    └── postprocessing/
        ├── process_pdb_residues.py
        ├── All_protein_RMSD.py
        ├── align_pdb_structures.py
        ├── LRMSD_calcRMS.py
        ├── AF3_scores.py
        ├── Chai_scores.py
        ├── calc_mpae.py
        ├── calculate_similarity.py
        ├── pca_plot.py
        ├── cluster_by_MCS.py
        ├── get_fingerprints_and_calc_tc_freechem.py
        ├── get_fingerprints_bfc.py
        ├── get_MCS_to_trained_set.py
        └── posebusters.py
```

- **config/** — Input templates and workflow config example.
- **scripts/preprocessing/** — Input generation (FASTA, JSON, YAML) and PDB prep; used before running co-folding.
- **scripts/postprocessing/** — All analysis after co-folding: PDB cleanup, RMSD, alignment, scores, similarity, PCA, clustering, PoseBusters.
- **scripts/run_analysis_workflow.py** — Runs the four-step analysis pipeline (see [WORKFLOW.md](WORKFLOW.md)).

## Quick start

### Inputs and co-folding

1. Set your protein sequence in **config/fold_input.json** (or **scripts/preprocessing/fold_input.json**).
2. For ligand-based runs, use an Excel file with columns **Dataset_ID** and **SMILES**.

Generate input JSONs for co-folding:

```bash
python scripts/preprocessing/input_json_generator.py \
  --input_json config/fold_input.json \
  --smiles_file ligands.xlsx \
  --output_dir output_directory
```

Run your co-folding jobs (Chai-1, AF3, or Boltz-2) with your own cluster/SLURM setup as needed.

### Analysis workflow (four steps)

From the repo root. Either use the orchestrator:

```bash
python scripts/run_analysis_workflow.py --config config/workflow_config.json
```

Or run steps manually (all scripts in **scripts/postprocessing/**):

```bash
# 1. Clean and renumber PDBs
python scripts/postprocessing/process_pdb_residues.py --base-dir /path/to/pdbs --new-start 3

# 2. Protein RMSD
python scripts/postprocessing/All_protein_RMSD.py --input_dir /path/to/input --output_csv protein_rmsd.csv

# 3. Align structures and extract ligands
python scripts/postprocessing/align_pdb_structures.py \
  --complexes-dir /path/to/reference \
  --predicted-dir /path/to/predicted \
  --aligned-dir /path/to/aligned \
  --ligand-dir /path/to/ligands

# 4. Ligand RMSD
python scripts/postprocessing/LRMSD_calcRMS.py \
  --excel ligands.xlsx --ref-dir /path/to/ref --pred-dir /path/to/pred --output ligand_rmsd.csv
```

Use **config/workflow_config.json.example** as a template for the orchestrator.

---

## Scripts by folder

### config/

| File | Purpose |
|------|--------|
| **fold_input.json** | Template for protein sequence input. |
| **workflow_config.json.example** | Example config for `run_analysis_workflow.py` (paths and step options). |

### scripts/preprocessing/

| Script | Purpose |
|--------|--------|
| **input_json_generator.py** | Build co-folding input JSON from a fold template and SMILES/Excel file. |
| **input_fasta_generator.py** | Generate FASTA inputs. |
| **input_yaml_generator.py** | Generate YAML inputs. |
| **make_of3_input.py** | Prepare input for OF3/co-folding. |
| **process_pdb_residues.py** | Clean and renumber PDBs; CIF → PDB. |
| **All_protein_RMSD.py** | Protein RMSD (PyMOL). |

### scripts/postprocessing/

| Script | Purpose |
|--------|--------|
| **process_pdb_residues.py** | Clean/renumber PDBs; CIF → PDB (workflow step 1). |
| **All_protein_RMSD.py** | Protein RMSD — full, pocket, backbone, side chain (workflow step 2). |
| **align_pdb_structures.py** | Align predicted to reference; extract ligands (workflow step 3). |
| **LRMSD_calcRMS.py** | Ligand RMSD and COM distance; needs Excel with Dataset_ID and SMILES (workflow step 4). |
| **AF3_scores.py** | Parse AF3 JSON → CSV (e.g. L-pLDDT, L-PAE). |
| **Chai_scores.py** | Parse Chai-1 output (e.g. scores.model_idx_*.npz) → Excel. |
| **calc_mpae.py** | mPAE calculation. |
| **calculate_similarity.py** | Pairwise similarity matrix (Tanimoto or MCS%) from a .smi file. |
| **pca_plot.py** | PCA on a similarity matrix CSV; plot and results CSV. |
| **cluster_by_MCS.py** | Clustering by MCS. |
| **get_fingerprints_and_calc_tc_freechem.py** | Fingerprints and Tanimoto (RDKit/FreeChem). |
| **get_fingerprints_bfc.py** | Fingerprints (BFc). |
| **get_MCS_to_trained_set.py** | MCS to trained set. |
| **posebusters.py** | PoseBusters checks (e.g. `bust pred.sdf -l ref.sdf protein.pdb --outfmt long`). |

Run any script with `--help` for arguments.

---

## Dependencies

- Python 3.x  
- **RDKit** — ligand handling, similarity, fingerprints  
- **PyMOL** — protein RMSD and alignment  
- pandas, numpy, matplotlib, seaborn, scipy, scikit-learn  

PoseBusters scripts need a separate PoseBusters environment (e.g. `conda activate <env>` and `pip show posebusters`).

---

## More info

- **Full pipeline and commands:** [WORKFLOW.md](WORKFLOW.md)  
- **Common issues:** [TROUBLESHOOTING.md](TROUBLESHOOTING.md)  
- **Excel:** Ligand workflows expect **Dataset_ID** and **SMILES** where noted.

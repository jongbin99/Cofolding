# Troubleshooting Guide

## Common Issues and Solutions

### RDKit Import Error

**Error Message:**
```
ImportError: No module named rdkit.Chem
```

**Solution:** Use a Python environment with RDKit installed:

```bash
# Activate your conda environment (e.g. py3.10)
conda activate py3.10

# Verify RDKit is available
python -c "import rdkit; print('RDKit version:', rdkit.__version__)"

# Run scripts from the repo root; workflow scripts live in scripts/postprocessing/
python scripts/postprocessing/get_fingerprints_bfc.py mac1.smi outname 0.35 557
```

**Alternative:** Install RDKit in your current environment:
```bash
conda install -c conda-forge rdkit
# OR
pip install rdkit-pypi
```

### Environment Setup Issues

**Problem:** Different scripts require different Python environments (e.g. RDKit, PyMOL, PoseBusters).

**Solution:** Use a single environment with all dependencies, or document which script needs which env:

- **Ligand/similarity (RDKit):** `scripts/postprocessing/calculate_similarity.py`, `LRMSD_calcRMS.py`, `get_fingerprints_*.py`, etc.
- **Protein RMSD / alignment (PyMOL):** `scripts/postprocessing/All_protein_RMSD.py`, `align_pdb_structures.py`
- **PoseBusters:** `scripts/postprocessing/posebusters.py` — requires PoseBusters and its env (e.g. `conda activate py3.10` and `pip show posebusters`)

Create a small wrapper if needed:
```bash
#!/bin/bash
# run_from_repo_root.sh
cd /path/to/Cofolding
conda activate py3.10
python scripts/postprocessing/your_script.py "$@"
```

### Path Not Found Errors

**Problem:** Scripts or job scripts reference paths that don't exist.

**Solution:**
1. Run workflow scripts from the **Cofolding repository root** so paths like `scripts/postprocessing/...` resolve.
2. Use absolute paths in config files (e.g. `config/workflow_config.json`) for base dirs and output CSV paths.
3. If you use **`jobs/`** (SLURM scripts for Chai-1, AF3, Boltz-2), edit the paths inside those scripts to match your cluster and project directories before submitting.

### Workflow Script Not Found

**Problem:** `run_analysis_workflow.py` or a step fails with "script not found".

**Solution:** All step scripts live under **`scripts/postprocessing/`**:
- Step 1: `scripts/postprocessing/process_pdb_residues.py`
- Step 2: `scripts/postprocessing/All_protein_RMSD.py`
- Step 3: `scripts/postprocessing/align_pdb_structures.py`
- Step 4: `scripts/postprocessing/LRMSD_calcRMS.py`

Ensure you have not removed or moved the `postprocessing/` folder. See [WORKFLOW.md](WORKFLOW.md) for full commands.

## Verification Checklist

Before running co-folding or analysis:

- [ ] RDKit is installed and importable (for ligand/similarity scripts)
- [ ] PyMOL is available (for protein RMSD and alignment)
- [ ] You are in or reference the Cofolding repo root so `scripts/postprocessing/` paths work
- [ ] Required inputs exist (SMILES, PDBs, Excel with `Dataset_ID`/`SMILES` where needed)
- [ ] Output directories are writable
- [ ] If using SLURM, job scripts in `jobs/` have paths updated for your environment
- [ ] Correct conda/environment is activated

## Getting Help

1. See [README.md](README.md) for setup and script overview.
2. See [WORKFLOW.md](WORKFLOW.md) for the four-step analysis pipeline and exact commands.
3. Run any script with `--help` for arguments (e.g. `python scripts/postprocessing/LRMSD_calcRMS.py --help`).
4. Check file permissions and that paths in your config file are absolute and valid.



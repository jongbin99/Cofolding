# Troubleshooting Guide

## Common Issues and Solutions

### RDKit Import Error

**Error Message:**
```
ImportError: No module named rdkit.Chem
```

**Solution:**
Ensure you're using the correct Python environment with RDKit installed:

```bash
# Activate the py3.10 environment (common on HPC clusters)
source /nfs/soft/anaconda3/bin/activate
conda activate py3.10

# Verify RDKit is available
python -c "import rdkit; print('RDKit version:', rdkit.__version__)"

# Then run your script
python ~ttummino/zzz.scripts/tc_and_clustering/get_fingerprints_bfc.py mac1.smi outname 0.35 557
```

**Alternative Solutions:**
```bash
# Install RDKit in current environment
conda install -c conda-forge rdkit
# OR
pip install rdkit-pypi
```

### Environment Setup Issues

**Problem:** Different scripts require different Python environments

**Solution:** Create a wrapper script or document required environments:

```bash
#!/bin/bash
# activate_py310.sh
source /nfs/soft/anaconda3/bin/activate
conda activate py3.10
export PATH="/nfs/home/.conda/envs/py3.10/bin:$PATH"
```

### Path Not Found Errors

**Problem:** Scripts reference paths that don't exist in your environment

**Solution:** 
1. Check the hardcoded paths in the script
2. Update paths to match your cluster/filesystem structure
3. For co-folding scripts in `jobs/`, update the directory paths before running

## Verification Checklist

Before running co-folding or analysis scripts:

- [ ] RDKit is installed and importable
- [ ] PyMOL is available (for protein RMSD calculations)
- [ ] Required data files exist (SMILES, PDB files, etc.)
- [ ] Output directories are writable
- [ ] SLURM job scripts have correct paths updated
- [ ] Correct conda environment is activated

## Getting Help

For additional issues:
1. Check the main README.md for usage instructions
2. Verify all dependencies are installed in your environment
3. Check file permissions and path accessibility
4. Review the script's command-line arguments with `--help`



# Cofolding
Co-folding for prospective pose prediction and rescoring (Chai-1, AF3, Boltz-2)

To run co-folding:
- Use fold_input.json as a template for how protein sequence will be saved – just replace the sequence with protein-of-interest.
- The excel file with SMILES must have two columns: “Dataset_ID” as a distinguishing ID, and “SMILES” for ligands.
- Automate input json file generation (can design something similar in .fasta for Chai-1 and .yaml for Boltz-2): input_json_generator.py
"python input_json_generator.py --input_json **.json --smiles_file **.xlsx --output_dir ../output_directory"

Then use any of the job scripts below to co-fold:
- Slurm job script for chai: chai-job.sh
- Slurm job script for AF3: af3-job.sh
- Slurm job script for Boltz-2: boltz-job.sh

To analyze co-folding
- Align the protein-ligand complex by proteins Cα and then conduct RMSD calculation for proteins, pockets, backbone, side chains: All_protein_RMSD.py
"python All_protein_RMSD.py --input_dir ../input_directory --output_csv **.csv"

- Ligand RMSD calculation post-alignment by protein backbone, including RDkit basic sanity checks for ligands: LRMSD_calcRMSD.py
"python LRMSD_calcRMS.py
--excel **.xlsx
--ref-dir ../reference_directory
--pred-dir ../predicted_directory
--output **.csv"

- Parse output for Chai-1 confidence metrics: Chai_scores.py
"python Chai_scores.py" (change directories inside the script)

- Parse output for AF3 confidence metrics: AF3_scores.py
"python AF3_scores.py --input_dir ../input_directory --csv_file **.csv"

- Tanimoto calculation: get_fingerprints_and_calc_tc_freechem.py (used the pre-coded Tanimoto calculation from the lab cluster)

- IFP calculation: ifp_interactions.py (used the pre-coded IFP script from the lab cluster)

- MCS calculation: get_MCS.py (input query smiles and trained set xlsx)
"python get_MCS.py **.smi **.xlsx "

- Random pKi value generator: random_pki_generator.py
"python random_pki_generator.py --upper # --num ## --out output.csv"

- Hit rate curves: hitrate.py
"python hitrate.py" (change directories inside the script)

- Slurm job script for IFP of multiple protein-ligand complexes (normally we find IFP for docked ligands on the protein of interest, but for co-folding we have different protein structure for EVERY model, so a separate job script was needed): interactions_csv.sh

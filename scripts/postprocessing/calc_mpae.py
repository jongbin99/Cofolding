import json
import os
import csv
from unittest import result
import numpy as np
import pandas as pd

def ligand_mpae(input_dir, output_dir, output_csv):
    os.makedirs(output_dir)
    
    results = []

    for filename in os.listdir(input_dir):
        if filename.endswith(".json"):
            filepath = os.path.join(input_dir, filename)

            try:
                with open(filepath, 'r') as f:
                    af_data = json.load(f)

                pae_matrix = np.array(af_data["pae"])
                res_ids = af_data["token_res_ids"]       # List of ints
                chain_ids = af_data["token_chain_ids"]   # List of "A", "B", etc.

                if not (len(res_ids) == len(chain_ids) == pae_matrix.shape[0]):
                    raise ValueError("Mismatch in dimensions of residue IDs, chain IDs, or PAE matrix")

                protein_indices = [i for i, chain in enumerate(chain_ids) if chain == "A"]
                ligand_indices = [i for i, chain in enumerate(chain_ids) if chain == "B"]

                if not protein_indices or not ligand_indices:
                    print(f"Skipping {filename}: Missing chain A or B")
                    continue

                # Find minimum PAE value between protein (A) and ligand (B)
                min_pair = min(((i, j) for i in protein_indices for j in ligand_indices),
                            key=lambda pair: pae_matrix[pair[0], pair[1]])
                min_pae_value = pae_matrix[min_pair[0], min_pair[1]]

                results.append({
                    "file": filename,
                    "min_PAE": float(min_pae_value),
                    "protein_index": int(min_pair[0]),
                    "ligand_index": int(min_pair[1]),
                    "protein_residue": res_ids[min_pair[0]],
                    "ligand_residue": res_ids[min_pair[1]]
                })
            except Exception as e:
                print(f"Error processing {filename}: {e}")
                continue

# Save to CSV
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)

# Define your directories here
# -----------------------------
input_dir  = "/nfs/home/jkim/work/projects/rotation/250115_alphafold3/conf_json"       # folder with input JSON files (confidence.json files)
output_dir = "/nfs/home/jkim/work/projects/rotation/250115_alphafold3/conf_json/conf_json_lig_mpae"  # an arbitrary output folder that contains JSON fil$
output_csv   = "/nfs/home/jkim/work/projects/rotation/250115_alphafold3/conf_json/mpae.csv"	 # CSV file summarizing average pLDDT

# Run the function
ligand_mpae(input_dir, output_dir, output_csv)
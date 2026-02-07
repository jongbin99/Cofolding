#!/usr/bin/env python3
import os
import copy
import argparse
import pandas as pd
import yaml

def load_table(path):
    if path.lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    return pd.read_csv(path)

def process_smiles(input_yaml, smiles_file, output_dir):
    with open(input_yaml, "r") as f:
        base_yaml = yaml.safe_load(f)

    df = load_table(smiles_file)
    os.makedirs(output_dir, exist_ok=True)

    for _, row in df.iterrows():
        smiles = row["SMILES"]
        dataset_id = row["Dataset_ID"]

        updated = copy.deepcopy(base_yaml)

        # name
        updated["name"] = dataset_id

        # sequences: keep protein as-is from template, append ligand
        if "sequences" not in updated:
            updated["sequences"] = []
        updated["sequences"].append({
            "ligand": {
                "id": ["B"],
                "smiles": smiles
            }
        })

        # complexes: define protein A with ligand B
        updated["complexes"] = [{
            "id": f"{dataset_id}_comp",
            "protein": ["A"],
            "ligand":  ["B"],
        }]

        # properties: affinity binder B (as in your example)
        updated["properties"] = [{
            "affinity": {
                "binder": "B"
            }
        }]

        out_file = os.path.join(output_dir, f"{dataset_id}.yaml")
        with open(out_file, "w") as out:
            yaml.safe_dump(updated, out, sort_keys=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate yaml input files for Boltz-2 co-folding")
    parser.add_argument("--input_yaml", required=True,
                        help="Template YAML containing only the protein (id [A]) and optionally other fixed fields.")
    parser.add_argument("--smiles_file", required=True,
                        help="Excel/CSV with columns 'SMILES' and 'Dataset_ID'.")
    parser.add_argument("--output_dir", required=True,
                        help="Directory where new YAML files will be saved.")

    args = parser.parse_args()
    process_smiles(args.input_yaml, args.smiles_file, args.output_dir)
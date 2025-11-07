import json
import os
import pandas as pd
import copy
import argparse
#input_json needs to be a template json that contains the sequence for the protein 'only' (no ligands) - I have included an example named fold_input.json
#smiles_file is an excel with two columns: 'SMILES' and 'Dataset_ID' (can be any name, just needs to match the code below)
#output_dir is the directory where you want to save the new json files

def process_smiles(input_json, smiles_file, output_dir):

    # Load the input JSON file
    with open(input_json, 'r') as infile:
        base_json = json.load(infile)

    df = pd.read_excel(smiles_file)  # Use read_csv if it's a CSV file

    # Iterate through the SMILES and create a JSON file for each
    for index, row in df.iterrows():
        smiles = row['SMILES']
        dataset_id = row['Dataset_ID'] 
        new_entry = {
            "smiles": smiles,
            "id": ["B"]
        }

        # Create a new JSON file with the updated data
        updated_json = copy.deepcopy(base_json)  # Deep copy for independence
        new_json = {"name": dataset_id}
        new_json.update(updated_json)
        updated_json = new_json

        if "sequences" not in updated_json:
            updated_json["sequences"] = [] 
        updated_json["sequences"].append({"ligand": new_entry})

        # Save the updated JSON to a new file
        output_file = os.path.join(output_dir, f"{dataset_id}.json")
        with open(output_file, 'w') as outfile:
            json.dump(updated_json, outfile, indent=4)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate per-ligand JSONs"
    )
    parser.add_argument("--input_json", required=True,
                        help="Template JSON containing only the protein sequence.")
    parser.add_argument("--smiles_file", required=True,
                        help="Excel or CSV file with 'SMILES' and 'Dataset_ID' columns.")
    parser.add_argument("--output_dir", required=True,
                        help="Directory where new JSON files will be saved.")

    args = parser.parse_args()

    process_smiles(args.input_json, args.smiles_file, args.output_dir)
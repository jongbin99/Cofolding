import os
import glob
import json
import csv
import argparse

def af3_score(input_dir, csv_file):
    results = []

    #loop through all the json files in the input directory
    for json_file in glob.glob(os.path.join(input_dir, "*.json")):
        with open(json_file, "r") as file:
            data = json.load(file)
            
        #get the atom_chain_ids, atom_plddts, token_chain_ids, and atom_pae from the json file
        atom_chain_ids = data.get("atom_chain_ids", [])
        atom_plddts   = data.get("atom_plddts", [])
        token_chain_ids = data.get("token_chain_ids", [])
        atom_pae = data.get("pae", [])
        ref_prefix    = os.path.splitext(os.path.basename(json_file))[0]

        filtered_plddts = [plddt for chain_id, plddt in zip(atom_chain_ids, atom_plddts) if chain_id == "B"]
        average_plddt = sum(filtered_plddts) / len(filtered_plddts) if filtered_plddts else 0
        filtered_pae = [pae for chain_id, pae in zip(token_chain_ids, atom_pae) if chain_id == "B"]
        flattened_pae = [value for row in filtered_pae for value in row]
        average_pae = sum(flattened_pae) / len(flattened_pae) if flattened_pae else 0
        results.append([ref_prefix, average_plddt, average_pae])

    with open(csv_file, "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Ligand", "L-pLDDT", "L-PAE"])
        writer.writerows(results)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Summarize AlphaFold3 ligand confidence (L-pLDDT, L-PAE) across JSON outputs."
    )
    parser.add_argument("--input_dir", required=True,
                        help="Directory containing AlphaFold3 output JSON files.")
    parser.add_argument("--csv_file", required=True,
                        help="Output CSV file path for summary results.")

    args = parser.parse_args()
    af3_score(args.input_dir, args.csv_file)

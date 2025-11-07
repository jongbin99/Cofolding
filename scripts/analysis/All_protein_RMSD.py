import os
import csv
from pymol import cmd
import argparse

def rmsd_calculation(input_dir, output_csv):
    reference_files = []
    predicted_files = []
    pocket_residues = "21+22+23+48+49+52+125+126+129+130+154+155+156+157+160" #residues in the binding pocket

    # Get reference and predicted files, when they are named so that reference files are in 'mac-x...' format and predicted files contain '_pred_chainA'
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".pdb"):
            if "_pred_chainA.pdb" in file_name:
                predicted_files.append(file_name)
            elif file_name.startswith("mac-x"):
                reference_files.append(file_name)

    rmsd_results = []

    # Iterate over reference files
    for ref_file in reference_files:
        prefix = os.path.splitext(ref_file)[0]
        pred_file = next((f for f in predicted_files if f.startswith(prefix) and "_pred_chainA" in f), None)

        if pred_file is None:
            print(f"[WARN] No predicted file found for {ref_file}; skipping.")
            continue

        ref_path = os.path.join(input_dir, ref_file)
        pred_path = os.path.join(input_dir, pred_file)

        ref_object = f"{prefix}_R"
        pred_object = f"{prefix}_P"
    
        try:
            cmd.load(ref_path, ref_object)
            cmd.load(pred_path, pred_object)
            print(f"[INFO] Loaded {ref_object} and {pred_object}")
        except Exception as e:
                print(f"[ERROR] Failed to load structures: {e}")
                continue

        ref_count = cmd.count_atoms(ref_object)
        pred_count = cmd.count_atoms(pred_object)

#Align and compute total RMSD
        try:
            alignment_result = cmd.align(pred_object, ref_object)
            full_rmsd = alignment_result[0]
        except Exception as e:
            print(f"[ERROR] Alignment failed: {e}")
            continue

#Compute pocket RMSD
        try:
            pocket_ref = f"{ref_object} and resi {pocket_residues}"
            pocket_pred = f"{pred_object} and resi {pocket_residues}"
            pocket_rmsd = cmd.rms_cur(pocket_pred, pocket_ref)

        except Exception as e:
            pocket_rmsd = "N/A"
    
#Compute pocket side-chain RMSD only
        try:
            sc_ref = f"{ref_object} and resi {pocket_residues} and not name CA+C+N+O"
            sc_pred = f"{pred_object} and resi {pocket_residues} and not name CA+C+N+O"
            side_chain_rmsd = cmd.rms_cur(sc_pred, sc_ref)

        except Exception as e:
            side_chain_rmsd = "N/A"
        
#Compute pocket backbone RMSD only
        try:
            backbone_selection_ref = f"{ref_object} and resi {pocket_residues} and name CA+C+N+O"
            backbone_selection_pred = f"{pred_object} and resi {pocket_residues} and name CA+C+N+O"
            ref_count_bb = cmd.count_atoms(backbone_selection_ref)
            pred_count_bb = cmd.count_atoms(backbone_selection_pred)

            if ref_count_bb == 0 or pred_count_bb == 0:
                backbone_rmsd = "N/A"
            else:
                backbone_rmsd = cmd.rms_cur(backbone_selection_pred, backbone_selection_ref)
        except Exception as e:
            backbone_rmsd = "N/A"

        rmsd_results.append([prefix, full_rmsd, pocket_rmsd, side_chain_rmsd, backbone_rmsd])
        cmd.delete(ref_object)
        cmd.delete(pred_object)

# Save results to CSV
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Protein", "Full Structure RMSD (Å)", "Pocket RMSD (Å)", "Side-Chain RMSD (Å)", "Backbone RMSD (Å)"])
        writer.writerows(rmsd_results)

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Calculate RMSD between reference and predicted protein structures.")
    parser.add_argument("--input_dir", required=True, help="Directory containing PDB files.")
    parser.add_argument("--output_csv", required=True, help="Output CSV file for RMSD results.")
    args=parser.parse_args()
    rmsd_calculation(args.input_dir, args.output_csv)
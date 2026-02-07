import os
import csv
from pymol import cmd
import argparse

def rmsd_calculation(input_dir, output_csv, reference_file="5SQW.pdb"):
    predicted_files = []
    pocket_residues = "130"
    #pocket_residues = "21+22+23+48+49+52+125+126+129+130+154+155+156+157+160" #residues in the binding pocket

    # Determine reference file path
    if os.path.isabs(reference_file):
        ref_path = reference_file
    else:
        ref_path = os.path.join(input_dir, reference_file)
    
    if not os.path.exists(ref_path):
        print(f"[ERROR] Reference file not found: {ref_path}")
        return
    
    ref_file_name = os.path.basename(reference_file)

    # Get predicted files (exclude the reference file)
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".pdb") and file_name != ref_file_name:
            predicted_files.append(file_name)

    print(f"[INFO] Using reference file: {ref_file_name}")
    print(f"[INFO] Found {len(predicted_files)} predicted file(s) to compare")

    # Load reference structure once
    ref_object = "reference"
    try:
        cmd.load(ref_path, ref_object)
        print(f"[INFO] Loaded reference structure: {ref_object}")
    except Exception as e:
        print(f"[ERROR] Failed to load reference structure: {e}")
        return

    rmsd_results = []

    # Iterate over all predicted files and compare to single reference
    for pred_file in predicted_files:
        pred_path = os.path.join(input_dir, pred_file)
        pred_name = os.path.splitext(pred_file)[0]
        pred_object = f"pred_{pred_name}"
    
        try:
            cmd.load(pred_path, pred_object)
            print(f"[INFO] Loaded predicted structure: {pred_object}")
        except Exception as e:
            print(f"[ERROR] Failed to load {pred_file}: {e}")
            continue

        # Align and compute total RMSD
        try:
            alignment_result = cmd.align(pred_object, ref_object)
            full_rmsd = alignment_result[0]
        except Exception as e:
            print(f"[ERROR] Alignment failed for {pred_file}: {e}")
            cmd.delete(pred_object)
            continue

        # Compute pocket RMSD
        try:
            pocket_ref = f"{ref_object} and resi {pocket_residues}"
            pocket_pred = f"{pred_object} and resi {pocket_residues}"
            pocket_rmsd = cmd.rms_cur(pocket_pred, pocket_ref)
        except Exception as e:
            pocket_rmsd = "N/A"
    
        # Compute pocket side-chain RMSD only
        try:
            sc_ref = f"{ref_object} and resi {pocket_residues} and not name CA+C+N+O"
            sc_pred = f"{pred_object} and resi {pocket_residues} and not name CA+C+N+O"
            side_chain_rmsd = cmd.rms_cur(sc_pred, sc_ref)
        except Exception as e:
            side_chain_rmsd = "N/A"
        
        # Compute pocket backbone RMSD only
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

        rmsd_results.append([pred_name, full_rmsd, pocket_rmsd, side_chain_rmsd, backbone_rmsd])
        cmd.delete(pred_object)

    # Clean up reference object
    cmd.delete(ref_object)

    # Save results to CSV
    with open(output_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Predicted File", "Full Structure RMSD (Å)", "Pocket RMSD (Å)", "Side-Chain RMSD (Å)", "Backbone RMSD (Å)"])
        writer.writerows(rmsd_results)
    
    print(f"[INFO] Results saved to {output_csv}")

if __name__ == "__main__":
    parser=argparse.ArgumentParser(description="Calculate RMSD between a single reference structure and multiple predicted protein structures.")
    parser.add_argument("--input_dir", required=True, help="Directory containing PDB files.")
    parser.add_argument("--output_csv", required=True, help="Output CSV file for RMSD results.")
    parser.add_argument("--reference_file", default="5SQW.pdb", help="Reference PDB file name (default: 5SQW.pdb). Can be filename or full path.")
    args=parser.parse_args()
    rmsd_calculation(args.input_dir, args.output_csv, args.reference_file)
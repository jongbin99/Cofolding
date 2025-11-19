
import argparse
import glob
import os
import sys
import numpy as np
import pandas as pd


def extract_chai_scores(base_directory, output_file, model_idx=0):
    # Validate base directory
    if not os.path.isdir(base_directory):
        print(f"[ERROR] Base directory not found: {base_directory}", file=sys.stderr)
        return 1
    
    # Find all .out directories
    out_folders = [
        folder for folder in glob.glob(os.path.join(base_directory, "*.out"))
        if os.path.isdir(folder)
    ]
    
    if len(out_folders) == 0:
        print(f"[WARN] No .out directories found in {base_directory}", file=sys.stderr)
        return 1
    
    print(f"[INFO] Found {len(out_folders)} .out directories")
    
    data_list = []
    processed = 0
    failed = 0
    
    for out_folder in out_folders:
        npz_file_path = os.path.join(
            out_folder, f"scores.model_idx_{model_idx}.npz"
        )
        
        if not os.path.exists(npz_file_path):
            print(f"[WARN] Scores file not found: {npz_file_path}", file=sys.stderr)
            failed += 1
            continue
        
        try:
            data = np.load(npz_file_path)
            
            aggregate_score = data.get('aggregate_score', np.nan)
            iptm = data.get('iptm', np.nan)
            ptm = data.get('ptm', np.nan)
            per_chain_ptm = data.get('per_chain_ptm', np.nan)
            per_chain_pair_iptm = data.get('per_chain_pair_iptm', np.nan)
            has_inter_chain_clashes = data.get('has_inter_chain_clashes', np.nan)
            chain_chain_clashes = data.get('chain_chain_clashes', np.nan)
            
            data_list.append({
                "Index": os.path.basename(out_folder),
                "Aggregate Score": aggregate_score,
                "iPTM": iptm,
                "pTM": ptm,
                "pTM per chain": per_chain_ptm,
                "iPTM per chain pair": per_chain_pair_iptm,
                "Interchain clash?": has_inter_chain_clashes,
                "Matrix chain-chain clashes": chain_chain_clashes
            })
            processed += 1
            
        except Exception as e:
            print(f"[ERROR] Failed to process {out_folder}: {e}", file=sys.stderr)
            failed += 1
            continue
    
    if len(data_list) == 0:
        print(f"[ERROR] No valid scores found to write", file=sys.stderr)
        return 1
    
    # Create DataFrame and save to Excel
    try:
        chai_output_scores = pd.DataFrame(data_list)
        chai_output_scores.to_excel(output_file, index=False)
        print(f"[INFO] Successfully processed {processed} directories")
        if failed > 0:
            print(f"[WARN] Failed to process {failed} directories")
        print(f"[INFO] Output written to: {output_file}")
    except Exception as e:
        print(f"[ERROR] Failed to write output file: {e}", file=sys.stderr)
        return 1
    
    return 0


def main():
    parser = argparse.ArgumentParser(
        description='Extract Chai-1 scores from .npz files in .out directories',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '--base-directory', required=True,
        help='Base directory containing .out subdirectories'
    )
    parser.add_argument(
        '--output-file', required=True,
        help='Output Excel file path'
    )
    parser.add_argument(
        '--model-idx', type=int, default=0,
        help='Model index for scores file (scores.model_idx_{model_idx}.npz)'
    )
    
    args = parser.parse_args()
    
    return extract_chai_scores(
        base_directory=args.base_directory,
        output_file=args.output_file,
        model_idx=args.model_idx
    )


if __name__ == "__main__":
    sys.exit(main())
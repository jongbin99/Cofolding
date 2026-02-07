#!/usr/bin/env python3
import os
import argparse
import pandas as pd

def read_protein_sequence(fasta_path):
    seq = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq:
                    break
                continue
            seq.append(line)
    return "".join(seq)

def load_smiles_table(path):
    if path.lower().endswith((".xlsx", ".xls")):
        return pd.read_excel(path)
    else:
        return pd.read_csv(path)

def process_smiles(template_protein_fasta, smiles_file, output_dir, protein_name="mac1"):
    protein_seq = read_protein_sequence(template_protein_fasta)
    df = load_smiles_table(smiles_file)

    os.makedirs(output_dir, exist_ok=True)

    for _, row in df.iterrows():
        smiles = row["SMILES"]
        dataset_id = row["Dataset_ID"]

        out_file = os.path.join(output_dir, f"{dataset_id}.fasta")
        with open(out_file, "w") as out:
            out.write(f">protein|name={protein_name}\n")
            out.write(protein_seq + "\n")
            out.write(f">ligand|name={dataset_id}\n")
            out.write(smiles + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate FASTA inputs for Chai-1 co-folding")
    parser.add_argument("--template_protein_fasta", required=True)
    parser.add_argument("--smiles_file", required=True)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--protein_name", default="mac1")

    args = parser.parse_args()

    process_smiles(
        args.template_protein_fasta,
        args.smiles_file,
        args.output_dir,
        args.protein_name,
    )
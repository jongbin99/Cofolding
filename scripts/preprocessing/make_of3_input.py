import json
import os
import pandas as pd
import argparse

def load_smiles_table(smiles_file: str) -> pd.DataFrame:
    if smiles_file.lower().endswith((".xlsx", ".xls")):
        df = pd.read_excel(smiles_file)
    else:
        df = pd.read_csv(smiles_file)
    # basic validation
    required = {"SMILES", "Dataset_ID"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required column(s): {sorted(missing)}. Found: {list(df.columns)}")
    return df

def process_smiles_openfold3(
    template_json: str,
    smiles_file: str,
    output_dir: str,
    protein_chain_id: str = "A",
    ligand_chain_id: str = "B",
    model_seeds=(1,),
    one_file_per_ligand: bool = True,
    output_filename: str = "of3_queries.json",
):
    """
    template_json must contain the protein sequence information.
    Two supported template styles:
      (1) Full OpenFold3-like template with a single query already present
      (2) Minimal template with top-level keys: {"protein_sequence": "..."} OR {"sequence": "..."}
    """

    os.makedirs(output_dir, exist_ok=True)

    with open(template_json, "r") as f:
        tmpl = json.load(f)

    # --- extract protein sequence from template ---
    protein_seq = None

    # Case 1: user gives a full OF3-ish JSON and we steal the first protein chain sequence
    if isinstance(tmpl, dict) and "queries" in tmpl and isinstance(tmpl["queries"], dict) and len(tmpl["queries"]) > 0:
        first_query = next(iter(tmpl["queries"].values()))
        if isinstance(first_query, dict) and "chains" in first_query:
            for ch in first_query["chains"]:
                if ch.get("molecule_type") == "protein" and "sequence" in ch:
                    protein_seq = ch["sequence"]
                    break

    # Case 2: minimal templates
    if protein_seq is None:
        protein_seq = tmpl.get("protein_sequence") or tmpl.get("sequence")

    if not protein_seq:
        raise ValueError(
            "Could not find protein sequence in template_json. "
            "Provide either an OpenFold3-like template with queries->...->chains including a protein sequence, "
            "or a minimal JSON with key 'protein_sequence' or 'sequence'."
        )

    # canonical protein chain object
    protein_chain = {
        "chain_ids": protein_chain_id,
        "molecule_type": "protein",
        "sequence": protein_seq,
    }

    df = load_smiles_table(smiles_file)

    if one_file_per_ligand:
        # write one OF3 json per ligand (each file has one query)
        for _, row in df.iterrows():
            smiles = str(row["SMILES"]).strip()
            dataset_id = str(row["Dataset_ID"]).strip()

            of3_json = {
                "queries": {
                    dataset_id: {
                        "chains": [
                            protein_chain,
                            {
                                "chain_ids": ligand_chain_id,
                                "molecule_type": "ligand",
                                "smiles": smiles,
                            },
                        ],
                        "model_seeds": list(model_seeds),
                    }
                }
            }

            outpath = os.path.join(output_dir, f"{dataset_id}.json")
            with open(outpath, "w") as out:
                json.dump(of3_json, out, indent=2)
    else:
        # write a single OF3 json containing all ligands as separate queries
        all_queries = {}
        for _, row in df.iterrows():
            smiles = str(row["SMILES"]).strip()
            dataset_id = str(row["Dataset_ID"]).strip()

            all_queries[dataset_id] = {
                "chains": [
                    protein_chain,
                    {
                        "chain_ids": ligand_chain_id,
                        "molecule_type": "ligand",
                        "smiles": smiles,
                    },
                ],
                "model_seeds": list(model_seeds),
            }

        of3_json = {"queries": all_queries}
        outpath = os.path.join(output_dir, output_filename)
        with open(outpath, "w") as out:
            json.dump(of3_json, out, indent=2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate OpenFold3 per-ligand query JSONs")
    parser.add_argument("--template_json", required=True,
                        help="Template JSON containing protein sequence (OpenFold3-like or minimal).")
    parser.add_argument("--smiles_file", required=True,
                        help="Excel/CSV with columns: SMILES, Dataset_ID")
    parser.add_argument("--output_dir", required=True,
                        help="Directory to write JSON(s)")

    parser.add_argument("--protein_chain_id", default="A")
    parser.add_argument("--ligand_chain_id", default="B")
    parser.add_argument("--model_seeds", default="1",
                        help="Comma-separated seeds, e.g. '1' or '1,2,3'")
    parser.add_argument("--single_file", action="store_true",
                        help="If set, write ONE JSON containing all ligands as separate queries.")
    parser.add_argument("--output_filename", default="of3_queries.json",
                        help="Used only with --single_file")

    args = parser.parse_args()
    seeds = tuple(int(x) for x in args.model_seeds.split(",") if x.strip())

    process_smiles_openfold3(
        template_json=args.template_json,
        smiles_file=args.smiles_file,
        output_dir=args.output_dir,
        protein_chain_id=args.protein_chain_id,
        ligand_chain_id=args.ligand_chain_id,
        model_seeds=seeds,
        one_file_per_ligand=(not args.single_file),
        output_filename=args.output_filename,
    )
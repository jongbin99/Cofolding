import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFMCS

def load_trained_smiles(trained_path, smiles_col="Ligand SMILES"):
    """
    Load trained-set SMILES. Supports:
      - .smi (whitespace-delimited; first column is SMILES)
      - .xlsx/.xls (Excel; uses smiles_col)
      - .csv/.tsv (uses smiles_col)
    """
    ext = os.path.splitext(trained_path)[1].lower()
    if ext == ".smi" or trained_path.lower().endswith(".smi.txt"):
        smiles = []
        with open(trained_path, "r") as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.strip().split()
                if len(parts) >= 1:
                    s = parts[0]
                    if s and s.lower() != "nan":
                        smiles.append(s)
        return smiles
    elif ext in [".xlsx", ".xls"]:
        df = pd.read_excel(trained_path)
    else:
        # allow csv/tsv fallback
        sep = "," if ext == ".csv" else "\t"
        df = pd.read_csv(trained_path, sep=sep)
    if smiles_col not in df.columns:
        raise ValueError(f"Column '{smiles_col}' not found in {trained_path}. Columns: {list(df.columns)}")
    # drop empty
    smiles = [s for s in df[smiles_col].astype(str).tolist() if s and s.lower() != "nan"]
    return smiles

def load_query_smiles(smiles_file):
    with open(smiles_file, "r") as f:
        lines = []
        for line in f:
            if line.strip():
                parts = line.strip().split()
                if len(parts) >= 1:
                    lines.append(parts[0])  # first column = SMILES
        return lines

def to_mol(s):
    try:
        m = Chem.MolFromSmiles(s)
        return m
    except Exception:
        return None

def valid_mols(smiles_list):
    out = []
    for s in smiles_list:
        m = to_mol(s)
        if m is not None:
            out.append((s, m))
    return out

def find_best_mcs_for_query(query_mol, trained_mols, timeout=10,
                            atom_compare=rdFMCS.AtomCompare.CompareElements,
                            bond_compare=rdFMCS.BondCompare.CompareOrder,
                            ring_matches_ring_only=True,
                            complete_rings_only=False,
                            match_valences=True):
    """
    Returns: (best_trained_smiles, mcs_smarts, mcs_num_atoms, mcs_num_bonds, coverage_query, coverage_trained)
    """
    best = (None, "", 0, 0, 0.0, 0.0)
    q_atoms = query_mol.GetNumAtoms()
    for tr_smiles, tr_mol in trained_mols:
        res = rdFMCS.FindMCS(
            [query_mol, tr_mol],
            timeout=timeout,
            ringMatchesRingOnly=ring_matches_ring_only,
            completeRingsOnly=complete_rings_only,
            matchValences=match_valences,
            atomCompare=atom_compare,
            bondCompare=bond_compare
        )
        if res.canceled:
            continue
        na, nb = res.numAtoms, res.numBonds
        # simple primary score: numAtoms; tie‑breaker: numBonds
        if (na > best[2]) or (na == best[2] and nb > best[3]):
            cov_q = na / max(1, q_atoms)
            cov_t = na / max(1, tr_mol.GetNumAtoms())
            best = (tr_smiles, res.smartsString, na, nb, cov_q, cov_t)
            # early exit if perfect coverage
            if na == q_atoms:
                break
    return best

def get_MCS_best_matches(smiles_file, trained_set, smiles_col="Ligand SMILES",
                         timeout=10, strict=True):
    # load inputs
    query_smiles = load_query_smiles(smiles_file)
    trained_smiles = load_trained_smiles(trained_set, smiles_col=smiles_col)

    # sanitize
    query = valid_mols(query_smiles)
    trained = valid_mols(trained_smiles)

    # choose comparison mode
    if strict:
        atom_cmp = rdFMCS.AtomCompare.CompareElements
        bond_cmp = rdFMCS.BondCompare.CompareOrder
        ring_only = True
        complete_rings = False
        match_val = True
    else:
        # looser, more scaffold‑like
        atom_cmp = rdFMCS.AtomCompare.CompareAny
        bond_cmp = rdFMCS.BondCompare.CompareAny
        ring_only = False
        complete_rings = False
        match_val = False

    rows = []
    for q_smiles, q_mol in query:
        best_tr, smarts, na, nb, cov_q, cov_t = find_best_mcs_for_query(
            q_mol, trained,
            timeout=timeout,
            atom_compare=atom_cmp,
            bond_compare=bond_cmp,
            ring_matches_ring_only=ring_only,
            complete_rings_only=complete_rings,
            match_valences=match_val
        )
        rows.append({
            "Query": q_smiles,
            "Best_Trained": best_tr,
            "MCS_SMARTS": smarts,
            "MCS_NumAtoms": na,
            "MCS_NumBonds": nb,
            "Coverage_Query": round(cov_q, 3),
            "Coverage_Trained": round(cov_t, 3)
        })

    return pd.DataFrame(rows)

if __name__ == "__main__":
    # Example usage:
    # python mcs_match.py queries.smi trained_set.smi
    import sys
    if len(sys.argv) < 3:
        sys.exit(1)
    df = get_MCS_best_matches(sys.argv[1], sys.argv[2], smiles_col="Ligand SMILES",
                              timeout=10, strict=True)
    out = "mcs_best_matches.csv"
    df.to_csv(out, index=False)
    print(f"Wrote {out}")

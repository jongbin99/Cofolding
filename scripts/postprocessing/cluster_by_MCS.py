import argparse
import os
import sys
import numpy as np
import pandas as pd

def load_mcs_matrix(xlsx_path, sheet_name="Sheet1", labels_source="auto"):
    """
    Load an NxN pairwise MCS% matrix from Excel.
    - If labels_source=='auto': use DataFrame index as labels; if columns match index, drop redundancy.
    - Ensures a float numpy array; converts percentages (>1) to 0..1.
    - Fills diagonal with 1.0 if missing/NaN.
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet_name, index_col=0)
    # If the first column is labels and columns are also labels, ensure square alignment
    # (Common pattern: index == columns)
    if not (df.shape[0] == df.shape[1]):
        raise ValueError(f"Matrix must be square; got {df.shape}")

    # Try to align columns to index if names are the same set but different order
    if set(df.index) == set(df.columns) and not df.index.equals(df.columns):
        df = df.loc[df.index, df.index]

    labels = df.index.astype(str).tolist()
    M = df.to_numpy(dtype=float)

    # Convert percent to fraction if needed
    if np.nanmax(M) > 1.0:
        M = M / 100.0

    # Diagonal to 1.0 (self-similarity)
    n = M.shape[0]
    if np.isnan(np.diag(M)).any():
        np.fill_diagonal(M, 1.0)
    else:
        np.fill_diagonal(M, 1.0)

    # Clip to [0,1] just in case
    M = np.clip(M, 0.0, 1.0)
    return M, labels


def best_first_clustering(sim_mat,
                          labels,
                          threshold=0.5,
                          max_clusters=None,
                          order="mean"):
    """
    Greedy Best-First Clustering on a similarity matrix (higher = more similar).
    - threshold: minimum similarity to a cluster center to join that cluster
    - max_clusters: cap on number of clusters (None = no cap)
    - order: 'mean' (default), 'max', or 'given' (original index order)

    Returns:
      assignments: np.array of length N with cluster IDs (0..K-1)
      centers: list of indices that are cluster centers
      details: list of tuples (cluster_id, mol_idx, mol_label, sim_to_center)
    """
    N = sim_mat.shape[0]
    assigned = np.full(N, -1, dtype=int)
    centers = []

    # Build order
    if order == "mean":
        # mean similarity (ignoring self on diagonal)
        S = sim_mat.copy()
        np.fill_diagonal(S, np.nan)
        prio = np.nanmean(S, axis=1)
        order_idx = np.argsort(-prio)  # descending
    elif order == "max":
        S = sim_mat.copy()
        np.fill_diagonal(S, np.nan)
        prio = np.nanmax(S, axis=1)
        order_idx = np.argsort(-prio)
    elif order == "given":
        order_idx = np.arange(N)
    else:
        raise ValueError("order must be one of {'mean','max','given'}")

    cluster_id = 0
    for i in order_idx:
        if assigned[i] != -1:
            continue
        if (max_clusters is not None) and (len(centers) >= max_clusters):
            break

        # Start a new cluster with i as center
        centers.append(i)
        assigned[i] = cluster_id

        # Assign unassigned members with similarity >= threshold to this cluster
        sims = sim_mat[i, :]
        members = np.where((assigned == -1) & (sims >= threshold))[0]
        assigned[members] = cluster_id

        cluster_id += 1

    # Any remaining unassigned (if max_clusters hit) become their own singletons
    # or remain -1; here we keep them -1, but you can promote them to a new cluster if desired.

    # Build details
    details = []
    for c_id, c_idx in enumerate(centers):
        # center row sims to members of that cluster
        in_cluster = np.where(assigned == c_id)[0]
        for m in in_cluster:
            sim_to_center = sim_mat[c_idx, m]
            details.append((c_id, int(m), labels[m], float(sim_to_center)))
        # ensure the center shows up even if alone
        if c_idx not in in_cluster:
            details.append((c_id, int(c_idx), labels[c_idx], 1.0))

    return assigned, centers, details


def write_outputs(out_prefix, assignments, centers, details, labels):
    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)

    # cluster_details.list: cluster_id, mol_index, mol_name, sim_to_center
    det_path = f"{out_prefix}_cluster_details.list"
    with open(det_path, "w", encoding="utf-8") as f:
        for c_id, idx, name, sim in details:
            f.write(f"{c_id},{idx},{name},{sim:.6f}\n")

    # counts
    valid = assignments[assignments >= 0]
    counts = np.bincount(valid) if valid.size else np.array([], dtype=int)
    cnt_path = f"{out_prefix}.count"
    with open(cnt_path, "w", encoding="utf-8") as f:
        for c in counts:
            f.write(f"{int(c)}\n")

    # centers
    cen_path = f"{out_prefix}_centers.txt"
    with open(cen_path, "w", encoding="utf-8") as f:
        for c in centers:
            f.write(f"{c}\t{labels[c]}\n")

    # assignments CSV (index,label,cluster)
    assign_path = f"{out_prefix}_assignments.csv"
    df = pd.DataFrame({
        "index": np.arange(len(labels)),
        "label": labels,
        "cluster": assignments
    })
    df.to_csv(assign_path, index=False)

    return {
        "details": det_path,
        "counts": cnt_path,
        "centers": cen_path,
        "assignments": assign_path
    }


def main():
    ap = argparse.ArgumentParser(
        description="Best-First Clustering from pairwise MCS% matrix in XLSX."
    )
    ap.add_argument("xlsx", help="Path to Excel with NxN MCS(%) matrix (index=rows, columns=labels).")
    ap.add_argument("-s", "--sheet", default=None, help="Sheet name (default: first).")
    ap.add_argument("-t", "--threshold", type=float, default=0.5,
                    help="Similarity cutoff (0..1). If your Excel is in %, we auto-scale.")
    ap.add_argument("-k", "--max-clusters", type=int, default=None,
                    help="Max clusters to form (default: no cap).")
    ap.add_argument("-o", "--out-prefix", default="bfc_mcs",
                    help="Output prefix for result files.")
    ap.add_argument("--order", choices=["mean", "max", "given"], default="mean",
                    help="Item processing order: mean (default), max, or given (original order).")
    args = ap.parse_args()

    sim, labels = load_mcs_matrix(args.xlsx, sheet_name=args.sheet)

    assignments, centers, details = best_first_clustering(
        sim_mat=sim,
        labels=labels,
        threshold=args.threshold,
        max_clusters=args.max_clusters,
        order=args.order
    )

    paths = write_outputs(args.out_prefix, assignments, centers, details, labels)

    # Brief summary
    valid = assignments[assignments >= 0]
    n_clusters = len(centers)
    total_assigned = int(valid.size)
    total = sim.shape[0]
    print(f"[INFO] Total items: {total}")
    print(f"[INFO] Assigned: {total_assigned} ({total_assigned/total:.1%})")
    print(f"[INFO] Clusters formed: {n_clusters}")
    if total_assigned:
        counts = np.bincount(valid)
        print(f"[INFO] Cluster size: mean={counts.mean():.2f}, median={np.median(counts):.0f}, "
              f"min={counts.min()}, max={counts.max()}")
    print("[INFO] Outputs:")
    for k, v in paths.items():
        print(f"  - {k}: {v}")


if __name__ == "__main__":
    main()
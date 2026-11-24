#!/usr/bin/env python3
"""
Cluster molecular structures using BitBIRCH-Lean (bblean) and identify
the most similar compound pairs assigned to different clusters.
"""

import argparse
import numpy as np
import pandas as pd
import bblean
import bblean.similarity as iSIM
from bblean_cluster import bblean_cluster


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Cluster molecules with bblean_cluster and find top cross-cluster similar pairs."
    )
    parser.add_argument(
        "--input",
        default="expansion_train.csv",
        help="Path to input CSV file (default: expansion_train.csv)",
    )
    parser.add_argument(
        "--smiles-col",
        default="SMILES",
        help="Name of the SMILES column (default: SMILES)",
    )
    parser.add_argument(
        "--name-col",
        default="Molecule Name",
        help="Name of the molecule identifier column (default: 'Molecule Name')",
    )
    parser.add_argument(
        "--clustered-output",
        default="expansion_train_clustered.csv",
        help="Path to save the clustered dataset CSV (default: expansion_train_clustered.csv)",
    )
    parser.add_argument(
        "--pairs-output",
        default="similar_compounds_different_clusters.csv",
        help="Path to save cross-cluster similar pairs CSV (default: similar_compounds_different_clusters.csv)",
    )
    parser.add_argument(
        "--min-sim",
        type=float,
        default=0.80,
        help="Minimum Tanimoto similarity threshold for saving cross-cluster pairs (default: 0.80)",
    )
    parser.add_argument(
        "--fp-kind",
        default="ecfp4",
        help="Fingerprint kind for bblean (default: ecfp4)",
    )
    parser.add_argument(
        "--n-features",
        type=int,
        default=2048,
        help="Number of fingerprint bits (default: 2048)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed BitBIRCH reclustering progress",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    # Load dataset
    df = pd.read_csv(args.input)
    print(f"Loaded {len(df)} molecules from {args.input}")

    # Run BitBIRCH-Lean clustering
    cluster_ids, _ = bblean_cluster(
        df,
        smiles_column=args.smiles_col,
        n_features=args.n_features,
        fp_kind=args.fp_kind,
        verbose=args.verbose,
    )
    df["cluster"] = cluster_ids

    n_clusters = df["cluster"].nunique()
    cluster_sizes = df["cluster"].value_counts().to_dict()
    n_singletons = sum(1 for size in cluster_sizes.values() if size == 1)

    print(f"Clustering complete: {n_clusters} clusters ({n_singletons} singletons)")
    df.to_csv(args.clustered_output, index=False)
    print(f"Saved clustered dataset to {args.clustered_output}")

    # Generate fingerprints and compute full pairwise Tanimoto similarity matrix
    fps = bblean.fps_from_smiles(
        df[args.smiles_col],
        pack=True,
        n_features=args.n_features,
        kind=args.fp_kind,
    )
    sim_matrix = iSIM.jt_sim_matrix_packed(fps)

    # Mask out lower triangle, diagonal, and pairs belonging to the same cluster
    n = len(df)
    c_arr = df["cluster"].values
    same_cluster_mask = c_arr[:, None] == c_arr[None, :]
    triu_mask = np.triu(np.ones((n, n), dtype=bool), k=1)
    valid_mask = triu_mask & (~same_cluster_mask)

    rows, cols = np.where(valid_mask & (sim_matrix >= args.min_sim))
    sims = sim_matrix[rows, cols]

    # Sort descending by similarity, breaking ties deterministically by row/col index
    order = np.lexsort((cols, rows, -sims))
    rows, cols, sims = rows[order], cols[order], sims[order]

    pairs_data = []
    for rank, (i, j, sim_val) in enumerate(zip(rows, cols, sims), 1):
        s1 = df.loc[i, args.smiles_col]
        s2 = df.loc[j, args.smiles_col]
        is_stereo = s1.replace("@", "") == s2.replace("@", "")
        pairs_data.append(
            {
                "Rank": rank,
                "Tanimoto_Similarity": round(float(sim_val), 6),
                "Molecule_1": df.loc[i, args.name_col],
                "SMILES_1": s1,
                "Cluster_1": int(c_arr[i]),
                "Cluster_1_Size": int(cluster_sizes[c_arr[i]]),
                "Molecule_2": df.loc[j, args.name_col],
                "SMILES_2": s2,
                "Cluster_2": int(c_arr[j]),
                "Cluster_2_Size": int(cluster_sizes[c_arr[j]]),
                "Is_Stereoisomer_Pair": is_stereo,
            }
        )

    pairs_df = pd.DataFrame(pairs_data)
    pairs_df.to_csv(args.pairs_output, index=False)
    print(
        f"Saved {len(pairs_df)} cross-cluster pairs (Tanimoto >= {args.min_sim:.2f}) to {args.pairs_output}"
    )

    if not pairs_df.empty:
        print("\nTop 10 Most Similar Pairs in Different Clusters:")
        cols_to_show = [
            "Rank",
            "Tanimoto_Similarity",
            "Molecule_1",
            "Cluster_1",
            "Molecule_2",
            "Cluster_2",
            "Is_Stereoisomer_Pair",
        ]
        print(pairs_df[cols_to_show].head(10).to_string(index=False))


if __name__ == "__main__":
    main()

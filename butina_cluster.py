import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
import bblean
import bblean.similarity as iSIM


def butina_cluster(
    dataframe: pd.DataFrame,
    smiles_column: str = "SMILES",
    n_features: int = 2048,
    radius: int = 2,
    dist_thresh: float | None = None,
    n_samples_for_std: int = 50,
    std_multiplier: int = 4,
    reordering: bool = True,
    verbose: bool = False,
) -> tuple[np.ndarray, list[tuple[int, ...]]]:
    """
    Clusters molecules in a DataFrame using RDKit's Taylor-Butina algorithm.

    If `dist_thresh` is None, the distance threshold is automatically set to
    `1.0 - (average_sim + std_multiplier * std)` using the same iSIM threshold
    estimation as `bblean_cluster.py`.

    Args:
        dataframe: DataFrame containing molecular data with SMILES strings.
        smiles_column: Name of the column containing SMILES strings. Defaults to "SMILES".
        n_features: Number of bits for Morgan fingerprint generation. Defaults to 2048.
        radius: Radius for Morgan fingerprint generation (2 = ECFP4). Defaults to 2.
        dist_thresh: Tanimoto distance cutoff (1 - similarity) for Butina clustering.
            If None, estimated automatically via iSIM.
        n_samples_for_std: Number of stratified samples for estimating similarity std dev
            when `dist_thresh` is None. Defaults to 50.
        std_multiplier: Multiplier for std dev when `dist_thresh` is None. Defaults to 4.
        reordering: Whether to recalculate neighbor counts for unassigned molecules after
            each cluster is formed (`reordering=True` in `Butina.ClusterData`). Defaults to True.
        verbose: If True, prints clustering threshold and summary statistics.

    Returns:
        A tuple containing:
        - numpy.ndarray: 1-indexed array of cluster IDs corresponding to rows in the input
          dataframe (ordered from largest cluster = 1 to smallest singletons).
        - list[tuple[int, ...]]: The raw Butina cluster tuples (sorted by cluster size descending).
    """
    if smiles_column not in dataframe.columns:
        raise ValueError(
            f"Column '{smiles_column}' not found in dataframe. Available columns: {list(dataframe.columns)}"
        )

    if len(dataframe) == 0:
        raise ValueError("DataFrame is empty. Cannot cluster empty dataset.")

    smiles_list = dataframe[smiles_column].tolist()
    n_mols = len(smiles_list)

    # Determine distance threshold via iSIM if not explicitly provided
    if dist_thresh is None:
        packed_fps = bblean.fps_from_smiles(
            smiles_list,
            pack=True,
            n_features=n_features,
            kind="ecfp4" if radius == 2 else "ecfp6",
        )
        average_sim = iSIM.jt_isim_packed(packed_fps)
        representative_samples = iSIM.jt_stratified_sampling(
            packed_fps, n_samples=n_samples_for_std
        )
        sim_matrix = iSIM.jt_sim_matrix_packed(packed_fps[representative_samples])
        sim_matrix = sim_matrix[~np.eye(sim_matrix.shape[0], dtype=bool)]
        std = float(np.std(sim_matrix))
        optimal_sim_threshold = float(average_sim + std_multiplier * std)
        dist_thresh = 1.0 - optimal_sim_threshold
        if verbose:
            print(
                f"iSIM estimated similarity threshold: {optimal_sim_threshold:.4f} "
                f"(Butina distance threshold: {dist_thresh:.4f})"
            )
    else:
        if verbose:
            print(
                f"Using explicit Butina distance threshold: {dist_thresh:.4f} "
                f"(similarity cutoff: {1.0 - dist_thresh:.4f})"
            )

    # Generate RDKit ExplicitBitVect Morgan fingerprints
    fpg = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_features)
    rd_fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            raise ValueError(f"Unable to parse SMILES: {smi}")
        rd_fps.append(fpg.GetFingerprint(mol))

    # Build lower-triangle pairwise Tanimoto distance list for Butina.ClusterData
    dists: list[float] = []
    for i in range(1, n_mols):
        sims = DataStructs.BulkTanimotoSimilarity(rd_fps[i], rd_fps[:i])
        dists.extend(1.0 - x for x in sims)

    # Run Taylor-Butina clustering
    clusters = Butina.ClusterData(
        dists, n_mols, dist_thresh, isDistData=True, reordering=reordering
    )
    clusters_sorted = sorted(clusters, key=len, reverse=True)

    # Assign 1-indexed cluster IDs sorted by cluster size descending
    assignments = np.zeros(n_mols, dtype=np.uint64)
    for cluster_id, member_indices in enumerate(clusters_sorted, start=1):
        for idx in member_indices:
            assignments[idx] = cluster_id

    if verbose:
        n_singletons = sum(1 for c in clusters_sorted if len(c) == 1)
        print(f"Final number of clusters: {len(clusters_sorted)}")
        print(f"Final number of singletons: {n_singletons}")

    return assignments, clusters_sorted

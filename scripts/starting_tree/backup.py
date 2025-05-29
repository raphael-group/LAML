from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, Optional, Union

import numpy as np
import pandas as pd
from numba import njit, prange
from skbio import DistanceMatrix, TreeNode

MissingType = Union[int, float]
MISSING_STATE: int = -1  # sentinel for missing characters

@njit  # type: ignore[nopython-warning]
def _whamming(a: np.ndarray, b: np.ndarray, w: np.ndarray, miss: int) -> float:
    """Numba-accelerated weighted-Hamming between *two* integer arrays."""
    numer = 0.0
    denom = 0.0
    for k in range(a.size):
        ai = a[k]
        bi = b[k]
        if ai == miss or bi == miss:
            continue
        denom += w[k]
        if ai != bi:
            numer += w[k]
    return np.nan if denom == 0.0 else numer / denom

def weighted_hamming_distance_matrix(
    char_matrix: pd.DataFrame,
    *,
    weights: Optional[Iterable[float]] = None,
    missing: MissingType = MISSING_STATE,
    parallel: bool = True,
) -> pd.DataFrame:
    """Return an *n × n* weighted-Hamming distance matrix (Numba-fast)."""

    data = char_matrix.to_numpy(dtype=np.int64, copy=True)
    n_samples, n_chars = data.shape

    w = np.ones(n_chars, dtype=np.float64) if weights is None else np.asarray(list(weights), dtype=np.float64)
    if w.size != n_chars:
        raise ValueError("`weights` must have one entry per column")

    dist = np.zeros((n_samples, n_samples), dtype=np.float64)
    outer = prange if parallel else range

    for i in outer(n_samples):
        for j in range(i + 1, n_samples):
            d = _whamming(data[i], data[j], w, int(missing))
            dist[i, j] = dist[j, i] = d

    return pd.DataFrame(dist, index=char_matrix.index, columns=char_matrix.index)

def _nj_skbio(dm: DistanceMatrix) -> TreeNode:
    from skbio.tree import nj  # local import keeps global import light

    return nj(dm)

def neighbor_joining_tree(
    char_matrix: pd.DataFrame,
    *,
    weights: Optional[Iterable[float]] = None,
    missing: MissingType = MISSING_STATE,
    nj_backend: str = "auto",  # {"auto", "ccphylo", "numba"}
) -> TreeNode:
    """Infer a **rooted** NJ tree and return it *without* the dummy `'root'` leaf."""

    cm = char_matrix.copy()
    cm.loc["root"] = 0

    dm_df = weighted_hamming_distance_matrix(cm, weights=weights, missing=missing)
    dm = DistanceMatrix(dm_df.values, ids=dm_df.index.tolist())
    tree = _nj_skbio(dm)

    tree = tree.root_at("root")
    tree = tree.shear([tip.name for tip in tree.tips() if tip.name != "root"])
    return tree

if __name__ == "__main__":
    # --- tiny illustrative character matrix --------------------------------
    cm = pd.DataFrame(
        [[1, 2, 3, 4], [1, 2, 4, 4], [3, 2, 4, -1], [3, -1, 4, 4]],
        index=["s1", "s2", "s3", "s4"],
        columns=["site1", "site2", "site3", "site4"],
        dtype=np.int64,
    )

    cm.loc["root"] = 0 

    # weight example: inverse of observed counts
    weights = 1 / cm.replace(MISSING_STATE, np.nan).count().values

    print("Distance matrix:\n", weighted_hamming_distance_matrix(cm, weights=weights), "\n")

    #tree = neighbor_joining_tree(cm, weights=weights)

    #newick_str = tree.to_newick() if hasattr(tree, "to_newick") else str(tree)
    #print("Rooted NJ tree (Newick):\n", newick_str)


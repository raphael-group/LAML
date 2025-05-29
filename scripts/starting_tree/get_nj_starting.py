from __future__ import annotations

import subprocess
import tempfile
from pathlib import Path
from typing import Iterable, Optional, Union, Dict
import sys

import numpy as np
import pandas as pd
from numba import njit, prange
from skbio import DistanceMatrix
from skbio.tree import nj, TreeNode

""" Read character matrix, produce rooted NJ tree topo (skbio) written to newick file. 
Rooting heuristic: Add an unedited sequence and use this as root.
Differs from cass in that it produces a binary tree, rather than collapsing mutationless edges. """

MissingType = Union[int, float]
MISSING_STATE: int = -1  # sentinel for missing characters

def nj_tree_from_distance_matrix(
    dm_df: pd.DataFrame,
    *,
    root_name: str = "root",
    prune_root: bool = False
) -> TreeNode:
    """
    Build a Neighbor-Joining tree from a square distance-matrix DataFrame and root it at `root_name`.
    Retruns a rooted skbio tree.
    """
    # ── force *string* IDs so scikit-bio never mistakes them for integers ──
    dm_df = dm_df.copy()
    dm_df.index = dm_df.index.map(str)
    dm_df.columns = dm_df.columns.map(str)
    root_name = str(root_name)
    dm = DistanceMatrix(dm_df.values, ids=dm_df.index.tolist())

    tree = nj(dm)                       # unrooted
    tree = tree.root_at(root_name)      # re-root on the dummy “root”

    # --- optionally remove the dummy root tip -----------------------------
    if prune_root:
        real_tips = [t.name for t in tree.tips() if t.name != root_name]
        tree = tree.shear(real_tips)

    return tree

def weighted_hamming_distance(s1, s2, miss=-1, weights: Optional[Dict[int, Dict[int, float]]] = None,) -> float:
    """
    Computes the weighted hamming distance between samples. 
    If weights are not given, then we increment dissimilarity by:
        - +2 if the states are different
        - +1 if one state is unedited and the other is an indel
        - +0 if the two states are identical 

    Logic same as in Cassiopeia: https://github.com/YosefLab/Cassiopeia/blob/master/cassiopeia/solver/dissimilarity_functions.py#L12
    """
    d, num_present = 0, 0
    for i in range(len(s1)):
        ai, bi = s1[i], s2[i]
        if ai == miss or bi == miss:
            continue
        num_present += 1

        if ai != bi:
            if ai == 0 or bi == 0:
                if weights:
                    if ai != 0:
                        d += weights[i][ai]
                    else:
                        d += weights[i][bi]
                else:
                    d += 1
            else:
                if weights:
                    d += weights[i][ai] + weights[i][bi]
                else:
                    d += 2

    if num_present == 0:
        return 0

    return d / num_present

def weighted_hamming_distance_matrix(
    cm: pd.DataFrame,
    *,
    miss: int = MISSING_STATE,
    weights: Optional[Union[Sequence[float], np.ndarray, Dict[int, Dict[int, float]]]] = None,
) -> pd.DataFrame:
    """Compute the full pair‑wise weighted‑Hamming distance matrix. Ignore weights."""

    taxa = cm.index.to_list()
    n = len(taxa)
    dmat = np.zeros((n, n), dtype=float)

    for i in range(n):
        si = cm.iloc[i].to_numpy()
        for j in range(i + 1, n):
            sj = cm.iloc[j].to_numpy()
            d = weighted_hamming_distance(si, sj, miss=miss, weights=None)
            dmat[i, j] = dmat[j, i] = d

    return pd.DataFrame(dmat, index=taxa, columns=taxa)

def _load_character_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, index_col=0)
    df = df.replace("?", MISSING_STATE).astype(np.int64)
    #total_missing = (df == MISSING_STATE).sum().sum()
    #print(f"{total_missing} entries contain the missing-state value ({MISSING_STATE}).")
    return df

def main():
    if len(sys.argv) != 3:
        sys.exit(
            "Usage: python get_nj_starting.py <character_matrix.csv> <out_tree.nwk>"
        )
    cm_path = Path(sys.argv[1])
    out_path = Path(sys.argv[2])
    if not cm_path.is_file():
        sys.exit(f"[ERROR] Character‑matrix file '{cm_path}' not found.")

    cm = _load_character_matrix(cm_path)
    cm.loc["root"] = 0
    
    dm_df = weighted_hamming_distance_matrix(cm)
    print(dm_df)
    tree = nj_tree_from_distance_matrix(dm_df, root_name="root", prune_root=True)

    # drop branch lengths and produce tree topology
    for node in tree.postorder():
        node.length = None

    # write to file
    tstr = str(tree)
    out_path.write_text(tstr + "\n")
    print(f"[✓] Rooted NJ topology with {tree.count(tips=True)} tips → {out_path}")

if __name__ == "__main__":
    main()

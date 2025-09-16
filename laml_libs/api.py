# laml_libs/api.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Literal, Dict, Any, Tuple, Union

from .runner import run_from_namespace  # shared core (CLI & API use this)

# ----------------------------- Types & Data ---------------------------------

Solver = Literal["Scipy", "EM", "fastEM-gpu", "fastEM-cpu"]
DelimiterKind = Literal["comma", "tab", "whitespace"]

PathLike = Union[str, Path]


@dataclass(frozen=True)
class LAMLParams:
    lambda_edit: float
    nu_silence: float
    phi_dropout: float
    neg_log_likelihood: float
    #output_status: Optional[str] = None


@dataclass(frozen=True)
class LAMLOutput:
    tree_newick: str
    params: LAMLParams
    annotations_path: Optional[Path]
    imputed_matrix_path: Optional[Path]
    log_path: Optional[Path]


# ------------------------------ Public API ----------------------------------
# note that the defaults are set in runner.py to match CLI defaults (run_laml.py)
def run_laml_infer(
    *,
    character_matrix: PathLike,
    tree_topology: Optional[PathLike] = None,
    priors: Optional[PathLike] = None,
    delimiter: Optional[DelimiterKind] = None,
    missing_data: Optional[str] = None,
    solver: Optional[Solver] = None,
    nInitials: Optional[int] = None,
    topology_search: bool = False,
    resolve_search: bool = False,
    keep_polytomies: bool = False,
    parallel: bool = False,
    randomreps: Optional[int] = None,
    maxIters: Optional[int] = None,
    timescale: Optional[float] = None, 
    noSilence: Optional[bool] = False,
    noDropout: Optional[bool] = False,
    compute_llh: Optional[bool] =  False,
    output_prefix: Optional[str] = None,
    verbose: Optional[bool] = False,
) -> LAMLOutput:
    """
    Programmatic entry point equivalent to the CLI (run_laml.py).
    Returns in-memory results while preserving the CLI's side effects (file outputs).

    Parameters mirror the CLI flags:
      - character_matrix (-c/--character_matrix): path to matrix file
      - tree_topology   (-t/--tree_topology): optional starting Newick
      - priors          (-p/--priors): optional priors file
      - delimiter       (--delimiter): 'comma' | 'tab' | 'whitespace'
      - missing_data    (--missing_data): token for missing values
      - solver          (--solver): 'Scipy' | 'EM' | 'fastEM-gpu' | 'fastEM-cpu'
      - nInitials       (--nInitials): EM inits
      - topology_search (--topology_search)
      - resolve_search  (--resolve_search)
      - keep_polytomies (--keep_polytomies)
      - parallel        (--parallel)
      - randomreps      (--randomreps)
      - maxIters        (--maxIters)
      - timescale       (--timescale)
      - noSilence       (--noSilence)
      - noDropout       (--noDropout)
      - compute_llh     (--compute_llh)
      - output_prefix   (-o/--output)
      - verbose         (--verbose)
    """

    def _s(p: Optional[PathLike]) -> Optional[str]:
        if p is None:
            return None
        return str(p)

    # Build a dict with canonical key names expected by runner.py.
    # (runner handles normalization & defaults too, but we keep names aligned.)
    args: Dict[str, Any] = {
        "character_matrix": _s(character_matrix),
        "tree_topology": _s(tree_topology),
        "priors": _s(priors),
        "delimiter": delimiter,
        "missing_data": missing_data,
        "solver": solver,
        "nInitials": int(nInitials),
        "topology_search": bool(topology_search),
        "resolve_search": bool(resolve_search),
        "keep_polytomies": bool(keep_polytomies),
        "parallel": bool(parallel),
        "randomreps": int(randomreps),
        "maxIters": maxIters if maxIters is None else int(maxIters),
        "timescale": float(timescale),
        "noSilence": bool(noSilence),
        "noDropout": bool(noDropout),
        "compute_llh": bool(compute_llh),
        "output": str(output_prefix),
        "verbose": bool(verbose),
    }

    tree_newick, params_dict, ann_path, imp_path, log_path = run_from_namespace(args)

    # Defensive conversions
    ann_p = Path(ann_path) if isinstance(ann_path, (str, Path)) else ann_path
    imp_p = Path(imp_path) if isinstance(imp_path, (str, Path)) else imp_path
    log_p = Path(log_path) if isinstance(log_path, (str, Path)) else log_path

    #print(params_dict.get("nll"))
    params = LAMLParams(
        lambda_edit=float(params_dict.get("lambda", float("nan"))),
        nu_silence=float(params_dict.get("nu", float("nan"))),
        phi_dropout=float(params_dict.get("phi", float("nan"))),
        neg_log_likelihood=float(params_dict.get("nllh", float("nan"))),
        #output_status=params_dict.get("nll")[1],
    )

    return LAMLOutput(
        tree_newick=str(tree_newick) if tree_newick is not None else "",
        params=params,
        annotations_path=ann_p,
        imputed_matrix_path=imp_p,
        log_path=log_p,
    )


__all__ = [
    "run_laml_infer",
    "LAMLOutput",
    "LAMLParams",
]


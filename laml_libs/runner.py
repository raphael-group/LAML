# laml_libs/runner.py
from __future__ import annotations
from pathlib import Path
from typing import Optional, Dict, Any, Tuple
import os
import re
import sys
import timeit
import datetime
from datetime import date
from copy import deepcopy

# Internal deps
import laml_libs as laml
from laml_libs.sequence_lib import read_sequences, read_priors, dedup, add_dup
from laml_libs.ML_solver import ML_solver
from laml_libs.EM_solver import EM_solver
from laml_libs.fastEM_solver import fastEM_solver, parse_data, parse_tree
from laml_libs.Topology_search_parallel import Topology_search_parallel as Topology_search_parallel
from laml_libs.Topology_search import Topology_search as Topology_search_sequential
from laml_libs.starting_tree import build_starting_tree
from treeswift import read_tree_newick
import jax

# -----------------------------------------------------------------------------
# Simple stdout logger (keeps your existing .log behavior)
# -----------------------------------------------------------------------------
class Logger(object):
    def __init__(self, output_prefix: str):
        self.terminal = sys.stdout
        self.log_path = f"{output_prefix}.log"
        self.log = open(self.log_path, "a")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    def flush(self):
        pass

# -----------------------------------------------------------------------------
# Shared core used by CLI and API
# -----------------------------------------------------------------------------
def run_from_namespace(args) -> Tuple[str, Dict[str, float], Optional[Path], Optional[Path], Optional[Path]]:
    """
    Returns:
      - tree_newick (str)
      - params_dict (keys: 'lambda', 'nu', 'phi', 'nll', 'mutation_rate')
      - annotations_path (Path|None)
      - imputed_matrix_path (Path|None)   [not produced here → None]
      - log_path (Path|None)
    """

    # 0) Normalize args to a dict (accept Namespace/_NS/dict)
    if not isinstance(args, dict):
        try:
            args = vars(args)
        except Exception:
            args = {k: getattr(args, k) for k in dir(args) if not k.startswith("_")}

    def as_bool(x, default=False):
        if x is None: return default
        if isinstance(x, bool): return x
        s = str(x).strip().lower()
        return s in {"1", "true", "t", "yes", "y"}
    # 1) Canonicalize expected keys (align with api.py)
    #    api.py passes: character_matrix, tree_topology, priors, delimiter, missing_data,
    #                   solver, nInitials, topology_search, resolve_search, keep_polytomies,
    #                   parallel, randomreps, maxIters, timescale, noSilence, noDropout,
    #                   compute_llh, output, verbose
    canon = {
        "character_matrix": args.get("character_matrix") or args.get("characters"),
        "tree_topology": args.get("tree_topology") or args.get("topology"),
        "priors": args.get("priors") or args.get("priorfile") or "uniform",
        "delimiter": (v if (v:=args.get("delimiter")) is not None else "comma"),
        "missing_data": (v if (v:=args.get("missing_data")) is not None else "?"),
        "solver": (v if (v:=args.get("solver")) is not None else "EM"),
        "nInitials": (int(v) if (v:=args.get('nInitials')) is not None else 20),
        "topology_search": as_bool(args.get("topology_search", False)),
        "resolve_search": as_bool(args.get("resolve_search", False)),
        "keep_polytomies": as_bool(args.get("keep_polytomies", False)),
        "parallel": as_bool(args.get("parallel", False)),
        "randomreps": (int(v) if (v:=args.get("randomreps")) is not None else 5),
        "maxIters": (int(v) if (v:=args.get("maxIters")) is not None else 500),
        "timescale": (float(v) if (v:=args.get("timescale")) is not None else 1.0),
        "noSilence": as_bool(args.get("noSilence") if args.get("noSilence") is not None else False),
        "noDropout": as_bool(args.get("noDropout", False)),
        "compute_llh": as_bool(args.get("compute_llh", False)),
        "output": args.get("output") or args.get("out") or "LAML_output",
        "verbose": as_bool(args.get("verbose", False)),
        # legacy/optional keys that may appear; give safe defaults
        "randseeds": args.get("randseeds", None),
        "ultra_constr": True, #bool(args.get("ultra_constr", False)),  # replaces noultrametric
    }

    prefix = str(canon["output"])
    sys.stdout = Logger(prefix)  # set up .log file early
    log_path = Path(f"{prefix}.log")

    # 2) MOSEK license checks (unchanged from your file)
    lic_file = os.path.join(os.path.expanduser("~"), 'mosek/mosek.lic')
    if 'MOSEKLM_LICENSE_FILE' not in os.environ and not os.path.isfile(lic_file):
        print("MOSEK license not found in environment variables. Please set the MOSEK license!")
        sys.exit(0)

    p0 = os.getenv("MOSEKLM_LICENSE_FILE", "").split(os.pathsep)[0]
    p = p0 if (p0 and os.path.isfile(p0)) else os.path.expanduser("~/mosek/mosek.lic")
    if not os.path.isfile(p): raise SystemExit("MOSEK license file not found")

    s = open(p, encoding="utf-8", errors="ignore").read()
    if "permanent" in s.lower():
        print("MOSEK license OK (permanent)")
        # Do NOT exit here; let the program continue
    else:
        ds = re.findall(r"\b\d{1,2}-[A-Za-z]{3}-\d{2,4}\b", s)
        if not ds: raise SystemExit("No expiry date found in license file")
        exp = min(datetime.datetime.strptime(d, "%d-%b-%Y" if len(d.split("-")[-1])==4 else "%d-%b-%y").date() for d in ds)
        if exp >= date.today():
            print(f"MOSEK license OK; expires {exp} ({(exp - date.today()).days} days left)")
        else:
            raise SystemExit(f"MOSEK license expired on {exp}")

    # 3) Inputs and basic validation
    cm_path = canon["character_matrix"]
    if not cm_path or not os.path.isfile(cm_path):
        print("Input files not found (character_matrix).")
        sys.exit(0)

    # Choose/prepare starting tree
    topo_path = canon["tree_topology"]
    if not topo_path:
        print("WARNING: No input topology found! Building NJ tree (heuristic root).")
        tree_file = f"{prefix}_input.nj.nwk"
        build_starting_tree(cm_path, tree_file)
    else:
        if not os.path.isfile(topo_path):
            print("WARNING: Provided topology not found! Building NJ tree (heuristic root).")
            tree_file = f"{prefix}_input.nj.nwk"
            build_starting_tree(cm_path, tree_file)
        else:
            tree_file = topo_path

    print("Launching " + laml.PROGRAM_NAME + " version " + laml.PROGRAM_VERSION)
    print(laml.PROGRAM_NAME + " was called as follows: " + " ".join(sys.argv))
    start_time = timeit.default_timer()

    # 4) Backend device selection based on solver
    solver_str = str(canon["solver"]).lower()
    if solver_str in {"fastem-cpu", "fastem-gpu"}:
        print(f"JAX version: {jax.__version__}")
        if solver_str == "fastem-gpu":
            os.environ['JAX_PLATFORM_NAME'] = 'gpu'
            # Respect user's JAX_PLATFORMS if set; otherwise expect 'cuda'
            if os.environ.get('JAX_PLATFORMS') != 'cuda':
                raise ValueError("JAX_PLATFORMS must be set to 'cuda' for fastEM-gpu.")
        else:
            os.environ['JAX_PLATFORM_NAME'] = 'cpu'
        print(f"Using JAX backend with {jax.devices()} devices.")
        print(f"Using device {jax.devices()[-1]} for computation.")
        jax.config.update("jax_default_device", jax.devices()[-1])
    else:
        os.environ['JAX_PLATFORM_NAME'] = 'cpu'

    if solver_str not in {"fastem-cpu", "fastem-gpu", "em", "scipy"}:
        print("Specified solver not recognized.")
        sys.exit(0)

    # 5) Preprocess: read inputs
    delim_map = {'tab': '\t', 'comma': ',', 'whitespace': ' '}
    delimiter = delim_map[canon["delimiter"]]
    msa, site_names = read_sequences(cm_path, filetype="charMtrx", delimiter=delimiter, masked_symbol=canon["missing_data"])

    # Load tree list
    with open(tree_file, 'r') as f:
        input_trees = [line.strip() for line in f]

    # Deduplicate sequences / trees
    has_dup, dup_map, msa, input_trees = dedup(msa, input_trees)

    k = len(msa[next(iter(msa.keys()))])
    # Fixed params / compute_llh mode
    if canon["compute_llh"]:
        # expecting "lambda phi nu" string; keep behavior as in your file
        fixed_lambda, fixed_phi, fixed_nu = [float(x) for x in str(canon["compute_llh"]).strip().split()]
    else:
        fixed_phi = 0 if canon["noDropout"] else None
        fixed_nu = 0 if canon["noSilence"] else None
        fixed_lambda = 1.0

    # Seeds (optional)
    if canon.get("randseeds") is None:
        random_seeds = None
    else:
        rs = [int(x) for x in str(canon["randseeds"]).strip().split()]
        random_seeds = rs[0] if canon["nInitials"] != 1 and len(rs) == 1 else rs

    # Priors
    if canon["priors"] == "uniform":
        print("No prior file detected, using uniform prior probabilities per site.")
        Q = []
        for i in range(k):
            M_i = set(msa[x][i] for x in msa if msa[x][i] not in [0, "?"])
            if len(M_i) == 0:
                m_i = 1
                q = {"1": 1.0}
            else:
                m_i = len(M_i)
                q = {x: 1/m_i for x in M_i}
            q[0] = 0
            Q.append(q)
    else:
        Q = read_priors(canon["priors"], msa, site_names=site_names)

    # Solver selection
    selected_solver = EM_solver
    em_selected = True
    fastem_selected = False
    if solver_str in {"fastem-cpu", "fastem-gpu"}:
        selected_solver = fastEM_solver
        fastem_selected = True
        if os.environ.get('XLA_PYTHON_CLIENT_PREALLOCATE') != 'false':
            print("Note: XLA_PYTHON_CLIENT_PREALLOCATE defaults to TRUE; recorded memory may look higher.")
    elif solver_str != "em":
        selected_solver = ML_solver
        em_selected = False

    # Pack data/prior/params
    data = {'charMtrx': msa}
    prior = {'Q': Q}
    params = {'nu': fixed_nu if fixed_nu is not None else laml.core_constants.eps,
              'phi': fixed_phi if fixed_phi is not None else laml.core_constants.eps}

    # fastEM: per-tree parsing and recoding
    if fastem_selected:
        data['charMtrx_recode'] = []
        data['ordered_leaf_labels'] = []
        prior['Q_recode'] = []
        for tstr in input_trees:
            swift_tree = read_tree_newick(tstr)
            parser_tree_out = parse_tree(swift_tree)
            out_phylogeny = parse_data(parser_tree_out, data, prior)
            ordered_leaf_labels = parser_tree_out['relabeling_vector'][:parser_tree_out['num_leaves']]
            data['charMtrx_recode'].append(out_phylogeny.character_matrix)
            data['ordered_leaf_labels'].append(ordered_leaf_labels)
            prior['Q_recode'].append(out_phylogeny.mutation_priors)

    # Topology search strategy
    Topology_search = Topology_search_sequential if not canon["parallel"] else Topology_search_parallel
    myTopoSearch = Topology_search(input_trees, selected_solver, data=data, prior=prior, params=params)

    # 6) Optimize / search
    if canon["compute_llh"]:
        print("Compute joint likelihood for provided trees and parameters (no optimization).")
        mySolver = myTopoSearch.get_solver()
        # rescale by lambda
        for tree in mySolver.trees:
            for node in tree.traverse_preorder():
                if node.edge_length is not None:
                    node.edge_length *= fixed_lambda
        nllh = mySolver.negative_llh()
        opt_trees = myTopoSearch.treeTopoList
        opt_params = myTopoSearch.params
        print("Tree negative log-likelihood:", nllh)
        print("Tree log-likelihood:", -nllh)
    else:
        my_strategy = deepcopy(laml.core_constants.DEFAULT_STRATEGY)
        # enforce ultrametric or not? (your original code referenced args["noultrametric"])
        my_strategy['ultra_constr'] = bool(canon.get("ultra_constr", False))
        resolve_polytomies = not canon["keep_polytomies"]
        my_strategy['resolve_search_only'] = canon["resolve_search"]
        if my_strategy['resolve_search_only'] and fastem_selected:
            print("Resolve-polytomies mode not implemented for fastEM_solver. Falling back to EM_solver.")
            selected_solver = EM_solver

        if not canon["resolve_search"] and not canon["topology_search"]:
            print("Optimizing branch lengths, phi, nu without topology search")
            if fastem_selected:
                print("Optimization by fast EM algorithm (JAX)")
            elif em_selected:
                print("Optimization by EM algorithm")
            else:
                print("Optimization by Scipy-SLSQP")
            mySolver = myTopoSearch.get_solver()
            nllh = mySolver.optimize(
                initials=canon["nInitials"],
                fixed_phi=fixed_phi, fixed_nu=fixed_nu,
                verbose=canon["verbose"],
                random_seeds=random_seeds,
                ultra_constr=True  # preserve prior behavior
            )
            myTopoSearch.update_from_solver(mySolver)
            opt_trees = myTopoSearch.treeTopoList
            opt_params = myTopoSearch.params
        else:
            if canon["resolve_search"]:
                if not resolve_polytomies:
                    print("WARNING: --resolve_search with --keep_polytomies → numeric optimization only.")
                else:
                    print("Starting local topology search to resolve polytomies")
                    if not myTopoSearch.has_polytomies:
                        print("No polytomy detected. Numeric optimization only.")
            else:
                if myTopoSearch.has_polytomy:
                    print("Detected polytomies in input trees.")
                    if not resolve_polytomies:
                        print("Flag --keep_polytomies on → keep all polytomies.")
                    else:
                        print("Flag --keep_polytomies off → resolve all polytomies.")
                else:
                    print("No polytomy detected.")
                print("Starting topology search")

            print("Running topology search {}...".format("in parallel" if canon["parallel"] else "sequentially"))
            checkpoint_file = f"{prefix}_ckpt.txt"
            opt_trees, max_score, opt_params = myTopoSearch.search(
                resolve_polytomies=resolve_polytomies,
                maxiter=canon["maxIters"],
                verbose=canon["verbose"],
                strategy=my_strategy,
                nreps=canon['randomreps'],
                checkpoint_file=checkpoint_file
            )
            nllh = -max_score

    # 7) Post-process and write outputs
    if has_dup:
        opt_trees = add_dup(opt_trees, dup_map)

    out_tree = f"{prefix}_trees.nwk"
    out_tree2 = f"{prefix}_trees.collapsed.nwk"
    out_annotate = f"{prefix}_annotations.txt"
    out_params = f"{prefix}_params.txt"

    # write scaled trees
    with open(out_tree, 'w') as fout:
        for tstr in opt_trees:
            tree = read_tree_newick(tstr)
            tree_height = tree.height(weighted=True)
            scaling_factor = tree_height / float(canon['timescale'])
            print(f"Tree height pre-scaling: {tree_height}, input timescale: {canon['timescale']}")
            for node in tree.traverse_preorder():
                node.edge_length = node.edge_length / scaling_factor
            tree_height = tree.height(weighted=True)
            mutation_rate = scaling_factor
            tparts = tree.__str__().split()
            if len(tparts) > 1:
                fout.write(''.join([tparts[0], "(", tparts[1][:-1], ");\n"]))
            else:
                fout.write(''.join(["(", tparts[0][:-1], ");\n"]))

    with open(out_tree2, 'w') as fout:
        for tstr in opt_trees:
            tree = read_tree_newick(tstr)
            tree.collapse_short_branches(threshold=laml.core_constants.dmin + (0.2 * laml.core_constants.dmin))
            tree_height = tree.height(weighted=True)
            scaling_factor = tree_height / float(canon['timescale'])
            for node in tree.traverse_preorder():
                if node.edge_length:
                    # problematic since we just collapsed short branches
                    node.edge_length = node.edge_length / scaling_factor
                else:
                    node.edge_length = laml.core_constants.dmin / scaling_factor
            tparts = tree.__str__().split()
            if len(tparts) > 1:
                fout.write(''.join([tparts[0], "(", tparts[1][:-1], ");\n"]))
            else:
                fout.write(''.join(["(", tparts[0][:-1], ");\n"]))

    # annotations (posterior formatting is unchanged)
    def format_posterior(p0, p_minus_1, p_alpha, alpha, q):
        if p0 == 1:
            return '0'
        elif p_minus_1 == 1:
            return '-1'
        elif p_alpha == 1:
            return alpha
        out = ''
        if p0 > 0:
            out += '0:' + str(p0)
        if p_minus_1 > 0:
            if out != '':
                out += '/'
            out += '-1:' + str(p_minus_1)
        if p_alpha > 0:
            if out != '':
                out += '/'
            if alpha == '?':
                out += "/".join([str(y) + ':' + str(round(p_alpha*q[y],3)) for y in q if round(p_alpha*q[y],3)>0])
            else:
                out += alpha + ":" + str(p_alpha)
        return out

    my_solver = EM_solver(opt_trees, {'charMtrx': msa}, {'Q': Q}, {'phi': opt_params['phi'], 'nu': opt_params['nu']})
    my_solver.az_partition()
    my_solver.Estep()
    idx = 0
    with open(out_annotate, 'w') as fout:
        for tree in my_solver.trees:
            # add root
            if len(tree.root.children) > 1:
                from treeswift import Node
                root = Node()
                root.label = 'I0'
                k2 = len(tree.root.alpha)
                root.alpha = ['z'] * k2
                root.post0 = [0] * k2
                root.post1 = [-float("inf")] * k2
                idx = 1
                root.add_child(tree.root)
                tree.root = root
            # branch length by expected #mutations
            all_labels = set()
            for node in tree.traverse_preorder():
                if node.is_root():
                    continue
                if node.label is None or node.label in all_labels:
                    node.label = 'I' + str(idx)
                    idx += 1
                all_labels.add(node.label)
                node.edge_length = round(sum(node.S1) + sum(node.S2) + sum(node.S4), 3)

            fout.write(tree.newick() + "\n")

            for node in tree.traverse_preorder():
                node.posterior = ''
                for j in range(len(node.alpha)):
                    if not node.is_root():
                        if node.alpha[j] == '?' and node.parent.alpha[j] != 'z':
                            node.alpha[j] = node.parent.alpha[j]
                    import math as _math
                    p0 = round(_math.exp(node.post0[j]), 2)
                    p_minus_1 = round(_math.exp(node.post1[j]), 2)
                    p_alpha = round(1 - p0 - p_minus_1, 2)
                    if node.posterior != '':
                        node.posterior += ','
                    node.posterior += format_posterior(p0, p_minus_1, p_alpha, str(node.alpha[j]), Q[j])
                fout.write(node.label + "," + str(node.posterior) + "\n")

    # params file
    with open(out_params, 'w') as fout:
        fout.write("Dropout rate: " + str(opt_params['phi']) + "\n")
        fout.write("Silencing rate: " + str(opt_params['nu']) + "\n")
        fout.write("Negative-llh: " + str(nllh) + "\n")
        fout.write("Mutation rate: " + str(mutation_rate) + "\n")

    stop_time = timeit.default_timer()
    print("Runtime (s):", stop_time - start_time)

    # Final return values
    solved_trees = [t.newick() for t in my_solver.trees]
    tree_newick_str = solved_trees[0] if solved_trees else ""
    params_dict = {'lambda': mutation_rate, 'phi': opt_params['phi'], 'nu': opt_params['nu'], 'nll': nllh}

    annotations_path = Path(out_annotate)
    imputed_matrix_path = None  # not produced here
    return tree_newick_str, params_dict, annotations_path, imputed_matrix_path, log_path


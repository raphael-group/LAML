"""
Author: Henri Schmidt
"""
import time
import timeit
import argparse
import math
import random
import sys
import jax
import fast_em.phylogeny
import os
import json

import optax
import optax.tree_utils as otu

import numpy as np
import pandas as pd
import networkx as nx
import jax.numpy as jnp
import loguru as lg
#import optimistix as optx
#import equinox.internal as eqxi

import fast_em.calculations as calc
from collections import defaultdict
from typing import Callable
from functools import partial

""" 
LAML: Lineage Analysis with Maximum Likelihood 

Code specific to LAML for computing the inside log likelihoods.
is included in the JIT-compiled functions:
    - `compute_internal_log_likelihoods`.
    - `initialize_leaf_inside_log_likelihoods`.

In the subsequent code, we assume that the alphabet
is ordered {0, 1, ..., A, -1}, where A is the size 
of the alphabet, 0 is the missing state, 1, ..., A
are non-missing states, and -1 is the unknown (?) state.
"""

M_STEP_DESCENT_STEPS  = 100
EM_STOPPING_CRITERION = 1e-5

def M_step_loss_fn(parameters, args):
    log_branch_lengths, logit_model_parameters = parameters
    edge_responsibilities, leaf_responsibilities, num_missing, num_not_missing = args[:4]
    branch_lengths = jnp.exp(log_branch_lengths)
    ν, ϕ           = jax.nn.sigmoid(logit_model_parameters)
    c1 = -edge_responsibilities[:, 0] * branch_lengths *  (1.0 + ν)
    c2 =  edge_responsibilities[:, 1] * (jnp.log(1 - jnp.exp(-branch_lengths)) - branch_lengths * ν)
    c3 =  edge_responsibilities[:, 2] * jnp.log(1 - jnp.exp(-branch_lengths * ν))
    c4 = -edge_responsibilities[:, 3] * branch_lengths * ν
    c5 =  edge_responsibilities[:, 4] * jnp.log(1 - jnp.exp(-branch_lengths * ν))
    c6 = num_not_missing * jnp.log(1 - ϕ) + (num_missing - jnp.sum(leaf_responsibilities)) * jnp.log(ϕ)
    return -(jnp.sum(c1 + c2 + c3 + c4 + c5) + c6)

jit_compute_log_likelihood = jax.jit(calc.compute_log_likelihood)
M_step_loss_fn_grad = jax.value_and_grad(M_step_loss_fn)
    
linesearch = optax.scale_by_backtracking_linesearch(max_backtracking_steps=15)
opt = optax.chain(optax.lbfgs(scale_init_precond=True), linesearch)

@jax.jit
def M_step_descent_step(params, state, args):
    branch_mask = args[-2]
    parameter_mask = args[-1]

    loss, grad = M_step_loss_fn_grad(params, args)
    grad = (grad[0] * branch_mask, grad[1] * parameter_mask)

    updates, state = opt.update(
        grad, state, params, value=loss, grad=grad, value_fn=M_step_loss_fn, args=args
    )

    params = optax.apply_updates(params, updates)
    return loss, params, state
    
@jax.jit
def M_step(params, args):
    def body_fun(carry, _):
        params, state = carry
        _, params, state = M_step_descent_step(params, state, args)
        return (params, state), None
    
    state = opt.init(params)
    (params, state), _ = jax.lax.scan(body_fun, (params, state), jnp.arange(M_STEP_DESCENT_STEPS))
    return params, M_STEP_DESCENT_STEPS

def optimize_parameters_expectation_maximization(
    leaves : jnp.array,
    internal_postorder : jnp.array,
    internal_postorder_children : jnp.array,
    parent_sibling : jnp.array,
    level_order : jnp.array,
    inside_log_likelihoods : jnp.array,
    model_parameters : jnp.array,
    model_parameters_mask : jnp.array,
    character_matrix : jnp.array,
    branch_lengths : jnp.array,
    branch_mask : jnp.array,
    mutation_priors : jnp.array,
    root : int,
    verbose = False
):
    model_parameters = jnp.maximum(model_parameters, calc.EPS)
    model_parameters = jnp.minimum(model_parameters, 1.0 - calc.EPS)

    logit_model_parameters = jnp.log(model_parameters / (1.0 - model_parameters))
    log_branch_lengths = jnp.log(jnp.maximum(branch_lengths, calc.EPS))

    # Assuming that missing states are represented by -1 in the character matrix
    num_missing = jnp.sum(character_matrix == -1, axis=[0,1]) # shape (num_leaves,)
    num_not_missing = jnp.sum(character_matrix != -1, axis=[0,1]) # shape (num_leaves,)

    # uses equation (12) in LAML manuscript to setup the loss function
    params = (log_branch_lengths, logit_model_parameters)

    compute_llh = lambda params: jit_compute_log_likelihood(
        jnp.exp(params[0]), mutation_priors, leaves, 
        internal_postorder, internal_postorder_children,
        parent_sibling, level_order, inside_log_likelihoods,
        jax.nn.sigmoid(params[1]), character_matrix, root
    )
    
    previous_nllh = compute_llh(params)

    current_nllh = jnp.inf
    iteration = 0

    while True:
        iteration += 1

        _, edge_responsibilities, leaf_responsibilities = calc.compute_E_step(
            jnp.exp(params[0]), mutation_priors, leaves, 
            internal_postorder, internal_postorder_children, 
            parent_sibling, level_order, inside_log_likelihoods, 
            jax.nn.sigmoid(params[1]), character_matrix, root
        )
        
        args = (
            edge_responsibilities, leaf_responsibilities, 
            num_missing, num_not_missing, 
            branch_mask, model_parameters_mask
        )

        params, its = M_step(params, args)
        current_nllh = compute_llh(params)

        if verbose:
            lg.logger.info(f"M-step descent steps: {its}")
            lg.logger.info(f"EM iteration {iteration}, Previous NLLH: {previous_nllh}, Current NLLH: {current_nllh}")
            lg.logger.info(f"Relative Improvement: {jnp.abs(current_nllh - previous_nllh) / jnp.abs(previous_nllh)}")

        if jnp.abs(current_nllh - previous_nllh) / jnp.abs(previous_nllh) < EM_STOPPING_CRITERION:
            break

        previous_nllh = current_nllh

    return current_nllh, jnp.exp(params[0]), jax.nn.sigmoid(params[1]), iteration

class EMOptimizer:
    def __init__(self, mutation_priors, character_matrix, verbose=False):
        self.mutation_priors  = mutation_priors   # (num_characters, max_alphabet_size) 
        self.character_matrix = character_matrix # (num_leaves, num_characters)
        self.num_leaves       = character_matrix.shape[0]
        self.alphabet_size    = mutation_priors.shape[1]
        self.num_characters   = character_matrix.shape[1]
        self.verbose          = verbose

        self.inside_log_likelihoods = jnp.zeros((self.num_characters, self.num_leaves * 2 - 1, self.alphabet_size + 2), dtype=jnp.float32)

    def optimize(
        self,
        tree : nx.DiGraph,
        branch_lengths : jnp.array, # (2 * num_leaves - 1,)
        nu : float,
        phi : float,
        branch_mask : jnp.array,    # (2 * num_leaves - 1,) binary mask
        fixed_nu : bool = False,
        fixed_phi : bool = False
    ):
        """
        Fits branch lengths and model parameters to a phylogenetic tree using EM. Uses the
        passed in `branch_lengths` and `model_parameters` as initial values for the optimization.

        Parameters:
            tree (nx.DiGraph): A directed graph representing the phylogenetic tree. The tree should 
                               be a binary tree with (2 * num_leaves - 1) nodes.
            branch_lengths (jnp.array): A 1D array representing branch lengths for all nodes in the tree.
                                        Its expected shape is (2 * num_leaves - 1,).
            nu (float): A parameter used in the optimization process.
            phi (float): A parameter used in the optimization process.
            branch_mask (jnp.array): A binary mask as a 1D array (with shape (2 * num_leaves - 1,)) 
                                     used to indicate which branches are considered for optimization.
        Returns:
            The optimized negative log likelihood, the optimized branch lengths, the optimized model parameters,
            and the number of EM iterations performed.
        """
        root = [n for n in tree.nodes() if tree.in_degree(n) == 0][0]
        leaves = jnp.array([n for n in tree.nodes() if tree.out_degree(n) == 0])
        level_order = nx.single_source_shortest_path_length(tree, root)
        internal_postorder = [[n, level_order[n]] for n in nx.dfs_postorder_nodes(tree, root) if tree.out_degree(n) > 0]
        internal_postorder = jnp.array(internal_postorder)
        internal_postorder_children = jnp.array([list(tree.successors(int(n))) for n in internal_postorder[:, 0]])
        level_order_jax = jnp.array([level_order[n] for n in range(2 * self.num_leaves - 1)])

        parent_sibling = []
        for i in range(2 * self.num_leaves - 1):
            if tree.in_degree(i) == 0:
                parent_sibling.append([-1, -1])
                continue
            parent = list(tree.predecessors(i))[0]
            siblings = list(tree.successors(parent))
            siblings.remove(i)
            parent_sibling.append([parent, siblings[0]])
        parent_sibling = jnp.array(parent_sibling)

        llh, branch_lengths, (nu, phi), em_iterations = optimize_parameters_expectation_maximization(
            leaves, 
            internal_postorder, 
            internal_postorder_children, 
            parent_sibling,
            level_order_jax,
            self.inside_log_likelihoods, 
            jnp.array([nu, phi]), 
            jnp.array([0.0 if fixed_nu else 1.0, 0.0 if fixed_phi else 1.0]),
            self.character_matrix, 
            branch_lengths, 
            branch_mask,
            self.mutation_priors, 
            root,
            verbose=self.verbose
        )

        return {
            "nllh": -llh,
            "branch_lengths": branch_lengths,
            "nu": nu,
            "phi": phi,
            "em_iterations": em_iterations
        }

def main(mode, phylo_opt):
    phylogeny = phylo_opt.phylogeny

    leaves = jnp.array([n for n in phylogeny.tree.nodes() if phylogeny.tree.out_degree(n) == 0])
    level_order = nx.single_source_shortest_path_length(phylogeny.tree, phylogeny.root)
    internal_postorder = [[n, level_order[n]] for n in nx.dfs_postorder_nodes(phylogeny.tree, phylogeny.root) if phylogeny.tree.out_degree(n) > 0]
    internal_postorder = jnp.array(internal_postorder)
    internal_postorder_children = jnp.array([list(phylogeny.tree.successors(int(n))) for n in internal_postorder[:, 0]])
    level_order_jax = jnp.array([level_order[n] for n in range(2 * phylogeny.num_leaves - 1)])
    
    parent_sibling = []
    for i in range(2 * phylogeny.num_leaves - 1):
        if phylogeny.tree.in_degree(i) == 0:
            parent_sibling.append([-1, -1])
            continue
        parent = list(phylogeny.tree.predecessors(i))[0]
        siblings = list(phylogeny.tree.successors(parent))
        siblings.remove(i)
        parent_sibling.append([parent, siblings[0]])
    parent_sibling = jnp.array(parent_sibling)

    if mode == "score":
        def llh_helper():
            return calc.compute_log_likelihood(
                phylo_opt.branch_lengths,
                phylogeny.mutation_priors, 
                leaves, 
                internal_postorder, 
                internal_postorder_children, 
                parent_sibling,
                level_order_jax,
                phylo_opt.inside_log_likelihoods, 
                phylo_opt.model_parameters,
                phylogeny.character_matrix, 
                phylogeny.root
            )

        llh_helper = jax.jit(llh_helper)
        llh_helper().block_until_ready()
        NUM_ITER = 200
        llh = llh_helper()
        runtime = timeit.timeit(lambda: llh_helper().block_until_ready(), number=NUM_ITER)
        avg_runtime = runtime / NUM_ITER

        root = [n for n in phylogeny.tree.nodes() if phylogeny.tree.in_degree(n) == 0][0]
        lg.logger.info(f"Log likelihood at root {root}: {llh}")
        lg.logger.info(f"Average runtime (s): {avg_runtime}")
    elif mode == "optimize-em":
        def optimize_helper(verbose):
            return optimize_parameters_expectation_maximization(
                leaves, 
                internal_postorder, 
                internal_postorder_children, 
                parent_sibling,
                level_order_jax,
                phylo_opt.inside_log_likelihoods, 
                phylo_opt.model_parameters, 
                jnp.array([1.0, 1.0]),
                phylogeny.character_matrix, 
                phylo_opt.branch_lengths, 
                jnp.ones(2 * phylogeny.num_leaves - 1),
                phylogeny.mutation_priors, 
                phylogeny.root,
                verbose=verbose
            )

        start = time.time()
        optimize_helper(False)[0] # warm start jit compiled functions
        end = time.time()
        compile_time = end - start

        start = time.time()
        nllh, post_branch_lengths, model_parameters, em_iterations = optimize_helper(True)
        nllh.block_until_ready()
        end = time.time()

        optimizer_specific_results = {"em_iterations": em_iterations}

        lg.logger.info(f"Compile time (s): {compile_time}, Optimization time (s): {end - start}")
        lg.logger.info(f"Optimized negative log likelihood(s): {nllh}")
        # lg.logger.info(f"Optimized branch lengths: {branch_lengths}")
        lg.logger.info(f"Optimized ν: {model_parameters[0]}, Optimized ϕ: {model_parameters[1]}")

    if "optimize" in mode:
        with open(f"{args.output}_results.json", "w") as f:
            res = {
                "nllh": -nllh.item(),
                "nu": model_parameters[0].item(),
                "phi": model_parameters[1].item(),
                "runtime": end - start,
                "compile_time": compile_time
            }

            res = res | optimizer_specific_results

            f.write(json.dumps(res))

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-c", "--character_matrix", help="Character matrix.", required=True)
    p.add_argument("-t", "--tree", help="Newick tree.", required=True)
    p.add_argument("-p", "--priors", help="Mutation priors CSV.")
    p.add_argument("-o", "--output", help="Prefix for output files.", default="output")
    p.add_argument("--nu", help="Heritable silencing rate (ν).", type=float, default=0.0)
    p.add_argument("--phi", help="Sequencing dropout rate (ϕ).", type=float, default=0.0)
    p.add_argument("--mode", help="Algorithm mode.", default="score", choices=["score", "optimize-em"])
    return p.parse_args()

if __name__ == "__main__":
    args = parse_args()

    lg.logger.remove()
    lg.logger.add(
        sys.stdout,
        format="<green>{time:YYYY-MM-DD HH:mm:ss.SSS}</green> | <level>{level: <8}</level> | {message}",
        colorize=True
    )

    tree, n = phylogeny.parse_newick(args.tree)

    character_matrix = pd.read_csv(args.character_matrix, sep=",", index_col=0)
    character_matrix.index = character_matrix.index.astype(str)
    character_matrix.replace("?", -1, inplace=True)
    character_matrix = character_matrix.astype(int)

    if args.priors is None:
        lg.logger.info("No mutation priors provided. Assuming uniform priors.")
        rows = []
        for i, c in enumerate(character_matrix.columns):
            states = set(character_matrix[c].unique()) - set([0, -1])
            num_states = len(states)
            for s in states:
                rows.append({"character": i, "state": s, "probability": 1.0 / num_states})
        priors = pd.DataFrame(rows)
        priors.character = priors.character.astype(int)
        priors.state = priors.state.astype(int)
        priors.set_index(["character", "state"], inplace=True)
    else:
        priors = pd.read_csv(args.priors, sep=",", header=None)
        priors.columns = ["character", "state", "probability"]
        priors.character = priors.character.astype(int)
        priors.state = priors.state.astype(int)
        priors.set_index(["character", "state"], inplace=True)

    if n != character_matrix.shape[0]:
        lg.logger.error("The tree and character matrix have different numbers of taxa.")
        sys.exit(1)

    phylo = phylogeny.build_phylogeny(tree, n, character_matrix, priors)

    if any(tree.nodes[i]["branch_length"] is None for i in range(2 * n - 1)) or "optimize" in args.mode:
        if args.mode == "optimize-em" or args.mode == "optimize-direct":
            lg.logger.info("Optimization mode. Initializing all branch lengths to 1.0.")
        else:
            lg.logger.error("Some branch lengths are missing. Initializing all branch lengths to 1.0.")

        if args.mode == "optimize-direct":
            branch_lengths = jnp.ones(2 * n - 1)
        else:
            branch_lengths = jnp.array([random.uniform(0.01, 0.5) for i in range(2 * n - 1)]) # TODO: give names to constants
    else:
        branch_lengths = jnp.array([tree.nodes[i]["branch_length"] for i in range(2 * n - 1)])

    lg.logger.info(f"Using JAX backend with {jax.devices()} devices.")
    lg.logger.info(f"Using device {jax.devices()[-1]} for computation.")
    jax.config.update("jax_default_device", jax.devices()[-1])

    lg.logger.info(f"Tree has {n} taxa and {2 * n - 1} nodes.")
    lg.logger.info(f"Character matrix has {character_matrix.shape[1]} characters and an alphabet size of {phylo.max_alphabet_size}.")
    model_parameters = jnp.array([args.nu, args.phi])
    phylo_opt = phylogeny.PhylogenyOptimization(
        phylogeny=phylo, 
        branch_lengths=branch_lengths, 
        model_parameters=model_parameters,
        inside_log_likelihoods=jnp.zeros((phylo.num_characters, phylo.num_leaves * 2 - 1, phylo.max_alphabet_size + 2), dtype=jnp.float32)
    )

    main(args.mode, phylo_opt)

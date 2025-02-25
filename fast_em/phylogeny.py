"""
Author: Henri Schmidt (w/ edits from Gillian)
"""
import treeswift
import networkx as nx
import pandas as pd
import sys
import jax.numpy as jnp
import loguru as lg
import random
from dataclasses import dataclass
from copy import deepcopy

@dataclass 
class Phylogeny:
    """
    Represents a phylogeny as a NetworkX DiGraph, with a root node,
    a mapping from leaf node names to row indices in the character matrix,
    the character matrix itself, and mutation priors.

    The Phylogeny is a binary tree, where each node has a name. Branch 
    lengths are optional.

    The character matrix is an integer matrix with the entries in each
    column corresponding to the character states of the taxa. The values
    are 0-indexed, with -1 indicating unknown data, 0 indiciating the
    missing data, and 1, 2, ... indicating the character states.
    """
    num_leaves : int
    num_characters : int
    max_alphabet_size : int
    root : int
    tree : nx.DiGraph
    character_matrix : jnp.array # shape (num_leaves, num_characters)
    mutation_priors : jnp.array  # shape (num_characters, max_alphabet_size)

@dataclass
class PhylogenyOptimization:
    """
    Represents a phylogeny optimization problem, with a Phylogeny object,
    branch lengths, and model parameters.
    """
    phylogeny : Phylogeny 
    branch_lengths : jnp.array # shape (2 * num_leaves - 1,)
    model_parameters : jnp.array
    inside_log_likelihoods : jnp.array # shape (num_characters, 2 * num_leaves - 1, max_alphabet_size + 2)

def build_phylogeny(tree: nx.DiGraph, num_leaves: int, character_matrix: pd.DataFrame, mutation_priors: pd.DataFrame) -> Phylogeny:
    """
    Constructs a Phylogeny object from the parsed components.
    """
    roots = [n for n in tree.nodes() if tree.in_degree(n) == 0]
    if len(roots) != 1:
        lg.logger.error("The tree has more than one root.")
        sys.exit(1)
    root = roots[0]

    character_matrix_recode = character_matrix.copy()
    original_to_new_mappings = []
    max_alphabet_size = 0

    for col in character_matrix.columns:
        column_data = character_matrix[col]
        states = column_data.unique()
        valid_states = sorted([s for s in states if s not in {-1, 0}])
        k = len(valid_states)
        if k > max_alphabet_size:
            max_alphabet_size = k
        mapping = {orig: new+1 for new, orig in enumerate(valid_states)}
        original_to_new_mappings.append(mapping)
        character_matrix_recode[col] = column_data.replace(mapping)

    character_matrix_jax = jnp.array(character_matrix_recode.values, dtype=jnp.int32)
    num_characters = character_matrix_jax.shape[1]

    state_mappings = []
    for mapping in original_to_new_mappings:
        state_mappings.append({v: k for k, v in mapping.items()})

    mutation_priors_jax = jnp.zeros((num_characters, max_alphabet_size), dtype=jnp.float32)
    for c in range(num_characters):
        char_name = character_matrix.columns[c]
        valid_states = sorted(original_to_new_mappings[c].keys())
        k = len(valid_states)
        for new_idx in range(k):
            orig_state = valid_states[new_idx]
            try:
                if orig_state == -1 or orig_state == 0:
                    continue 
                prior = mutation_priors.loc[(c, orig_state), 'probability']
            except KeyError:
                lg.logger.error(f"Missing mutation prior for character {char_name}, state {orig_state}")
                sys.exit(1)
            mutation_priors_jax = mutation_priors_jax.at[c, new_idx].set(prior)

    leaf_to_row = jnp.zeros(num_leaves, dtype=jnp.int32)
    for node in tree.nodes():
        if tree.out_degree(node) == 0:
            leaf_name = tree.nodes[node]["label"]
            row_idx = character_matrix_recode.index.get_loc(leaf_name)
            leaf_to_row = leaf_to_row.at[node].set(row_idx)

    # reorder the character matrix to match the leaf order in the tree
    character_matrix_jax = character_matrix_jax[leaf_to_row, :]

    return Phylogeny(
        num_leaves=num_leaves,
        num_characters=num_characters,
        max_alphabet_size=max_alphabet_size,
        root=root,
        tree=tree,
        character_matrix=character_matrix_jax,
        mutation_priors=mutation_priors_jax
    )

def parse_newick(newick):
    """
    Parses a Newick tree into a NetworkX DiGraph, labeling each 
    node with an integer index and storing the node name and 
    branch length into each node as attributes. 

    Requirement: Each leaf node must have a label.

    Guarantee: The leaf nodes will be labeled as {0, 1, 2, ..., n - 1},
    where n is the number of taxa in the tree.
    """
    swift_tree = treeswift.read_tree_newick(newick)
    tree = nx.DiGraph()
    for idx, v in enumerate(swift_tree.traverse_preorder()):
        label = v.get_label()
        if label is None and v.is_leaf():
            lg.logger.error("Leaf node has no label.")
            sys.exit(1)
        elif label is None:
            label = f"node_{str(idx)}"

        tree.add_node(idx, label=label, branch_length=v.get_edge_length())
        v.set_label(idx)

        if v.is_root(): continue
        tree.add_edge(v.get_parent().get_label(), idx) 

    # input tree has a synthetic root with out-degree 1
    root = [n for n in tree.nodes() if tree.in_degree(n) == 0][0]
    if tree.out_degree(root) == 1:
        tree.remove_node(root)

    n, relabeling_map = 0, {}
    for node in tree.nodes():
        if tree.out_degree(node) == 0:
            relabeling_map[node] = n
            n += 1

    num_leaves = n
    for node in tree.nodes():
        if tree.out_degree(node) > 0:
            relabeling_map[node] = n
            n += 1

    tree = nx.relabel_nodes(tree, relabeling_map)
    return tree, num_leaves

def parse_swift_tree(swift_tree, has_branch_mask=False):
    """
    Parses Treeswift tree object into a NetworkX DiGraph, labeling each 
    node with an integer index and storing the node name and 
    branch length into each node as attributes. 
    
    Returns: 
    - NetworkX DiGraph tree w index numbers as labels
    - treeswift tree w/ new node names 
    - branch_lengths: floats vector in order of the index numbers. Initializes to [0.01,0.5] if not given.
    - branch_mask: boolean vector in order of the index numbers
    - relabeling_map (swift_tree node.label : networkx idx)

    Requirement: Each leaf node must have a label.

    Guarantee: The leaf nodes will be labeled as {0, 1, 2, ..., n - 1},
    where n is the number of taxa in the tree.
    """
    tree = nx.DiGraph()
    swift_tree = deepcopy(swift_tree)
    for idx, v in enumerate(swift_tree.traverse_preorder()):
        label = v.get_label()
        if label is None and v.is_leaf():
            lg.logger.error("Leaf node has no label.")
            sys.exit(1)
        elif label is None:
            label = f"node_{str(idx)}"

        if has_branch_mask:
            tree.add_node(label, label=label, branch_length=v.get_edge_length(), branch_mask=v.mark_fixed)
        else:
            tree.add_node(label, label=label, branch_length=v.get_edge_length())
        v.set_label(label)

        if v.is_root(): continue
        tree.add_edge(v.get_parent().get_label(), label) 

    # input tree has a synthetic root with out-degree 1
    root = [n for n in tree.nodes() if tree.in_degree(n) == 0][0]
    if tree.out_degree(root) == 1:
        tree.remove_node(root)

    branch_lengths = []
    branch_mask = []
    n, relabeling_map = 0, {}
    for node in tree.nodes():
        if tree.out_degree(node) == 0: # leaf nodes
            nx_label = tree.nodes[node]['label']
            relabeling_map[nx_label] = n
            n += 1
            bl = tree.nodes[node]['branch_length'] 
            bl = bl if bl is not None else random.uniform(0.01, 0.5)
            branch_lengths.append(bl)
            if has_branch_mask:
                branch_mask.append(not node['branch_mask'])
            else:
                branch_mask.append(True)

    # how should i initialize the branch lengths?

    num_leaves = n
    for node in tree.nodes():
        if tree.out_degree(node) > 0:
            nx_label = tree.nodes[node]['label']
            relabeling_map[nx_label] = n
            n += 1
            bl = tree.nodes[node]['branch_length'] 
            bl = bl if bl is not None else random.uniform(0.01, 0.5)
            branch_lengths.append(bl)
            if has_branch_mask:
                branch_mask.append(not node['branch_mask'])
            else:
                branch_mask.append(True)

    tree = nx.relabel_nodes(tree, relabeling_map)

    return tree, num_leaves, swift_tree, jnp.array(branch_lengths), jnp.array(branch_mask), relabeling_map

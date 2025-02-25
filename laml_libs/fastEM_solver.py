from fast_em.laml import *
from fast_em.phylogeny import build_phylogeny
from laml_libs.EM_solver import *
from math import exp,log
from treeswift import * 
import cvxpy as cp
from laml_libs import min_llh, conv_eps, eps
import numpy as np
import time
import networkx as nx
import pandas as pd
import os
import random
os.environ['JAX_PLATFORM_NAME'] = 'cpu'

def parse_data(parser_tree_out, data, prior):
    """
    Take in result object from parse_tree, the LAML data with data['charMtrx'] entry, and prior with prior['Q'] entry.
    """
    #nx_tree, num_leaves, swift_tree, branch_lengths, branch_mask, relabeling_map = parse_swift_tree(swift_tree)
    #leaf_labels = set([n.label for n in swift_tree.traverse_leaves()])
    ordered_leaf_labels = parser_tree_out['relabeling_vector'][:parser_tree_out['num_leaves']]
    #sorted(leaf_labels, key=lambda x: relabeling_map[x])

    # pre-process character matrix and mutation priors
    cmat_tmp = pd.DataFrame(data['charMtrx']).T.replace('?', -1).reindex(ordered_leaf_labels) # missing data is ? 
    # prior is a dictionary of dictionaries
    mutation_priors = prior['Q']
    data = []
    for column_index, state_dict in enumerate(mutation_priors): 
        for orig_state, probability in state_dict.items():
            data.append((column_index, orig_state, probability))

    mutation_priors = pd.DataFrame(data, columns=['c', 'orig_state', 'probability'])
    mutation_priors.set_index(['c', 'orig_state'], inplace=True) 
    mutation_priors['probability'] = pd.to_numeric(mutation_priors['probability'])
    prior_tmp = mutation_priors

    out_phylogeny = build_phylogeny(parser_tree_out['nxtree'], parser_tree_out['num_leaves'], cmat_tmp, prior_tmp)
    return out_phylogeny

def parse_tree(swift_tree, has_branch_mask=False):
    """ 
    Take in a treeswift tree and output:
    - networkx tree object
    - branch_lengths vector
    - branch_mask vector (optionally): True means to reoptimize branch lengths. By default, applies to all.
    - relabeling_vector: corresponds to swift_tree labels in the index order of branch_lengths

    Note: The only change to the original treeswift object is to 
    add internal node labels if they did not previously exist.
    Assumes all leaves are labeled.

    If branch lengths are none on the input swift_tree, we will initialize to random.uniform(0.01, 0.5).
    """
    nxtree = nx.DiGraph()
    for idx, v in enumerate(swift_tree.traverse_preorder()):
        label = v.get_label()
        if label is None:
            label = f"node_{str(idx)}"
        if has_branch_mask:
            nxtree.add_node(label, label=label, branch_length=v.get_edge_length(), branch_mask=v.mark_fixed)
        else:
            nxtree.add_node(label, label=label, branch_length=v.get_edge_length())
        v.set_label(label)

        if v.is_root(): continue
        nxtree.add_edge(v.get_parent().get_label(), label)
    # input tree has a synthetic root with out-degree 1
    root = [n for n in nxtree.nodes() if nxtree.in_degree(n) == 0][0]
    if nxtree.out_degree(root) == 1:
        nxtree.remove_node(root)

    # build branch_lengths, branch_mask, relabeling_vector: leaves first
    branch_mask_vector = []
    branch_lengths = []
    n, relabeling_map = 0, {}
    for node in nxtree.nodes():
        if nxtree.out_degree(node) == 0: # leaf nodes
            nx_label = nxtree.nodes[node]['label']
            relabeling_map[nx_label] = n #(nx_label)
            n += 1
            bl = nxtree.nodes[node]['branch_length']
            bl = bl if bl is not None else random.uniform(0.01, 0.5)
            branch_lengths.append(bl)
            if has_branch_mask:
                branch_mask_vector.append(not node['branch_mask'])
            else:
                branch_mask_vector.append(True)

    num_leaves = n
    # build branch_lengths, branch_mask, relabeling_vector: internal nodes next
    for node in nxtree.nodes():
        if nxtree.out_degree(node) > 0:
            nx_label = nxtree.nodes[node]['label']
            relabeling_map[nx_label] = n
            n += 1
            bl = nxtree.nodes[node]['branch_length']
            bl = bl if bl is not None else random.uniform(0.01, 0.5)
            branch_lengths.append(bl)
            if has_branch_mask:
                branch_mask_vector.append(not node['branch_mask'])
            else:
                branch_mask_vector.append(True)
    nxtree = nx.relabel_nodes(nxtree, relabeling_map)
    relabeling_vector = sorted([x for x in relabeling_map], key=lambda x: relabeling_map[x])

    parsed = {}
    parsed['nxtree'] = nxtree
    parsed['num_leaves'] = num_leaves
    #parsed['swift_tree'] = swift_tree
    parsed['branch_lengths'] = jnp.array(branch_lengths)
    parsed['branch_mask'] = jnp.array(branch_mask_vector)
    parsed['relabeling_vector'] = relabeling_vector

    #print("after parse_tree:", swift_tree.newick())

    return parsed


class fastEM_solver(EM_solver):

    def __init__(self,treeList,data,prior,params={'nu':0,'phi':0,'sigma':0}):
        super().__init__(treeList,data,prior,params={'nu':0,'phi':0,'sigma':0}) #, **kwargs)
        self.character_matrix_recode = None
        self.mutation_priors_recode = None
        self.myEMOptimizers = []
        self.nx_trees = []
        self.swift_trees = []

        for tidx, swift_tree in enumerate(self.trees):  # treeList
            # order it according to the relabeling map of nodes and drop the header and cell names
            parser_tree_out = parse_tree(swift_tree)
            out_phylogeny = parse_data(parser_tree_out, data, prior)

            # set up an object of class EMOptimizer
            self.myEMOptimizers.append(EMOptimizer(out_phylogeny.mutation_priors, out_phylogeny.character_matrix))
  
    def score_tree(self, strategy={'ultra_constr':False,'fixed_phi':None,'fixed_nu':None,'fixed_brlen':None}, compare=False): 
        if compare: 
            trees = deepcopy(self.trees)
            ultra_constr = strategy['ultra_constr']
            fixed_phi = strategy['fixed_phi']
            fixed_nu = strategy['fixed_nu']
            fixed_brlen = strategy['fixed_brlen']

            # take the best of the results

            #print(f"input to oldEM (should have internal node labels now): {[t.newick() for t in self.trees]}")
            start_time = time.time()
            oldEM_score , _ = self.score_tree_oldEM(strategy)
            oldEM_phi = self.params.phi
            oldEM_nu = self.params.nu
            end_time = time.time()
            oldEM_elapsed_time = end_time - start_time
            #if fixed_brlen: 
            oldEM_trees = deepcopy(self.trees)

            #print(f"after old EM: {[t.newick() for t in self.trees]}")
            ultra_constr = strategy['ultra_constr']
            fixed_phi = strategy['fixed_phi']
            fixed_nu = strategy['fixed_nu']
            fixed_brlen = strategy['fixed_brlen']
            
            #print(f"input to fastEM: {[t.newick() for t in trees]}")
            all_relabeling_vectors = []
            all_branch_lengths = []
            results = []
            # pass in tree from topology search
            start_time = time.time()
            for tidx, swift_tree in enumerate(trees): #self.trees):
                # assumes swift_tree already has branches marked?
                #nx_tree, num_leaves, _, branch_lengths, branch_mask, relabeling_map = parse_swift_tree(swift_tree)
                parser_tree_out = parse_tree(swift_tree)
                nu = fixed_nu if fixed_nu else 0.5 
                phi = fixed_phi if fixed_phi else 0.5
                
                # call laml API from fast_em
                out = self.myEMOptimizers[tidx].optimize(parser_tree_out['nxtree'], parser_tree_out['branch_lengths'], 
                                                         nu, phi, parser_tree_out['branch_mask'], fixed_nu, fixed_phi)

                #print("relabeling_vector", parser_tree_out['relabeling_vector'])
                all_relabeling_vectors.append(parser_tree_out['relabeling_vector'])

                nllh = out['nllh']
                branch_lengths = out['branch_lengths']
                nu = out['nu']
                phi = out['phi']
                em_iterations = out['em_iterations']
                self.params.phi = phi
                self.params.nu = nu

                all_branch_lengths.append([float(x) for x in branch_lengths])

            end_time = time.time()
            fastEM_elapsed_time = end_time - start_time

            score = -nllh
            results.append(score)

            #print(f"after fastEM: {[t.newick() for t in trees]}, {[list(x) for x in all_branch_lengths]}")
            print(f"Current fastEM_solver score: {score}, {self.params.phi}, {self.params.nu}, {fastEM_elapsed_time:.6f} seconds")

            oldEM_branch_lengths = []
            for tidx, t in enumerate(oldEM_trees):
                l2n = t.label_to_node(selection='all')
                oldEM_branch_lengths.append([l2n[n].get_edge_length() for n in all_relabeling_vectors[tidx]])

            print(f"avg brlen diff between EMsolvers (old - new): {[np.mean(np.array(oldEM_branch_lengths[i]) - np.array(all_branch_lengths[i])) for i in range(len(self.trees))]}") #[t.newick() for t in self.trees]}")
            print(f"llh diff between EMsolvers (old - new): {oldEM_score - max(results)}")
            print(f"time diff between EMsolvers (old - new): {oldEM_elapsed_time - fastEM_elapsed_time}")
            print(f"param (phi) diff between EMsolvers (old - new): {oldEM_phi - self.params.phi}")
            print(f"param (nu) diff between EMsolvers (old - new): {oldEM_nu - self.params.nu}")

            status = "optimal"
            # should return best nllh and status over tree list
            #score = None if nllh is None else -nllh
            return max(results), status
        else:
            ultra_constr = strategy['ultra_constr']
            fixed_phi = strategy['fixed_phi']
            fixed_nu = strategy['fixed_nu']
            fixed_brlen = strategy['fixed_brlen']

            if 'final_optimization' in strategy.keys():
                oldEM_score, oldEM_status = self.score_tree_oldEM(strategy)
                return oldEM_score, oldEM_status
            
            #print(f"input to fastEM: {[t.newick() for t in trees]}")
            all_relabeling_vectors = []
            all_branch_lengths = []
            results = []
            # pass in tree from topology search
            start_time = time.time()
            for tidx, swift_tree in enumerate(self.trees): #self.trees):
                # assumes swift_tree already has branches marked?
                parser_tree_out = parse_tree(swift_tree)
                nu = fixed_nu if fixed_nu else 0.5 
                phi = fixed_phi if fixed_phi else 0.5
                
                # call laml API from fast_em
                out = self.myEMOptimizers[tidx].optimize(parser_tree_out['nxtree'], parser_tree_out['branch_lengths'], 
                                                         nu, phi, parser_tree_out['branch_mask'], fixed_nu, fixed_phi)
                all_relabeling_vectors.append(parser_tree_out['relabeling_vector'])

                nllh = out['nllh']
                branch_lengths = out['branch_lengths']
                nu = out['nu']
                phi = out['phi']
                em_iterations = out['em_iterations']
                self.params.phi = phi
                self.params.nu = nu

                all_branch_lengths.append([float(x) for x in branch_lengths])
            end_time = time.time()
            fastEM_elapsed_time = end_time - start_time

            score = -nllh
            results.append(score)
            #print(f"Current fastEM_solver score: {score}, {self.params.phi}, {self.params.nu}, {fastEM_elapsed_time:.6f} seconds")
            status = "optimal"
            return max(results), status


    def score_tree_oldEM(self, strategy): # if final_optimization 
        #start_time = time.time()
        score,status = super().score_tree(strategy)
        oldEM_params = super().get_params()
        #end_time = time.time()
        #elapsed_time = end_time - start_time
        #print(f"Current oldEM_solver score: {score}, {oldEM_params['phi']}, {oldEM_params['nu']}, {elapsed_time:.6f} seconds")
        return score,status
        

import unittest
from laml_libs.sequence_lib import read_sequences, read_priors
from laml_libs.EM_solver import EM_solver 
from laml_libs.fastEM_solver import fastEM_solver, parse_data, parse_tree
from laml_libs.ML_solver import ML_solver
from treeswift import *
from math import log
from random import random
from laml_libs import DEFAULT_STRATEGY
from copy import deepcopy
import pkg_resources

class fastEMTest(unittest.TestCase):
    # test optimize_one
    def test_1(self):
        treedata_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test4.tre')
        msa_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test4_charMtrx.txt')
        prior_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test4_prior.csv')
        T = read_tree_newick(treedata_path) #"test_data/test_EM/test4.tre")
        msa, site_names = read_sequences(msa_path,filetype="charMtrx",delimiter=",",masked_symbol='?')
        k = len(msa[next(iter(msa.keys()))])
        Q = read_priors(prior_path,msa,site_names=site_names)
        true_nllh = 1007.009158767873 
        true_phi = 0
        true_nu = 0

        parser_tree_out = parse_tree(T)
        out_phylogeny = parse_data(parser_tree_out, {'charMtrx':msa}, {'Q': Q}) # obj used to init the EMoptimizer
        ordered_leaf_labels = parser_tree_out['relabeling_vector'][:parser_tree_out['num_leaves']]
        charMtrx_recode = out_phylogeny.character_matrix
        ordered_leaf_labels = ordered_leaf_labels
        Q_recode = out_phylogeny.mutation_priors

        mySolver = fastEM_solver([T],{'charMtrx':msa, 'ordered_leaf_labels': [ordered_leaf_labels], 'charMtrx_recode': [charMtrx_recode]},{'Q':Q, "Q_recode": [Q_recode]},{'phi':0,'nu':0})

        nllh,status = mySolver.score_tree(strategy={'ultra_constr':False,'fixed_phi':None,'fixed_nu':None,'fixed_brlen':None,'nodes_to_recompute':None})
        nllh = - nllh
        phi = mySolver.params.phi
        nu = mySolver.params.nu
        threshold = 10
        abs_diff = abs(true_nllh - nllh)
        self.assertLess(abs_diff, threshold, msg="EMTest: test_47 failed.")
        self.assertAlmostEqual(true_phi,phi,places=4,msg="EMTest: test_47 failed.")
        self.assertAlmostEqual(true_nu,nu,places=4,msg="EMTest: test_47 failed.")
    # test score_tree
    def test_48(self):
        treedata_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test4.tre')
        msa_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test4_charMtrx.txt')
        prior_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test4_prior.csv')
        T = read_tree_newick(treedata_path)
        msa,_ = read_sequences(msa_path,filetype="charMtrx",delimiter=",",masked_symbol='?')
        k = len(msa[next(iter(msa.keys()))])
        Q = [{0:0} for i in range(k)]
        seen_sites = set()
        with open(prior_path, 'r') as fin:
            lines = fin.readlines()
            for line in lines:
                site_idx,char_state,prob = line.strip().split(',')
                site_idx = int(site_idx)
                if site_idx not in seen_sites:
                    seen_sites.add(site_idx)
                char_state = int(char_state)
                prob = float(prob)
                Q[len(seen_sites) - 1][char_state] = prob
        true_nllh = 1007.009158767873
        true_phi = 0
        true_nu = 0
        
        parser_tree_out = parse_tree(T)
        out_phylogeny = parse_data(parser_tree_out, {'charMtrx':msa}, {'Q': Q}) # obj used to init the EMoptimizer
        ordered_leaf_labels = parser_tree_out['relabeling_vector'][:parser_tree_out['num_leaves']]
        charMtrx_recode = out_phylogeny.character_matrix
        ordered_leaf_labels = ordered_leaf_labels
        Q_recode = out_phylogeny.mutation_priors

        mySolver = fastEM_solver([T],{'charMtrx':msa, 'ordered_leaf_labels': [ordered_leaf_labels], 'charMtrx_recode': [charMtrx_recode]},{'Q':Q, "Q_recode": [Q_recode]},{'phi':0,'nu':0})

        my_strategy = deepcopy(DEFAULT_STRATEGY) 
        my_strategy['fixed_brlen'] = None
        score,status = mySolver.score_tree(strategy=my_strategy)
        nllh = -score
        phi = mySolver.params.phi
        nu = mySolver.params.nu
        threshold = 10
        abs_diff = abs(true_nllh - nllh)
        self.assertLess(abs_diff, threshold, msg="EMTest: test_48 failed.")
        self.assertAlmostEqual(true_phi,phi,places=4,msg="EMTest: test_48 failed.")
        self.assertAlmostEqual(true_nu,nu,places=4,msg="EMTest: test_48 failed.")


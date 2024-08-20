import unittest
from laml_libs.IO_handler.sequence_lib import read_sequences
from laml_libs.Count_model.PMM_base import PMM_model, Alphabet,AlleleTable
from treeswift import *
from math import log
from random import random
from laml_libs import DEFAULT_STRATEGY
from copy import deepcopy
import pkg_resources
from .utils import *

class VirtualUnitTest(unittest.TestCase):
    def get_reduced_trees(self,tree_str):
        tree_obj = read_tree_newick(tree_str)
        tree_reduced = {}
        for node in tree_obj.traverse_leaves():
            tree_reduced[node.label] = tree_str
        node_list = [node for node in tree_obj.traverse_postorder() if not node.is_root() and not node.is_leaf()]
        for node in node_list:
            u,v = node.children
            u_len = u.edge_length
            v_len = v.edge_length
            node.remove_child(u)
            node.remove_child(v)
            tree_reduced[node.label] = tree_obj.newick()
            node.add_child(u)
            node.add_child(v)
            u.edge_length = u_len
            v.edge_length = v_len
        return tree_reduced    

    # computing the outside likelihood
    def check_outllh(self,T,Q,charMtrx,test_no=0,give_label=False,**params):
        # generic function to test the computed out0 and out1 after calling Estep_out_llh
        K = len(list(charMtrx.values())[0])
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
        allele_table = charMtrx_2_alleleTable(charMtrx,alphabet)
        myModel = PMM_model([T],{'alleleTable':allele_table},{'Q':Q},**params)
        if give_label:
            currIdx = 0
            for node in myModel.trees[0].traverse_preorder():
                if not node.is_leaf():
                    node.label = "I" + str(currIdx)
                    currIdx += 1
        myModel.Estep_in_llh()
        myModel.Estep_out_llh()
        out_llh = {}
        for node in myModel.trees[0].traverse_postorder():
            out_llh[node.label] = node.out_llh
        tree_reduced = self.get_reduced_trees(myModel.trees[0].newick())
        for x in tree_reduced:    
            charMtrx0 = {y:charMtrx[y] for y in charMtrx}
            charMtrx0[x]= [0]*K
            allele_table0 = charMtrx_2_alleleTable(charMtrx0,alphabet)
            tree_str = tree_reduced[x]
            myModel0 = PMM_model([tree_str],{'alleleTable':allele_table0},{'Q':Q},**params)
            myModel0.Estep_in_llh()
            for k in range(K):
                true = myModel0.trees[0].root.in_llh[k][tuple([0])]
                est = out_llh[x][k][tuple([0])]
                self.assertAlmostEqual(true,est+log(1-params['phi']),places=5,msg="PMMTest out llh: test_" + str(test_no) + " failed.")                

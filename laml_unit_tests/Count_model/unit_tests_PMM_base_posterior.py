import os 
import unittest
from laml_libs.Count_model.PMM_base import *
from treeswift import *
from laml_libs.IO_handler.sequence_lib import read_sequences
from random import random
#from .utils import *
from laml_libs.Count_model.utils import *
from .virtual_unit_tests import VirtualUnitTest
from math import *

class PMMTest_posterior(VirtualUnitTest):
    # test in_llh computation
    def test_1(self): 
        Q = [[{1:1}]]
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1]]])
        counts_a = countgen([(0,),(1,)],(1,))
        counts_b = countgen([(0,),(1,)],(1,))
        counts_c = countgen([(0,),(1,)],(1,))
        allele_table = AlleleTable(K,J,{'a':[counts_a],'b':[counts_b],'c':[counts_c]},alphabet)
        
        myModel = PMM_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,nu=0,phi=0)
        myModel.Estep()

        true_node_post0 = {'abc':-1.5016140490560352,'ab':-3.16936964867503}
        true_edge_post00 = {'abc':-1.5016140490560352,'ab':-3.16936964867503}
        true_edge_post01 = {'abc': -0.2520193579008661,'c':-1.5016140490560352,'ab':-1.710694503287948,'b':-3.16936964867503,'a':-3.16936964867503}
        true_edge_post11 = {'c':-0.2520193579008661,'ab':-0.2520193579008661,'b':-0.04293890366895339,'a':-0.04293890366895339}

        for node in myModel.trees[0].traverse_preorder():
            # test node_post0
            if node.label in true_node_post0:
                true = true_node_post0[node.label]
                est = node.log_node_posterior[0][(0,)]
                self.assertAlmostEqual(true,est,places=5,msg="PMMTest posterior: test_1 failed.")
            # test edge_post00
            if node.label in true_edge_post00:
                true = true_edge_post00[node.label]
                est = node.log_edge_posterior[0][((0,),(0,))]
                self.assertAlmostEqual(true,est,places=5,msg="PMMTest posterior: test_1 failed.")
            # test edge_post01
            if node.label in true_edge_post01:
                true = true_edge_post01[node.label]
                est = node.log_edge_posterior[0][((0,),(1,))]
                self.assertAlmostEqual(true,est,places=5,msg="PMMTest posterior: test_1 failed.")
            # test edge_post11
            if node.label in true_edge_post11:
                true = true_edge_post11[node.label]
                est = node.log_edge_posterior[0][((1,),(1,))]
                self.assertAlmostEqual(true,est,places=5,msg="PMMTest posterior: test_1 failed.")
    '''
    def test_2(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        
        K = len(list(charMtrx.values())[0])
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])        
        allele_table = charMtrx_2_alleleTable(charMtrx,alphabet)

        T = '((((b:1,c:1)bc:1,d:1)bcd:1,a:1)abcd:1)r;'                
            
        myModel = PMM_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,phi=0,nu=0)
        myModel.Estep()

        for node in myModel.trees[0].traverse_preorder():
            print(node.label,node.log_edge_posterior) '''

import os 
import unittest
from laml_libs.Count_model.PMMN_model import *
from treeswift import *
from laml_libs.IO_handler.sequence_lib import read_sequences
from random import random
from laml_libs.Count_model.utils import *
from .virtual_unit_tests import VirtualUnitTest
from math import *
from laml_libs.Count_model.CharMtrx import CharMtrx

class PMMNTest_posterior(VirtualUnitTest):
    # test in_llh computation
    def test_1(self): 
        Q = [[{1:1}]]
        charMtrx = {'a':[1],'b':[1],'c':[1]}
        T = "(((a:1,b:1)ab:1,c:1)abc:1)r;"
        K = 1
        J = 1
        alphabet = Alphabet(K,J,[[[0,1,-1]]])
        DLT_data = CharMtrx(charMtrx,alphabet)
        
        myModel = PMMN_model([T],{'DLT_data':DLT_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
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
                self.assertAlmostEqual(true,est,places=5,msg="PMMNTest posterior: test_1 failed.")
            # test edge_post00
            if node.label in true_edge_post00:
                true = true_edge_post00[node.label]
                est = node.log_edge_posterior[0][((0,),(0,))]
                self.assertAlmostEqual(true,est,places=5,msg="PMMNTest posterior: test_1 failed.")
            # test edge_post01
            if node.label in true_edge_post01:
                true = true_edge_post01[node.label]
                est = node.log_edge_posterior[0][((0,),(1,))]
                self.assertAlmostEqual(true,est,places=5,msg="PMMNTest posterior: test_1 failed.")
            # test edge_post11
            if node.label in true_edge_post11:
                true = true_edge_post11[node.label]
                est = node.log_edge_posterior[0][((1,),(1,))]
                self.assertAlmostEqual(true,est,places=5,msg="PMMNTest posterior: test_1 failed.")

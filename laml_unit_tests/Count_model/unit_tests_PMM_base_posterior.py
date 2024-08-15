import os 
import unittest
from laml_libs.Count_model.PMM_base import *
from treeswift import *
from laml_libs.IO_handler.sequence_lib import read_sequences
from random import random
from .utils import *
from .virtual_unit_tests import VirtualUnitTest

class PMM_Test_posterior(VirtualUnitTest):
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

        true_post0 = {'abc':-1.5016140490560352,'ab':-3.16936964867503}
        for node in myModel.trees[0].traverse_preorder():
            if node.label in true_post0:
                true = true_post0[node.label]
                est = node.log_posterior[0][(0,)]
                self.assertAlmostEqual(true,est,places=5,msg="PMMTest posterior: test_1 failed.")

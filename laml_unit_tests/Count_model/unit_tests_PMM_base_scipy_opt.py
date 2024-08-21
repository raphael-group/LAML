import os 
import unittest
from laml_libs.Count_model.PMM_base import *
from treeswift import *
from laml_libs.IO_handler.sequence_lib import read_sequences
from random import random
#from laml_libs.Count_model.utils import *
from .utils import *
from .virtual_unit_tests import VirtualUnitTest

class PMMTest_scipy_opt(VirtualUnitTest):
    # test optimization using scipy 
    def test_1(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        
        K = len(list(charMtrx.values())[0])
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])        
        allele_table = charMtrx_2_alleleTable(charMtrx,alphabet)

        tree_list = ['((((b,c),d),a));', '(((b,c),(d,a)));', '((d,((b,c),a)));', '((d,(b,(c,a))));', '((d,(c,(b,a))));', '((((b,d),c),a));', '(((b,d),(c,a)));', '((c,((b,d),a)));', '((c,(b,(d,a))));', '((c,(d,(b,a))));', '((((c,d),b),a));', '(((c,d),(b,a)));', '((b,((c,d),a)));', '((b,(c,(d,a))));', '((b,(d,(c,a))));']  
        true_nllh_list = [11.809140931727208, 11.809141253319298, 11.80914097392261, 11.809140974410134, 10.338626804278578, 11.809141098008908, 11.809140926336006, 11.809141047672332, 11.809141363494154, 10.33862520405755, 9.322029697571756, 7.851513459377366, 9.32202949252424, 11.809140996738018, 11.809140939639644]
        
        for T,true_nllh in zip(tree_list,true_nllh_list):
            myModel = PMM_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,phi=0,nu=0)
            test_nllh,_ = myModel.optimize('Scipy',initials=1,verbose=-1,ultra_constr=False)
            self.assertAlmostEqual(true_nllh,test_nllh,places=4,msg="PMMTest Scipy-opt: test_1 failed.")

    def test_2(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        
        K = len(list(charMtrx.values())[0])
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])        
        allele_table = charMtrx_2_alleleTable(charMtrx,alphabet)

        tree_list = ['((((b,c),d),a));', '(((b,c),(d,a)));', '((d,((b,c),a)));', '((d,(b,(c,a))));', '((d,(c,(b,a))));', '((((b,d),c),a));', '(((b,d),(c,a)));', '((c,((b,d),a)));', '((c,(b,(d,a))));', '((c,(d,(b,a))));', '((((c,d),b),a));', '(((c,d),(b,a)));', '((b,((c,d),a)));', '((b,(c,(d,a))));', '((b,(d,(c,a))));']  
        true_nllh_list = [7.595936888280069, 5.078900745505063, 7.595936937640891, 10.0830486451717, 10.083048625619265, 10.08305235692635, 10.083048544279704, 10.083048478217442, 7.566011711889361, 10.083048503399493, 10.083048490449752, 10.083048551149702, 10.083048664742668, 7.566011487530715, 10.083048480708316]
        
        for T,true_nllh in zip(tree_list,true_nllh_list):
            myModel = PMM_model([T],{'alleleTable':allele_table},{'Q':Q},mu=1,phi=0,nu=0)
            test_nllh,_ = myModel.optimize('Scipy',initials=1,verbose=-1,ultra_constr=False)
            self.assertAlmostEqual(true_nllh,test_nllh,places=4,msg="PMMTest Scipy-opt: test_2 failed.")

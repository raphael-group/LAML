import os 
import unittest
from laml_libs.Count_model.PMMN_model import PMMN_model
from laml_libs.Count_model.CharMtrx import CharMtrx
from treeswift import *
from laml_libs.IO_handler.sequence_lib import read_sequences
from random import random
from .utils import *
from .virtual_unit_tests import VirtualUnitTest
import pkg_resources
from laml_libs import DEFAULT_STRATEGY
from copy import deepcopy

class PMMNTest_EM_opt(VirtualUnitTest):
    # test optimization using scipy 
    def test_1(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        
        K = len(list(charMtrx.values())[0])
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])        
        charMtrx = CharMtrx(charMtrx,alphabet)

        tree_list = ['((((b,c),d),a));', '(((b,c),(d,a)));', '((d,((b,c),a)));', '((d,(b,(c,a))));', '((d,(c,(b,a))));', '((((b,d),c),a));', '(((b,d),(c,a)));', '((c,((b,d),a)));', '((c,(b,(d,a))));', '((c,(d,(b,a))));', '((((c,d),b),a));', '(((c,d),(b,a)));', '((b,((c,d),a)));', '((b,(c,(d,a))));', '((b,(d,(c,a))));']  
        true_nllh_list = [11.809140931727208, 11.809141253319298, 11.80914097392261, 11.809140974410134, 10.338626804278578, 11.809141098008908, 11.809140926336006, 11.809141047672332, 11.809141363494154, 10.33862520405755, 9.322029697571756, 7.851513459377366, 9.32202949252424, 11.809140996738018, 11.809140939639644]
        
        for T,true_nllh in zip(tree_list,true_nllh_list):
            myModel = PMMN_model([T],{'DLT_data':charMtrx},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
            test_nllh,_ = myModel.optimize('EM',initials=1,verbose=-1,ultra_constr=False,fixed_params={'rho':1})
            self.assertAlmostEqual(true_nllh,test_nllh,places=4,msg="PMMNTest EM-opt: test_1 failed.")

    def test_2(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        
        K = len(list(charMtrx.values())[0])
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])        
        charMtrx = CharMtrx(charMtrx,alphabet)

        tree_list = ['((((b,c),d),a));', '(((b,c),(d,a)));', '((d,((b,c),a)));', '((d,(b,(c,a))));', '((d,(c,(b,a))));', '((((b,d),c),a));', '(((b,d),(c,a)));', '((c,((b,d),a)));', '((c,(b,(d,a))));', '((c,(d,(b,a))));', '((((c,d),b),a));', '(((c,d),(b,a)));', '((b,((c,d),a)));', '((b,(c,(d,a))));', '((b,(d,(c,a))));']  
        true_nllh_list = [7.595936888280069, 5.078900745505063, 7.595936937640891, 10.0830486451717, 10.083048625619265, 10.08305235692635, 10.083048544279704, 10.083048478217442, 7.566011711889361, 10.083048503399493, 10.083048490449752, 10.083048551149702, 10.083048664742668, 7.566011487530715, 10.083048480708316]
        
        for T,true_nllh in zip(tree_list,true_nllh_list):
            myModel = PMMN_model([T],{'DLT_data':charMtrx},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
            test_nllh,_ = myModel.optimize('EM',initials=1,verbose=-1,ultra_constr=False,fixed_params={'rho':1})
            self.assertAlmostEqual(true_nllh,test_nllh,places=4,msg="PMMNTest EM-opt: test_2 failed.")
    
    def test_3(self):
        treedata_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test4.tre')
        charMtrx_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test4_charMtrx.txt')
        prior_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test4_prior.csv')
        
        T = read_tree_newick(treedata_path)
        charMtrx,_ = read_sequences(charMtrx_path,filetype="charMtrx",delimiter=",",masked_symbol='?')
        
        K = len(charMtrx[next(iter(charMtrx.keys()))])
        J = 1
        
        Q = []
        for k in range(K):
            Q_k = [{}]
            Q.append(Q_k)
        seen_cassette = set()
        
        with open(prior_path,'r') as fin:
            lines = fin.readlines()
            for line in lines:
                cassette_idx,char_state,prob = line.strip().split(',')
                cassette_idx = int(cassette_idx)
                if cassette_idx not in seen_cassette:
                    seen_cassette.add(cassette_idx)
                char_state = int(char_state)
                prob = float(prob)
                Q[len(seen_cassette) - 1][0][char_state] = prob
        
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
        charMtrx = CharMtrx(charMtrx,alphabet)

        true_nllh = 1007.009158767873 
        true_phi = 0
        true_nu = 0

        myModel = PMMN_model([T],{'DLT_data':charMtrx},{'Q':Q})
        randseed = 1221
        nllh,status = myModel.optimize('EM',initials=1,random_seeds=randseed,verbose=-1,ultra_constr=False,fixed_params={'rho':1})
        phi = myModel.params.get_value('phi')
        nu = myModel.params.get_value('nu')
        self.assertAlmostEqual(0,abs(true_nllh-nllh)/true_nllh,places=4,msg="PMMNTest EM-opt: test_3 failed.")
        self.assertAlmostEqual(true_phi,phi,places=4,msg="PMMNTest EM-opt: test_3 failed.")
        self.assertAlmostEqual(true_nu,nu,places=4,msg="PMMNTest EM-opt: test_3 failed.")
    
    def test_4(self):
        treedata_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test4.tre')
        charMtrx_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test4_charMtrx.txt')
        prior_path = pkg_resources.resource_filename('laml_unit_tests', 'test_data/test_EM/test4_prior.csv')
        
        T = read_tree_newick(treedata_path)
        charMtrx,_ = read_sequences(charMtrx_path,filetype="charMtrx",delimiter=",",masked_symbol='?')
        
        K = len(charMtrx[next(iter(charMtrx.keys()))])
        J = 1
        
        Q = []
        for k in range(K):
            Q_k = [{}]
            Q.append(Q_k)
        seen_cassette = set()
        
        with open(prior_path,'r') as fin:
            lines = fin.readlines()
            for line in lines:
                cassette_idx,char_state,prob = line.strip().split(',')
                cassette_idx = int(cassette_idx)
                if cassette_idx not in seen_cassette:
                    seen_cassette.add(cassette_idx)
                char_state = int(char_state)
                prob = float(prob)
                Q[len(seen_cassette) - 1][0][char_state] = prob
        
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
        charMtrx = CharMtrx(charMtrx,alphabet)

        true_nllh = 1007.009158767873 
        true_phi = 0
        true_nu = 0

        myModel = PMMN_model([T],{'DLT_data':charMtrx},{'Q':Q})
        randseed = 1221
        my_strategy = deepcopy(DEFAULT_STRATEGY) 
        my_strategy['fixed_brlen'] = None
        my_strategy['fixed_params'] = {'rho':1}
        score,status = myModel.score_tree(my_strategy)
        nllh = -score
        phi = myModel.params.get_value('phi')
        nu = myModel.params.get_value('nu')
        self.assertAlmostEqual(0,abs(true_nllh-nllh)/true_nllh,places=4,msg="PMMNTest EM-opt: test_4 failed.")
        self.assertAlmostEqual(true_phi,phi,places=4,msg="PMMNTest EM-opt: test_4 failed.")
        self.assertAlmostEqual(true_nu,nu,places=4,msg="PMMNTest EM-opt: test_4 failed.")

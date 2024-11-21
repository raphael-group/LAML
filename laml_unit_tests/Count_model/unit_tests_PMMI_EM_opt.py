import numpy as np 
import os 
import unittest
from laml_libs.Count_model.PMMI_model import PMMI_model
from laml_libs.Count_model.IntensityData import IntensityData
from treeswift import *
from laml_libs.IO_handler.sequence_lib import read_sequences
from random import random
from .utils import *
from .virtual_unit_tests import VirtualUnitTest
import pkg_resources
from laml_libs import DEFAULT_STRATEGY
from copy import deepcopy

class PMMITest_EM_opt(VirtualUnitTest):
    # test optimization using scipy 
    def test_1(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        dim = 2

        means = [[0]*dim, [5]*dim]
        covariance_matrix = [[1 if i==y else 0 for i in range(dim)] for y in range(dim)]
        training_data = {}
        num_samples = 1000
        for i in range(2):
            training_data[i] = list(np.random.multivariate_normal(means[0], covariance_matrix, num_samples))
        intensityData = {x:[means[y] for y in charMtrx[x]] for x in charMtrx}
        #intensityData = {x:[list(np.random.multivariate_normal(means[y], covariance_matrix, 1)[0]) for y in charMtrx[x]] for x in charMtrx}

        K = len(list(intensityData.values())[0])
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])        
        intensityData = IntensityData(intensityData,alphabet)

        tree_list = ['((((b,c),d),a));', '(((b,c),(d,a)));', '((d,((b,c),a)));', '((d,(b,(c,a))));', '((d,(c,(b,a))));', '((((b,d),c),a));', '(((b,d),(c,a)));', '((c,((b,d),a)));', '((c,(b,(d,a))));', '((c,(d,(b,a))));', '((((c,d),b),a));', '(((c,d),(b,a)));', '((b,((c,d),a)));', '((b,(c,(d,a))));', '((b,(d,(c,a))));']  
        true_nllh_list = [11.809140931727208, 11.809141253319298, 11.80914097392261, 11.809140974410134, 10.338626804278578, 11.809141098008908, 11.809140926336006, 11.809141047672332, 11.809141363494154, 10.33862520405755, 9.322029697571756, 7.851513459377366, 9.32202949252424, 11.809140996738018, 11.809140939639644]
       
        # check that the relative ordering of the trees by likelihood is the same
        true_sorted_indices = np.argsort(true_nllh_list)
        print(true_sorted_indices)
        print(true_nllh_list)
        test_nllh_list = []
        for T,true_nllh in zip(tree_list,true_nllh_list):
            myModel = PMMI_model([T],{'DLT_data':intensityData, 'training_data': training_data},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
            test_nllh,_ = myModel.optimize('EM',initials=1,verbose=-1,ultra_constr=False,fixed_params={'rho':1})
            test_nllh_list.append(test_nllh)
            #self.assertAlmostEqual(true_nllh,test_nllh,places=4,msg="PMMITest EM-opt: test_1 failed.")
        test_sorted_indices = np.argsort(test_nllh_list)
        print(test_sorted_indices)
        print(test_nllh_list)

        order_matches = np.array_equal(true_sorted_indices, test_sorted_indices)
        print(order_matches)
        self.assertTrue(order_matches, "PMMITest EM-opt: test_1 failed. The relative order of tree likelihoods does not match.")



    """
    def test_2(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        intensityData = {x:[{y:1} for y in charMtrx[x]] for x in charMtrx}

        K = len(list(intensityData.values())[0])
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])        
        intensityData = IntensityData(intensityData,alphabet)

        tree_list = ['((((b,c),d),a));', '(((b,c),(d,a)));', '((d,((b,c),a)));', '((d,(b,(c,a))));', '((d,(c,(b,a))));', '((((b,d),c),a));', '(((b,d),(c,a)));', '((c,((b,d),a)));', '((c,(b,(d,a))));', '((c,(d,(b,a))));', '((((c,d),b),a));', '(((c,d),(b,a)));', '((b,((c,d),a)));', '((b,(c,(d,a))));', '((b,(d,(c,a))));']  
        true_nllh_list = [7.595936888280069, 5.078900745505063, 7.595936937640891, 10.0830486451717, 10.083048625619265, 10.08305235692635, 10.083048544279704, 10.083048478217442, 7.566011711889361, 10.083048503399493, 10.083048490449752, 10.083048551149702, 10.083048664742668, 7.566011487530715, 10.083048480708316]
        
        for T,true_nllh in zip(tree_list,true_nllh_list):
            myModel = PMMI_model([T],{'DLT_data':intensityData},{'Q':Q},{'mu':1,'nu':0,'phi':0,'rho':1})
            test_nllh,_ = myModel.optimize('EM',initials=1,verbose=-1,ultra_constr=False,fixed_params={'rho':1})
            self.assertAlmostEqual(true_nllh,test_nllh,places=4,msg="PMMITest EM-opt: test_2 failed.")
    
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
        intensityData = {x:[{y:1} for y in charMtrx[x]] for x in charMtrx}
        DLT_data = IntensityData(intensityData,alphabet)

        true_nllh = 1007.009158767873 
        true_phi = 0
        true_nu = 0

        myModel = PMMI_model([T],{'DLT_data':DLT_data},{'Q':Q})
        randseed = 1221
        nllh,status = myModel.optimize('EM',initials=1,random_seeds=randseed,verbose=-1,ultra_constr=False,fixed_params={'rho':1})
        phi = myModel.params.get_value('phi')
        nu = myModel.params.get_value('nu')
        self.assertAlmostEqual(0,abs(true_nllh-nllh)/true_nllh,places=4,msg="PMMITest EM-opt: test_3 failed.")
        self.assertAlmostEqual(true_phi,phi,places=4,msg="PMMITest EM-opt: test_3 failed.")
        self.assertAlmostEqual(true_nu,nu,places=4,msg="PMMITest EM-opt: test_3 failed.")
    
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
        intensityData = {x:[{y:1} for y in charMtrx[x]] for x in charMtrx}
        DLT_data = IntensityData(intensityData,alphabet)

        true_nllh = 1007.009158767873 
        true_phi = 0
        true_nu = 0

        myModel = PMMI_model([T],{'DLT_data':DLT_data},{'Q':Q})
        randseed = 1221
        my_strategy = deepcopy(DEFAULT_STRATEGY) 
        my_strategy['fixed_brlen'] = None
        my_strategy['fixed_params'] = {'rho':1}
        score,status = myModel.score_tree(my_strategy)
        nllh = -score
        phi = myModel.params.get_value('phi')
        nu = myModel.params.get_value('nu')
        self.assertAlmostEqual(0,abs(true_nllh-nllh)/true_nllh,places=4,msg="PMMITest EM-opt: test_4 failed.")
        self.assertAlmostEqual(true_phi,phi,places=4,msg="PMMITest EM-opt: test_4 failed.")
        self.assertAlmostEqual(true_nu,nu,places=4,msg="PMMITest EM-opt: test_4 failed.")
    """

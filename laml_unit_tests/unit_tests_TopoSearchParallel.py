import os 
import unittest
from scmail_libs import *
from scmail_libs.EM_solver import EM_solver
from scmail_libs.ML_solver import ML_solver
from scmail_libs.Topology_search_parallel import Topology_search_parallel as Topology_search
from treeswift import *
from copy import deepcopy

class TopoSearchParallelTest(unittest.TestCase):
    # topology search with EM_solver
    def test_1(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        #best_nllh,best_tree = self.__brute_force_search__(msa,Q,['a','b','c','d'],solver=ML_solver)
        nllh_bf = 7.877269958604131 # precomputed from brute-force
         
        T0 = '((a,c),(b,d));'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        # topology search with EM_solver
        myTopoSearch_EM = Topology_search([T0],EM_solver,data=data,prior=prior,params=params)
        best_tree,max_score,best_params = myTopoSearch_EM.search(maxiter=200,nreps=1,verbose=False)
        nllh_nni_EM = -max_score
        
        self.assertAlmostEqual(nllh_bf,nllh_nni_EM,places=4,msg="TopoSearchTest: test_1 failed.")
    
    # topology search with ML_solver
    def test_2(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        nllh_bf = 7.877269958604131 # precomputed from brute-force
        
        T0 = '((a,c),(b,d));'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'nu':0,'phi':0}
        
        # topology search with ML_solver
        myTopoSearch_ML = Topology_search([T0],ML_solver,data=data,prior=prior,params=params)
        best_tree,max_score,best_params = myTopoSearch_ML.search(maxiter=200,verbose=False,nreps=1)
        nllh_nni_ML = -max_score
        
        self.assertAlmostEqual(nllh_bf,nllh_nni_ML,places=4,msg="TopoSearchTest: test_2 failed.")
    
    # topology search on a tree with polytomies
    def test_3(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        nllh_bf = 7.877269958604131 # precomputed from brute-force
        
        T0 = '(a,b,c,d);'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'nu':0,'phi':0}
        
        myTopoSearch = Topology_search([T0],EM_solver,data=data,prior=prior,params=params)
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,nreps=1)
        nllh_nni = -max_score
        
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_3 failed.")

    # only resolve polytomies
    def test_4(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        nllh_bf = 7.877269958604131 # precomputed from brute-force
        
        T0 = '((a,b),c,d);'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'nu':0,'phi':0}
        
        myTopoSearch = Topology_search([T0],EM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['resolve_search_only'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)
        nllh_nni = -max_score
        
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_4 failed.")
    
    # resolve polytomies set to true on starting tree without polytomies
    def test_5(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        #best_nllh,best_tree = self.__brute_force_search__(msa,Q,['a','b','c','d'],solver=ML_solver)
        nllh_bf = 12.48204978089449 # precomputed from brute-force
        
        T0 = '((a,c),(b,d));'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        # topology search with EM_solver
        myTopoSearch = Topology_search([T0],EM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['resolve_search_only'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)
        nllh_nni = -max_score
        
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_5 failed.")
    
    # full topology search
    def test_6(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        #best_nllh,best_tree = self.__brute_force_search__(msa,Q,['a','b','c','d'],solver=ML_solver)
         
        nllh_bf = 6.042549046654788
        
        T0 = '((a,b),(c,d));'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search([T0],EM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        #my_strategy['resolve_search_only'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)

        nllh_nni = -max_score

        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_6 failed.") 
        
    # resolve polytomies only
    def test_7(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        #best_nllh,best_tree = self.__brute_force_search__(msa,Q,['a','b','c','d'],solver=ML_solver)
        nllh_bf = 12.010496077449595
        
        T0 = '((a,b),c,d);'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search([T0],EM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['resolve_search_only'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)
        
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_7 failed.")
    
    # resolve polytomies only
    def test_8(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf = 12.010496077449595
        
        T0 = '((c,d),a,b);'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search([T0],EM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['resolve_search_only'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_8 failed.")
    
    # enforce ultrametric
    def test_9(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        #best_nllh,best_tree = self.__brute_force_search__(msa,Q,['a','b','c','d'],solver=ML_solver,ultra_constr=True)
        nllh_bf = 7.0063474330891236 # pre-computed using brute-force search
        
        T0 = '((a,b),(c,d));'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search([T0],EM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['ultra_constr'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_9 failed.")
    
    # enforce ultrametric and starting tree is a star-tree
    def test_10(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf = 7.0063474330891236 # pre-computed using brute-force search
        T0 = '(a,b,c,d);'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search([T0],EM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['ultra_constr'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)

        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_10 failed.")

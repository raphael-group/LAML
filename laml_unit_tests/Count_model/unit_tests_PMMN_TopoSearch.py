import os 
import unittest
from laml_libs import *
from laml_libs.PMM_original.EM_solver import EM_solver
from laml_libs.PMM_original.ML_solver import ML_solver
from laml_libs.Count_model.PMMN_model import PMMN_model
from laml_libs.Count_model.CharMtrx import CharMtrx
from laml_libs.Count_model.Alphabet import Alphabet
from laml_libs.TopoSearch.Topology_search import Topology_search
from treeswift import *
from copy import deepcopy
from .virtual_unit_tests import VirtualUnitTest

class PMMNTest_TopoSearch(VirtualUnitTest):
    def __list_topologies(self,leafset):
        def __triplet(a,b,c):
            return "((("+a+","+b+"),"+c+"));"
        def __add_one(tree,new_leaf):    
            new_topos = []
            nodes = [node for node in tree.traverse_preorder() if not node.is_root()]
            for node in nodes:
                # create new topology
                new_node1 = Node()
                new_node2 = Node()
                new_node2.label = new_leaf
                pnode = node.parent
                pnode.remove_child(node)
                pnode.add_child(new_node1)
                new_node1.add_child(node)
                new_node1.add_child(new_node2)
                new_topos.append(tree.newick())
                # turn back
                pnode.remove_child(new_node1)
                pnode.add_child(node)
            return new_topos            

        # initialization
        a,b,c = leafset[-3:]
        T1 = __triplet(a,b,c)
        T2 = __triplet(a,c,b)
        T3 = __triplet(b,c,a)
        topos = [T1,T2,T3]
        L = leafset[:-3]
        # elaborate
        while L:
            new_leaf = L.pop()
            new_topos = []
            for T in topos:
                tree = read_tree_newick(T)
                new_topos += __add_one(tree,new_leaf)
            topos = new_topos    
        out_topos = [topo[1:-2]+";" for topo in topos]
        return out_topos    

    def __brute_force_search(self,charMtrx,Q,L,solver=EM_solver,ultra_constr=False,initials=1):
        topos = self.__list_topologies(L)
        best_nllh = float("inf")
        best_tree = ""
        for T in topos:
            mySolver = solver([T],{'charMtrx':charMtrx},{'Q':Q})
            nllh,_ = mySolver.optimize(initials=initials,verbose=-1,ultra_constr=ultra_constr)
            if nllh < best_nllh:
                best_nllh = nllh
                best_tree = mySolver.trees[0].newick()
        return best_nllh,best_tree 
    
    def test_1(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        nllh_bf,bf_tree = self.__brute_force_search(charMtrx,[Q_k[0] for Q_k in Q],['a','b','c','d'],solver=ML_solver)
        
        K = 5
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
         
        T0 = '((a,c),(b,d));'
        data = {'DLT_data':CharMtrx(charMtrx,alphabet)}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0,'mu':1,'rho':1}
        
        myTopoSearch = Topology_search([T0],PMMN_model,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY) 
        my_strategy['fixed_params'] = {'rho':1}
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,nreps=1,verbose=False,strategy=my_strategy)
        nllh_nni = -max_score
        
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="PMMNTest_TopoSearch: test_1 failed.")
    
    # topology search on a tree with polytomies
    def test_2(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        nllh_bf,bf_tree = self.__brute_force_search(charMtrx,[Q_k[0] for Q_k in Q],['a','b','c','d'],solver=ML_solver)
        
        K = 5
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
        
        T0 = '(a,b,c,d);'
        data = {'DLT_data':CharMtrx(charMtrx,alphabet)}
        prior = {'Q':Q}
        params = {'nu':0,'phi':0,'mu':1,'rho':1}
        
        myTopoSearch = Topology_search([T0],PMMN_model,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY) 
        my_strategy['fixed_params'] = {'rho':1}
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,nreps=1,strategy=my_strategy)
        nllh_nni = -max_score
        
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="PMMNTest_TopoSearch: test_2 failed.")

    # only resolve polytomies
    def test_3(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        nllh_bf,bf_tree = self.__brute_force_search(charMtrx,[Q_k[0] for Q_k in Q],['a','b','c','d'],solver=ML_solver)
        
        K = 5
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
        
        T0 = '((a,b),c,d);'
        data = {'DLT_data':CharMtrx(charMtrx,alphabet)}
        prior = {'Q':Q}
        params = {'nu':0,'phi':0,'mu':1,'rho':1}
        
        myTopoSearch = Topology_search([T0],PMMN_model,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['resolve_search_only'] = True
        my_strategy['fixed_params'] = {'rho':1}
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)
        nllh_nni = -max_score
        
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="PMMNTest_TopoSearch: test_3 failed.")
    
    # resolve polytomies set to true on starting tree without polytomies
    def test_4(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        nllh_bf = 11.809140958825493 # precomputed from brute-force
        
        K = 5
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
        
        T0 = '((a,c),(b,d));'
        data = {'DLT_data':CharMtrx(charMtrx,alphabet)}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0,'mu':1,'rho':1}
        
        # topology search with EM_solver
        myTopoSearch = Topology_search([T0],PMMN_model,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['resolve_search_only'] = True
        my_strategy['fixed_params'] = {'rho':1}
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="PMMNTest_TopoSearch: test_4 failed.")
    
    # full topology search
    def test_5(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf,bf_tree = self.__brute_force_search(charMtrx,[Q_k[0] for Q_k in Q],['a','b','c','d'],solver=ML_solver)
        
        K = 5
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
        
        T0 = '((a,b),(c,d));'
        data = {'DLT_data':CharMtrx(charMtrx,alphabet)}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0,'mu':1,'rho':1}
        
        myTopoSearch = Topology_search([T0],PMMN_model,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['fixed_params'] = {'rho':1}
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)
        nllh_nni = -max_score

        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="PMMNTest_TopoSearch: test_5 failed.") 
        
    # resolve polytomies only
    def test_6(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf = 10.083048489886435 # precomputed from brute-force
 
        K = 5
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
        
        T0 = '((a,b),c,d);'
        data = {'DLT_data':CharMtrx(charMtrx,alphabet)}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0,'mu':1,'rho':1}
        
        myTopoSearch = Topology_search([T0],PMMN_model,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['resolve_search_only'] = True
        my_strategy['fixed_params'] = {'rho':1}
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)

        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="PMMNTest_TopoSearch: test_6 failed.")
    
    # resolve polytomies only
    def test_7(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf = 10.08306053772235
        
        K = 5
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
        
        T0 = '((c,d),a,b);'
        data = {'DLT_data':CharMtrx(charMtrx,alphabet)}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0,'mu':1,'rho':1}
        
        myTopoSearch = Topology_search([T0],PMMN_model,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['resolve_search_only'] = True
        my_strategy['fixed_params'] = {'rho':1}
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="PMMNTest_TopoSearch: test_7 failed.")
    
    # enforce ultrametric
    def test_8(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf,bf_tree = self.__brute_force_search(charMtrx,[Q_k[0] for Q_k in Q],['a','b','c','d'],solver=ML_solver,ultra_constr=True)
        
        K = 5
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])

        T0 = '((a,b),(c,d));'
        data = {'DLT_data':CharMtrx(charMtrx,alphabet)}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0,'mu':1,'rho':1}
        
        myTopoSearch = Topology_search([T0],PMMN_model,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['ultra_constr'] = True
        my_strategy['fixed_params'] = {'rho':1}
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="PMMNTest_TopoSearch: test_8 failed.")
    
    # enforce ultrametric and starting tree is a star-tree
    def test_9(self):
        Q = [[{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}], [{1:1.0}]]
        charMtrx = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf,bf_tree = self.__brute_force_search(charMtrx,[Q_k[0] for Q_k in Q],['a','b','c','d'],solver=ML_solver,ultra_constr=True)
        
        K = 5
        J = 1
        alphabet = Alphabet(K,J,[[[0,-1]+list(Q[k][0].keys())] for k in range(K)])
        
        T0 = '(a,b,c,d);'
        data = {'DLT_data':CharMtrx(charMtrx,alphabet)}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0,'mu':1,'rho':1}
        
        myTopoSearch = Topology_search([T0],PMMN_model,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['ultra_constr'] = True
        my_strategy['fixed_params'] = {'rho':1}
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)

        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="PMMNTest_TopoSearch: test_9 failed.")

import os 
import unittest
from laml_libs import *
from laml_libs.fastEM_solver import fastEM_solver, parse_data, parse_tree
from laml_libs.EM_solver import EM_solver
from laml_libs.ML_solver import ML_solver
from laml_libs.Topology_search import Topology_search
from treeswift import *
from copy import deepcopy

class TopoSearchTest_wfastEMSolver(unittest.TestCase):
    def __list_topologies__(self,leafset):
        def __triplet__(a,b,c):
            return "((("+a+","+b+"),"+c+"));"
        def __add_one__(tree,new_leaf):    
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
        T1 = __triplet__(a,b,c)
        T2 = __triplet__(a,c,b)
        T3 = __triplet__(b,c,a)
        topos = [T1,T2,T3]
        L = leafset[:-3]
        # elaborate
        while L:
            new_leaf = L.pop()
            new_topos = []
            for T in topos:
                tree = read_tree_newick(T)
                new_topos += __add_one__(tree,new_leaf)
            topos = new_topos    
        out_topos = [topo[1:-2]+";" for topo in topos]
        return out_topos    

    def __brute_force_search__(self,msa,Q,L,solver=EM_solver,ultra_constr=False,initials=1):
        topos = self.__list_topologies__(L)
        best_nllh = float("inf")
        best_tree = ""
        for T in topos:
            mySolver = solver([T],{'charMtrx':msa},{'Q':Q})
            nllh,_ = mySolver.optimize(initials=initials,verbose=-1,ultra_constr=ultra_constr)
            #print(T,nllh)
            if nllh < best_nllh:
                best_nllh = nllh
                best_tree = mySolver.trees[0].newick()
        return best_nllh,best_tree 
    
    # topology search with fastEM_solver
    def test_1(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        nllh_bf,bf_tree = self.__brute_force_search__(msa,Q,['a','b','c','d'],solver=ML_solver)
         

        T0 = '((a,c),(b,d));'
        
        parser_tree_out = parse_tree(read_tree_newick(T0))
        out_phylogeny = parse_data(parser_tree_out, {'charMtrx':msa}, {'Q': Q}) # obj used to init the EMoptimizer
        Q_recode = out_phylogeny.mutation_priors
        ordered_leaf_labels = parser_tree_out['relabeling_vector'][:parser_tree_out['num_leaves']]
        charMtrx_recode = out_phylogeny.character_matrix

        data = {'charMtrx':msa, 'ordered_leaf_labels': [ordered_leaf_labels], 'charMtrx_recode': [charMtrx_recode]}
        prior = {'Q':Q, 'Q_recode': [Q_recode]}
        params = {'phi':0,'nu':0}
       
        # topology search with EM_solver
        myTopoSearch_EM = Topology_search([T0],fastEM_solver,data=data,prior=prior,params=params)
        best_tree,max_score,best_params = myTopoSearch_EM.search(maxiter=200,nreps=10,verbose=False)
        nllh_nni_EM = -max_score
        
        self.assertAlmostEqual(nllh_bf,nllh_nni_EM,places=4,msg="TopoSearchTest: test_1 failed.")
  
    # topology search on a tree with polytomies
    def test_3(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        nllh_bf,bf_tree = self.__brute_force_search__(msa,Q,['a','b','c','d'],solver=ML_solver)
        
        T0 = '(a,b,c,d);'
        
        parser_tree_out = parse_tree(read_tree_newick(T0))
        out_phylogeny = parse_data(parser_tree_out, {'charMtrx':msa}, {'Q': Q}) # obj used to init the EMoptimizer
        Q_recode = out_phylogeny.mutation_priors
        ordered_leaf_labels = parser_tree_out['relabeling_vector'][:parser_tree_out['num_leaves']]
        charMtrx_recode = out_phylogeny.character_matrix

        data = {'charMtrx':msa, 'ordered_leaf_labels': [ordered_leaf_labels], 'charMtrx_recode': [charMtrx_recode]}
        prior = {'Q':Q, 'Q_recode': [Q_recode]}
        params = {'phi':0,'nu':0}
       
        myTopoSearch = Topology_search([T0],fastEM_solver,data=data,prior=prior,params=params)
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,nreps=10)
        nllh_nni = -max_score
        
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_3 failed.")
    
    # only resolve polytomies [NOT IMPLEMENTED]
    """
    def test_4(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        nllh_bf,bf_tree = self.__brute_force_search__(msa,Q,['a','b','c','d'],solver=ML_solver)
        
        T0 = '((a,b),c,d);'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'nu':0,'phi':0}
        
        myTopoSearch = Topology_search([T0],fastEM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['resolve_search_only'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)
        nllh_nni = -max_score
        
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_4 failed.")
    """

    # resolve polytomies set to true on starting tree without polytomies
    def test_5(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 0, 0, 0], 'b':[1, 1, 1, 0, 0], 'c':[0, 0, 0, 1, 0], 'd':[0, 0, 0, 1, 0]}
        nllh_bf = 11.809140958825493 # precomputed from brute-force
        
        T0 = '((a,c),(b,d));'
        parser_tree_out = parse_tree(read_tree_newick(T0))
        out_phylogeny = parse_data(parser_tree_out, {'charMtrx':msa}, {'Q': Q}) # obj used to init the EMoptimizer
        Q_recode = out_phylogeny.mutation_priors
        ordered_leaf_labels = parser_tree_out['relabeling_vector'][:parser_tree_out['num_leaves']]
        charMtrx_recode = out_phylogeny.character_matrix

        data = {'charMtrx':msa, 'ordered_leaf_labels': [ordered_leaf_labels], 'charMtrx_recode': [charMtrx_recode]}
        prior = {'Q':Q, 'Q_recode': [Q_recode]}
        params = {'phi':0,'nu':0}
        
        # topology search with EM_solver
        myTopoSearch = Topology_search([T0],fastEM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['resolve_search_only'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=10)
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_5 failed.")
    
    # full topology search
    def test_6(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf,tree_bf = self.__brute_force_search__(msa,Q,['a','b','c','d'],solver=ML_solver)
        
        T0 = '((a,b),(c,d));'
        parser_tree_out = parse_tree(read_tree_newick(T0))
        out_phylogeny = parse_data(parser_tree_out, {'charMtrx':msa}, {'Q': Q}) # obj used to init the EMoptimizer
        Q_recode = out_phylogeny.mutation_priors
        ordered_leaf_labels = parser_tree_out['relabeling_vector'][:parser_tree_out['num_leaves']]
        charMtrx_recode = out_phylogeny.character_matrix

        data = {'charMtrx':msa, 'ordered_leaf_labels': [ordered_leaf_labels], 'charMtrx_recode': [charMtrx_recode]}
        prior = {'Q':Q, 'Q_recode': [Q_recode]}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search([T0],fastEM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        #my_strategy['resolve_search_only'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=10)
        nllh_nni = -max_score

        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_6 failed.") 
     
    """
    # resolve polytomies only [NOT IMPLEMENTED]
    def test_7(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf = 10.083048489886435 # precomputed from brute-force
 
        T0 = '((a,b),c,d);'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search([T0],fastEM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['resolve_search_only'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)

        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_7 failed.")
    
    # resolve polytomies only [NOT IMPLEMENTED]
    def test_8(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf = 10.08306053772235
        
        T0 = '((c,d),a,b);'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search([T0],fastEM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['resolve_search_only'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=1)
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_8 failed.")
    """ 
    
    # enforce ultrametric
    def test_9(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf,bf_tree = self.__brute_force_search__(msa,Q,['a','b','c','d'],solver=ML_solver,ultra_constr=True)
        T0 = '((a,b),(c,d));'
        parser_tree_out = parse_tree(read_tree_newick(T0))
        out_phylogeny = parse_data(parser_tree_out, {'charMtrx':msa}, {'Q': Q}) # obj used to init the EMoptimizer
        Q_recode = out_phylogeny.mutation_priors
        ordered_leaf_labels = parser_tree_out['relabeling_vector'][:parser_tree_out['num_leaves']]
        charMtrx_recode = out_phylogeny.character_matrix

        data = {'charMtrx':msa, 'ordered_leaf_labels': [ordered_leaf_labels], 'charMtrx_recode': [charMtrx_recode]}
        prior = {'Q':Q, 'Q_recode': [Q_recode]}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search([T0],fastEM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['ultra_constr'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=10)
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_9 failed.")
    
    # enforce ultrametric and starting tree is a star-tree
    def test_10(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf,bf_tree = self.__brute_force_search__(msa,Q,['a','b','c','d'],solver=ML_solver,ultra_constr=True)
        #nllh_bf = 7.0063474330891236 # pre-computed using brute-force search
        T0 = '(a,b,c,d);'
        parser_tree_out = parse_tree(read_tree_newick(T0))
        out_phylogeny = parse_data(parser_tree_out, {'charMtrx':msa}, {'Q': Q}) # obj used to init the EMoptimizer
        Q_recode = out_phylogeny.mutation_priors
        ordered_leaf_labels = parser_tree_out['relabeling_vector'][:parser_tree_out['num_leaves']]
        charMtrx_recode = out_phylogeny.character_matrix

        data = {'charMtrx':msa, 'ordered_leaf_labels': [ordered_leaf_labels], 'charMtrx_recode': [charMtrx_recode]}
        prior = {'Q':Q, 'Q_recode': [Q_recode]}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search([T0],fastEM_solver,data=data,prior=prior,params=params)
        my_strategy = deepcopy(DEFAULT_STRATEGY)
        my_strategy['ultra_constr'] = True
        best_tree,max_score,best_params = myTopoSearch.search(maxiter=200,verbose=False,strategy=my_strategy,nreps=10)

        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=4,msg="TopoSearchTest: test_10 failed.")

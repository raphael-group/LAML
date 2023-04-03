import os 
import unittest
from problin_libs.EM_solver import EM_solver
from problin_libs.ML_solver import ML_solver
from problin_libs.Topology_search import Topology_search
from treeswift import *

class TopoSearchTest(unittest.TestCase):
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

    def __brute_force_search__(self,msa,Q,L,ultra_constr=False,initials=1):
        topos = self.__list_topologies__(L)
        best_nllh = float("inf")
        best_tree = ""
        for T in topos:
            #mySolver = EM_solver(msa,Q,T)
            mySolver = EM_solver(T,{'charMtrx':msa},{'Q':Q})
            nllh = mySolver.optimize(initials=initials,verbose=-1,ultra_constr=ultra_constr)
            if nllh < best_nllh:
                best_nllh = nllh
                best_tree = mySolver.tree.newick()
        return best_nllh,best_tree 

    def test_1(self):
        #L = ['a','b','c','d']
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 1, 1, 1], 'b':[0, 0, 0, 0, 0], 'c':[0, 0, 0, 0, 0], 'd':[1, 1, 1, 1, 1]}
        #nllh_bf,_ = self.__brute_force_search__(msa,Q,L)
        nllh_bf = 4.581468106634933 # pre-computed using brute-force search
        
        T0 = '((a,b)e,(c,d)f)g;'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        mySolver = EM_solver(T0,data,prior)
        mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
        T0_brlen = mySolver.tree.newick()
        params = {'nu':mySolver.params.nu,'phi':mySolver.params.phi}

        # topology search with EM_solver
        myTopoSearch_EM = Topology_search(T0_brlen,EM_solver,data=data,prior=prior,params=params)
        nni_replicates = myTopoSearch_EM.search(maxiter=1000, verbose=False)
        T1,_ = nni_replicates[0][1][-1]
        mySolver = EM_solver(T1,data,prior,params)
        nllh_nni_EM = mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
        
        self.assertAlmostEqual(nllh_bf,nllh_nni_EM,places=5,msg="TopoSearchTest: test_1 failed.")
        
    def test_2(self):
        #L = ['a','b','c','d']
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 1, 1, 1], 'b':[0, 0, 0, 0, 0], 'c':[0, 0, 0, 0, 0], 'd':[1, 1, 1, 1, 1]}
        #nllh_bf,_ = self.__brute_force_search__(msa,Q,L)
        nllh_bf = 4.581468106634933 # pre-computed using brute-force search
        
        T0 = '((a,b)e,(c,d)f)g;'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        mySolver = EM_solver(T0,data,prior)
        mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
        T0_brlen = mySolver.tree.newick()
        params = {'nu':mySolver.params.nu,'phi':mySolver.params.phi}
        
        # topology search with ML_solver
        myTopoSearch_ML = Topology_search(T0_brlen,ML_solver,data=data,prior=prior,params=params)
        nni_replicates = myTopoSearch_ML.search(maxiter=1000, verbose=False)
        T2,_ = nni_replicates[0][1][-1]
        mySolver = EM_solver(T2,data,prior,params)
        nllh_nni_ML = mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
        
        self.assertAlmostEqual(nllh_bf,nllh_nni_ML,places=5,msg="TopoSearchTest: test_1 failed.")

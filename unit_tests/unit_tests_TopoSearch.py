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
            mySolver = EM_solver(T,{'charMtrx':msa},{'Q':Q})
            nllh = mySolver.optimize(initials=initials,verbose=-1,ultra_constr=ultra_constr)
            #print(T,nllh)
            if nllh < best_nllh:
                best_nllh = nllh
                best_tree = mySolver.tree.newick()
        return best_nllh,best_tree 
    
    # topology search with EM_solver
    def test_1(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 1, 1, 1], 'b':[0, 0, 0, 0, 0], 'c':[0, 0, 0, 0, 0], 'd':[1, 1, 1, 1, 1]}
        nllh_bf = 4.581468106634933 # pre-computed using brute-force search
        
        T0 = '((a,b),(c,d));'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        mySolver = EM_solver(T0,data,prior)
        mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
        T0_brlen = mySolver.get_tree_newick()
        params = {'nu':mySolver.params.nu,'phi':mySolver.params.phi}

        # topology search with EM_solver
        myTopoSearch_EM = Topology_search(T0_brlen,EM_solver,data=data,prior=prior,params=params)
        nni_replicates = myTopoSearch_EM.search(maxiter=200,nreps=5,verbose=False)

        max_score = -float("inf")
        T1 = ""
        for score,tree_topos in nni_replicates:
            if score > max_score:
                max_score = score
                T1,_ = tree_topos[-1]
        mySolver = EM_solver(T1,data,prior,params)
        nllh_nni_EM = mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
        
        self.assertAlmostEqual(nllh_bf,nllh_nni_EM,places=5,msg="TopoSearchTest: test_5 failed.")
    
    # topology search with ML_solver
    def test_2(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 1, 1, 1], 'b':[0, 0, 0, 0, 0], 'c':[0, 0, 0, 0, 0], 'd':[1, 1, 1, 1, 1]}
        nllh_bf = 4.581468106634933 # pre-computed using brute-force search
        
        T0 = '((a,b),(c,d));'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        mySolver = EM_solver(T0,data,prior)
        mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
        T0_brlen = mySolver.get_tree_newick()
        params = {'nu':mySolver.params.nu,'phi':mySolver.params.phi}
        
        # topology search with ML_solver
        myTopoSearch_ML = Topology_search(T0_brlen,ML_solver,data=data,prior=prior,params=params)
        nni_replicates = myTopoSearch_ML.search(maxiter=200,verbose=False,nreps=5)

        max_score = -float("inf")
        T1 = ""
        for score,tree_topos in nni_replicates:
            if score > max_score:
                max_score = score
                T1,_ = tree_topos[-1]

        mySolver = EM_solver(T1,data,prior,params)
        nllh_nni_ML = mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
        
        self.assertAlmostEqual(nllh_bf,nllh_nni_ML,places=5,msg="TopoSearchTest: test_2 failed.")
    
    # resolve polytomies and continue beyond that
    def test_3(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 1, 1, 1], 'b':[0, 0, 0, 0, 0], 'c':[0, 0, 0, 0, 0], 'd':[1, 1, 1, 1, 1]}
        nllh_bf = 4.581468106634933 # pre-computed using brute-force search
        
        T0 = '(a:0.5,b:0.5,c:0.5,d:0.5):0.4;'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'nu':0,'phi':0}
        
        # topology search with ML_solver
        myTopoSearch = Topology_search(T0,EM_solver,data=data,prior=prior,params=params)
        nni_replicates = myTopoSearch.search(maxiter=200,verbose=False,nreps=5)
        max_score = -float("inf")
        T1 = ""
        for score,tree_topos in nni_replicates:
            if score > max_score:
                max_score = score
                T1,_ = tree_topos[-1]
        mySolver = EM_solver(T1,data,prior,params)
        nllh_nni = mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
        
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=5,msg="TopoSearchTest: test_3 failed.")
     
    # only resolve polytomies
    def test_4(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 1, 1, 1], 'b':[0, 0, 0, 0, 0], 'c':[0, 0, 0, 0, 0], 'd':[1, 1, 1, 1, 1]}
        nllh_bf = 9.720780114038451 # # pre-computed using brute-force search
        
        T0 = '((a:0.5,b:0.5):0.1,c:0.5,d:0.5):0.4;'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'nu':0,'phi':0}
        
        myTopoSearch = Topology_search(T0,EM_solver,data=data,prior=prior,params=params)
        nni_replicates = myTopoSearch.search(maxiter=200,verbose=False,strategy={'resolve_polytomies':True,'only_marked':True,'optimize':False,'ultra_constr':False},nreps=5)
        max_score = -float("inf")
        T1 = ""
        for score,tree_topos in nni_replicates:
            if score > max_score:
                max_score = score
                T1,_ = tree_topos[-1]
        mySolver = EM_solver(T1,data,prior,params)
        nllh_nni = mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
        
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=5,msg="TopoSearchTest: test_4 failed.")
    
    def test_5(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 1, 1, 1], 'b':[0, 0, 0, 0, 0], 'c':[0, 0, 0, 0, 0], 'd':[1, 1, 1, 1, 1]}
        nllh_bf = 4.581468106634933 # pre-computed using brute-force search
        
        T0 = '((a,b),(c,d));'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        # topology search with EM_solver
        myTopoSearch = Topology_search(T0,EM_solver,data=data,prior=prior,params=params)
        nni_replicates = myTopoSearch.search(maxiter=200,verbose=False,strategy={'resolve_polytomies':True,'only_marked':False,'optimize':True,'ultra_constr':False},nreps=1)

        max_score = -float("inf")
        T1 = ""
        for score,tree_topos in nni_replicates:
            if score > max_score:
                max_score = score
                T1,_ = tree_topos[-1]
        mySolver = EM_solver(T1,data,prior,params)
        nllh_nni_EM = mySolver.optimize(initials=1,verbose=-1,ultra_constr=False)
        
        self.assertAlmostEqual(nllh_bf,nllh_nni_EM,places=5,msg="TopoSearchTest: test_5 failed.")
    
    def test_6(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf = 7.552170051824473 # pre-computed using brute-force search
        
        T0 = '((a,b),(c,d));'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search(T0,EM_solver,data=data,prior=prior,params=params)
        nni_replicates = myTopoSearch.search(maxiter=200,verbose=False,strategy={'resolve_polytomies':True,'only_marked':False,'optimize':True,'ultra_constr':False},nreps=1)

        max_score = -float("inf")
        T1 = ""
        for score,tree_topos in nni_replicates:
            if score > max_score:
                max_score = score
                T1,_ = tree_topos[-1]
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=5,msg="TopoSearchTest: test_6 failed.")
        
    def test_7(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf = 7.552170051824473 # pre-computed using brute-force search
        
        T0 = '((a,b),c,d);'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search(T0,EM_solver,data=data,prior=prior,params=params)
        nni_replicates = myTopoSearch.search(maxiter=200,verbose=False,strategy={'resolve_polytomies':True,'only_marked':False,'optimize':True,'ultra_constr':False},nreps=1)

        max_score = -float("inf")
        T1 = ""
        for score,tree_topos in nni_replicates:
            if score > max_score:
                max_score = score
                T1,_ = tree_topos[-1]
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=5,msg="TopoSearchTest: test_7 failed.")
    
    def test_8(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf = 13.609086537509608 # pre-computed using brute-force search
        #self.__brute_force_search__(msa,Q,['a','b','c','d'],ultra_constr=False,initials=1)
        
        T0 = '((a,b),c,d);'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search(T0,EM_solver,data=data,prior=prior,params=params)
        nni_replicates = myTopoSearch.search(maxiter=200,verbose=False,strategy={'resolve_polytomies':True,'only_marked':True,'optimize':True,'ultra_constr':False},nreps=1)

        max_score = -float("inf")
        T1 = ""
        for score,tree_topos in nni_replicates:
            if score > max_score:
                max_score = score
                T1,_ = tree_topos[-1]
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=5,msg="TopoSearchTest: test_8 failed.")
    
    def test_9(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf = 8.54902142585155 # pre-computed using brute-force search
        #self.__brute_force_search__(msa,Q,['a','b','c','d'],ultra_constr=True,initials=1)
        
        #T0 = '((a,b),c,d);'
        T0 = '((a,b),(c,d));'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search(T0,EM_solver,data=data,prior=prior,params=params)
        nni_replicates = myTopoSearch.search(maxiter=200,verbose=False,strategy={'resolve_polytomies':True,'only_marked':False,'optimize':True,'ultra_constr':True},nreps=1)

        max_score = -float("inf")
        T1 = ""
        for score,tree_topos in nni_replicates:
            if score > max_score:
                max_score = score
                T1,_ = tree_topos[-1]
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=5,msg="TopoSearchTest: test_9 failed.")
    
    def test_10(self):
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[0, 1, 1, 1, 1], 'b':[1, 0, 0, 0, 0], 'c':[1, 0, 0, 0, 0], 'd':[0, 1, 1, 1, 1]}
        nllh_bf = 8.54902142585155 # pre-computed using brute-force search
        
        T0 = '((a,b),c,d);'
        data = {'charMtrx':msa}
        prior = {'Q':Q}
        params = {'phi':0,'nu':0}
        
        myTopoSearch = Topology_search(T0,EM_solver,data=data,prior=prior,params=params)
        nni_replicates = myTopoSearch.search(maxiter=200,verbose=False,strategy={'resolve_polytomies':True,'only_marked':False,'optimize':True,'ultra_constr':True},nreps=1)

        max_score = -float("inf")
        T1 = ""
        for score,tree_topos in nni_replicates:
            if score > max_score:
                max_score = score
                T1,_ = tree_topos[-1]
        nllh_nni = -max_score
        self.assertAlmostEqual(nllh_bf,nllh_nni,places=5,msg="TopoSearchTest: test_10 failed.")

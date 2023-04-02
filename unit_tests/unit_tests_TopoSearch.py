import os 
import unittest
from problin_libs.EM_solver import EM_solver
from problin_libs.Topology_search import Topology_search

class TopoSearchTest(unittest.TestCase):
    def test_1(self):
        T = "[&R] ((a:0.5,b:1)e:2,(c:1,d:0.5)f:1)g:1;"
        Q = [{0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}, {0:0, 1:1.0}]
        msa = {'a':[1, 1, 1, 1, 1], 'b':[0, 0, 0, 0, 0], 'c':[0, 0, 0, 0, 0], 'd':[1, 1, 1, 1, 1]}
        correct_topology = "[&R] (((a:0.5,d:0.5)f:1,c:1)e:2,b:1)g:1;"
        
        mySolver = EM_solver(msa,Q,T)
        myTopoSearch = Topology_search(mySolver)
        myTopoSearch.search(maxiter=10, verbose=False, strategy="vanilla", trynextbranch=False, conv=0.2)

        opt_topo = mySolver.tree.newick()
        self.assertEqual(opt_topo, correct_topology, msg="TopoSearchTest: test_1 failed.")

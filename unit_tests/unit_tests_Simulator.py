from itertools import combinations
import treeswift
import unittest
from problin_libs.sim_lib import *
from unit_tests.utils import count_all, count_missing, setup, calc_expected

class SimulatorTest(unittest.TestCase):

    def test_1(self):
        randomseed = 1984
        nreps = 20
        m, mu, k, d = 30, 1, 30, 0.0
        
        tree = treeswift.read_tree_newick("unit_tests/test_data/test_simulator/test.tre") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                Q = sim_Q(k, m)
                cmtx, char_full = simulate_seqs(tree, Q, mu, dropout_rate=d)
                allreps[rep] = cmtx
        nzeros, nmissing, total = count_all(allreps)
        propzeros = float(nzeros)/(total - nmissing)
        self.assertAlmostEqual(d,propmissing,places=2,msg="SimulatorTest: test_1 failed, " + str(d) + str(propmissing))

    def test_2(self):
        randomseed = 1984
        nreps = 50
        m, mu, k, d = 30, 1, 30, 0.01
        tree = treeswift.read_tree_newick("unit_tests/test_data/test_simulator/test.tre") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                Q = sim_Q(k, m)
                cmtx, char_full = simulate_seqs(tree, Q, mu, dropout_rate=d)
                allreps[rep] = cmtx

        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(d, propmissing,places=2,msg="SimulatorTest: test_2 failed, " + str(d) + " " + str(propmissing))

    def test_3(self):
        randomseed = 1984
        nreps = 50
        m, mu, k, d = 30, 1, 30, 0.3
        tree = treeswift.read_tree_newick("unit_tests/test_data/test_simulator/test.tre") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                Q = sim_Q(k, m)
                cmtx, char_full = simulate_seqs(tree, Q, mu, dropout_rate=d)
                allreps[rep] = cmtx

        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(d,propmissing,places=2,msg="SimulatorTest: test_3 failed, " + str(d) + " " + str(propmissing))
    
    def test_4(self):
        randomseed = 1984
        nreps = 50
        m, mu, k, d = 30, 1, 200, 0.01
        tree = treeswift.read_tree_newick("unit_tests/test_data/test_simulator/test.tre") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                Q = sim_Q(k, m)
                cmtx, char_full = simulate_seqs(tree, Q, mu, dropout_rate=d)
                allreps[rep] = cmtx

        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(d, propmissing,places=2,msg="SimulatorTest: test_4 failed, " + str(d) + " " + str(propmissing))

    def test_5(self):
        randomseed = 1984
        nreps = 50
        m, mu, k, s = 30, 1, 30, 0.0
        tree = treeswift.read_tree_newick("unit_tests/test_data/test_simulator/test.tre") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                Q = sim_Q(k, m)
                random.seed(randomseed+rep)
                cmtx, char_full = simulate_seqs(tree, Q, mu, with_heritable=True, silencing_rate=s)
                allreps[rep] = cmtx
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(propmissing, s, places=2, msg="SimulatorTest: test_5 failed, " + str(s) + " " + str(propmissing))

    def test_6A(self):
        randomseed = 1984
        nreps = 75 
        for s in [0.0, 0.06, 0.13, 0.20]:
            m, mu, k = 30, 1, 200 
            height, bl = 6, 0.2
            tree = get_balanced_tree(height, bl) # set balanced tree 
            tree = treeswift.read_tree_newick(tree)
            d = height * bl # tree.height(weighted=True)

            allreps = dict()
            for rep in range(nreps):
                    Q = sim_Q(k, m)
                    random.seed(randomseed+rep)
                    cmtx, char_full = simulate_seqs(tree, Q, mu, silencing_rate=s)
                    allreps[rep] = cmtx
            nzeros, nmissing, total = count_all(allreps)
            propmissing = float(nmissing)/total
            exp_propmissing = 1 - exp(-d*s) 
            self.assertAlmostEqual(propmissing, exp_propmissing, places=2, msg="SimulatorTest: test_6A failed, silencing rate:" + str(s) + " propmissing: " + str(propmissing) + " exp_propmissing:" + str(exp_propmissing))
    
    def test_6B(self):
        randomseed = 1984
        nreps = 75 
        for s in [0.0, 0.06, 0.13, 0.20]:
            m, mu, k = 30, 1, 200
            height, bl = 6, 0.2
            tree = get_balanced_tree(height, bl) # set balanced tree 
            tree = treeswift.read_tree_newick(tree)
            d = height * bl # tree.height(weighted=True)

            allreps = dict()
            for rep in range(nreps):
                    Q = sim_Q(k, m)
                    random.seed(randomseed+rep)
                    cmtx, char_full = simulate_seqs(tree, Q, mu, silencing_rate=s)
                    allreps[rep] = cmtx
            nzeros, nmissing, total = count_all(allreps)
            propzeros = float(nzeros)/(total - nmissing)
            exp_propzeros = exp(-d*(1+s)) 
            self.assertAlmostEqual(propzeros, exp_propzeros, places=2, msg="SimulatorTest: test_6B failed, propzeros: " + str(propzeros) + " exp_propzeros:" + str(exp_propzeros))

    
    def test_7(self):
        # test that different replicates give different cmtxs
        nreps = 10 
        m, mu, k = 30, 1, 30
        tree = treeswift.read_tree_newick("unit_tests/test_data/test_simulator/test.tre") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
            Q = sim_Q(k, m)
            allreps[rep], _ = simulate_seqs(tree, Q, mu)
        for x, y in list(combinations(range(nreps), 2)):
            self.assertFalse(allreps[x] == allreps[y], msg="SimulatorTest: test_7 failed, produces equal character matrices over different replicates.")

    def test_8(self):
        # write chi squared test on tree with 
        nreps = 10
        m, mu, k, s = 30, 1, 200, 0.05
        tree = treeswift.read_tree_newick("unit_tests/test_data/test_simulator/test.tre") #set unbalanced cass tree 
        allreps = dict()
        Q = sim_Q(k, m)
        for rep in range(nreps):
            allreps[rep], _ = simulate_seqs(tree, Q, mu, with_heritable=True, silencing_rate=s)
        for c in list(range(1,m)) + ["?"]:
            for node_label in allreps[rep]:
                self.assertTrue(calc_expected(node_label, 3.98, k, c, Q, allreps, s), msg="SimulatorTest: test_8 failed, unexpected character distribution for node_label: " + str(node_label) + " and char: " + str(c) ) 

    # TODO Add tests for the full character matrix
    # TODO Test both types of missing data


import treeswift
import unittest
from problin_libs.sim_lib import *

class SimulatorTest(unittest.TestCase):

    def test_1(self):
        randomseed = 1984
        nreps = 20
        m, mu, k, d = 30, 0.1, 30, 0.0
        
        tree = treeswift.read_tree_newick("/n/fs/ragr-research/projects/problin/unit_tests/cass_n150m30_nomissing.nwk") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                cmtx = setup(tree, m, mu, k)
                allreps[rep] = simulate_dropout(cmtx, d, randomseed + rep)
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(d,propmissing,places=2,msg="SimulatorTest: test_1 failed, " + str(d) + str(propmissing))

    def test_2(self):
        randomseed = 1984
        nreps = 100
        m, mu, k, d = 30, 0.1, 30, 0.01
        tree = treeswift.read_tree_newick("/n/fs/ragr-research/projects/problin/unit_tests/cass_n150m30_nomissing.nwk") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                cmtx = setup(tree, m, mu, k)
                allreps[rep] = simulate_dropout(cmtx, d, randomseed + rep)

        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(d, propmissing,places=2,msg="SimulatorTest: test_2 failed, " + str(d) + " " + str(propmissing))

    def test_3(self):
        randomseed = 1984
        nreps = 100 
        m, mu, k, d = 30, 0.1, 30, 0.1
        tree = treeswift.read_tree_newick("/n/fs/ragr-research/projects/problin/unit_tests/cass_n150m30_nomissing.nwk") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                cmtx = setup(tree, m, mu, k)
                allreps[rep] = simulate_dropout(cmtx, d, randomseed + rep)
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(d, propmissing,places=1,msg="SimulatorTest: test_3 failed, " + str(d) + " " + str(propmissing))
    
    def test_4(self):
        randomseed = 1984
        nreps = 100 
        m, mu, k, d = 30, 0.1, 30, 0.3
        tree = treeswift.read_tree_newick("/n/fs/ragr-research/projects/problin/unit_tests/cass_n150m30_nomissing.nwk") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                cmtx = setup(tree, m, mu, k)
                allreps[rep] = simulate_dropout(cmtx, d, randomseed + rep)
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(d,propmissing,places=1,msg="SimulatorTest: test_4 failed, " + str(d) + " " + str(propmissing))
    
    def test_5(self):
        randomseed = 1984
        nreps = 100
        m, mu, k, d = 30, 0.1, 200, 0.01
        tree = treeswift.read_tree_newick("/n/fs/ragr-research/projects/problin/unit_tests/cass_n150m30_nomissing.nwk") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                cmtx = setup(tree, m, mu, k)
                allreps[rep] = simulate_dropout(cmtx, d, randomseed + rep)

        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(d, propmissing,places=2,msg="SimulatorTest: test_5 failed, " + str(d) + " " + str(propmissing))

    def test_6(self):
        randomseed = 1984
        nreps = 100 
        m, mu, k, d = 30, 0.1, 200, 0.1
        tree = treeswift.read_tree_newick("/n/fs/ragr-research/projects/problin/unit_tests/cass_n150m30_nomissing.nwk") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                cmtx = setup(tree, m, mu, k)
                allreps[rep] = simulate_dropout(cmtx, d, randomseed + rep)
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(d, propmissing,places=1,msg="SimulatorTest: test_6 failed, " + str(d) + " " + str(propmissing))

    def test_7(self):
        randomseed = 1984
        nreps = 100 
        m, mu, k, s = 30, 0.1, 200, 0.0
        tree = treeswift.read_tree_newick("/n/fs/ragr-research/projects/problin/unit_tests/cass_n150m30_nomissing.nwk") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                Q = sim_Q(k, m)
                allreps[rep] = simulate_seqs(tree, Q, mu, with_heritable=True, silencing_rate=s, s=randomseed+rep)
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(propmissing, s, places=2, msg="SimulatorTest: test_7 failed, " + str(s) + " " + str(propmissing))

    def test_8(self):
        randomseed = 1984
        nreps = 100 
        m, mu, k, s = 30, 0.1, 200, 0.1e-3
        tree = treeswift.read_tree_newick("/n/fs/ragr-research/projects/problin/unit_tests/cass_n150m30_nomissing.nwk") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                Q = sim_Q(k, m)
                allreps[rep] = simulate_seqs(tree, Q, mu, with_heritable=True, silencing_rate=s, s=randomseed+rep)
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertLess(propmissing, s,msg="SimulatorTest: test_8 failed, " + str(s) + " " + str(propmissing))
    
    def test_9(self):
        randomseed = 1984
        nreps = 100 
        m, mu, k, s = 30, 0.1, 200, 1e-4
        tree = get_balanced_tree(8, 0.2) # set balanced tree 
        tree = treeswift.read_tree_newick(tree)

        allreps = dict()
        for rep in range(nreps):
                randomseed = randomseed + rep
                Q = sim_Q(k, m)
                allreps[rep] = simulate_seqs(tree, Q, mu, with_heritable=True, silencing_rate=s, s=randomseed)
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        propzeros = float(nzeros)/total
        # print("nmissing", nmissing, "propmissing", propmissing)
        self.assertGreater(nmissing,0,msg="SimulatorTest: test_9A failed. No missing data when silencing rate was " + str(s)) 
        self.assertLess(propmissing, s, msg="SimulatorTest: test_9B failed.")

    def test_10(self):
        randomseed = 1984
        nreps = 100 
        m, mu, k = 30, 0.1, 200 
        height, bl = 8, 0.2
        tree = get_balanced_tree(height, bl) # set balanced tree 
        tree = treeswift.read_tree_newick(tree)

        allreps = dict()
        for rep in range(nreps):
                randomseed = randomseed + rep
                Q = sim_Q(k, m)
                allreps[rep] = simulate_seqs(tree, Q, mu, s=randomseed)
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        propzeros = float(nzeros)/total
        truezeros = 1 - exp(-height * bl) # 1 - probability of mutating from root to leaf
        self.assertAlmostEqual(truezeros,propzeros,places=1,msg="SimulatorTest: test_10A failed.")
        
        '''
        tree = treeswift.read_tree_newick("/n/fs/ragr-research/projects/problin/unit_tests/cass_n150m30_nomissing.nwk") #set unbalanced cass tree 
        allreps = dict()
        for rep in range(nreps):
                randomseed = randomseed + rep
                Q = sim_Q(k, m)
                allreps[rep] = simulate_seqs(tree, Q, mu, s=randomseed)
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        propzeros = float(nzeros)/total
        truezeros = 1 - exp(-tree.height()) # 1 - probability of mutating from root to leaf
        self.assertAlmostEqual(truezeros,propzeros,places=1,msg="SimulatorTest: test_10B failed.")
        '''
    
    def test_11(self):
        randomseed = 1984
        nreps = 100 
        m, mu, k, s = 30, 0.1, 200, 0.1e-3
        tree = treeswift.read_tree_newick("/n/fs/ragr-research/projects/problin/unit_tests/cass_n150m30_nomissing.nwk") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                Q = sim_Q(k, m)
                allreps[rep] = simulate_seqs(tree, Q, mu, with_heritable=True, silencing_rate=s, s=randomseed+rep)
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertGreater(nmissing,0,msg="SimulatorTest: test_11A failed. No missing data when silencing rate was " + str(s)) 
        self.assertLess(propmissing, s,msg="SimulatorTest: test_11B failed, " + str(s) + " " + str(propmissing))
    
    def test_12(self):
        randomseed = 1984
        nreps = 100 
        m, mu, k, s = 30, 0.1, 200, 0.1e-2
        tree = treeswift.read_tree_newick("/n/fs/ragr-research/projects/problin/unit_tests/cass_n150m30_nomissing.nwk") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                Q = sim_Q(k, m)
                allreps[rep] = simulate_seqs(tree, Q, mu, with_heritable=True, silencing_rate=s, s=randomseed+rep)
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertGreater(nmissing,0,msg="SimulatorTest: test_12A failed. No missing data when silencing rate was " + str(s)) 
        self.assertLess(propmissing, s,msg="SimulatorTest: test_12B failed, " + str(s) + " " + str(propmissing))
    
    def test_13(self):
        randomseed = 1984
        nreps = 100 
        m, mu, k, s = 30, 0.1, 200, 0.1e-1
        tree = treeswift.read_tree_newick("/n/fs/ragr-research/projects/problin/unit_tests/cass_n150m30_nomissing.nwk") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                Q = sim_Q(k, m)
                allreps[rep] = simulate_seqs(tree, Q, mu, with_heritable=True, silencing_rate=s, s=randomseed+rep)
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertGreater(nmissing,0,msg="SimulatorTest: test_13A failed. No missing data when silencing rate was " + str(s)) 
        self.assertLess(propmissing, s, msg="SimulatorTest: test_13B failed, " + str(s) + " " + str(propmissing))


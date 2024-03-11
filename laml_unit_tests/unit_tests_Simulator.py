from itertools import combinations
import treeswift
import unittest
from scmail_libs.sim_lib import *
from unit_tests.utils import count_all, count_missing, setup, calc_expected

class SimulatorTest(unittest.TestCase):

    # TODO: Change tests to use significance/chi tests instead of assertAlmostEqual

    def test_1(self):
        randomseed = 1984
        nreps = 20
        m, mu, k, d = 30, 1, 500, 0.0
        
        tree = treeswift.read_tree_newick("unit_tests/test_data/test_simulator/test.tre") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                Q = sim_Q(k, m)
                cmtx, char_full = simulate_seqs(tree, Q, mu, dropout_rate=d)
                allreps[rep] = cmtx
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(d,propmissing,places=2,msg="SimulatorTest: test_1 failed, " + str(d) + str(propmissing))

    def test_2(self):
        randomseed = 1984
        nreps = 50
        m, mu, k, d = 30, 1, 500, 0.01
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
        m, mu, k, d = 30, 1, 500, 0.3
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
        m, mu, k, d = 30, 1, 500, 0.01
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
        m, mu, k, s = 30, 1, 500, 0.0
        tree = treeswift.read_tree_newick("unit_tests/test_data/test_simulator/test.tre") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
                Q = sim_Q(k, m)
                random.seed(randomseed+rep)
                cmtx, char_full = simulate_seqs(tree, Q, mu, silencing_rate=s)
                allreps[rep] = cmtx
        nzeros, nmissing, total = count_all(allreps)
        propmissing = float(nmissing)/total
        self.assertAlmostEqual(propmissing, s, places=2, msg="SimulatorTest: test_5 failed, " + str(s) + " " + str(propmissing))
    
    def test_6(self):
        #randomseed = 1984
        nreps = 100
        s = 0.0 #[0.0, 0.06, 0.13, 0.20]:
        m, mu, k = 30, 1, 500 
        height, bl = 6, 0.2
        tree = get_balanced_tree(height, bl) # set balanced tree 
        tree = treeswift.read_tree_newick(tree)
        d = (height + 1)* bl # tree.height(weighted=True)

        allreps = dict()
        leaves = [l.label for l in tree.traverse_leaves()]
        sampled_leaves = random.sample(leaves, 10)

        counter = dict() 
        for rep in range(nreps):
            Q = sim_Q(k, m)
            cmtx, char_full = simulate_seqs(tree, Q, mu, silencing_rate=s)
            allreps[rep] = cmtx
            for leaf in sampled_leaves:
                if leaf not in counter:
                    counter[leaf] = []
                counter[leaf].append(len(list(c for c in cmtx[leaf] if c == 0))/k)
        nzeros, nmissing, total = count_all(allreps)
        
        propmissing = float(nmissing)/total
        exp_propmissing = 1 - exp(-d*s) 
        self.assertAlmostEqual(propmissing, exp_propmissing, places=2, msg="SimulatorTest: test_6A failed, silencing rate:" + str(s) + " propmissing: " + str(propmissing) + " exp_propmissing:" + str(exp_propmissing))

        for leaf in sampled_leaves:
            propzeros = sum(counter[leaf])/len(counter[leaf]) #float(nzeros)/(total - nmissing)
            exp_propzeros = exp(-d*(1+s)) 
            self.assertAlmostEqual(propzeros, exp_propzeros, places=2, msg="SimulatorTest: test_6B failed, propzeros: " + str(propzeros) + " exp_propzeros:" + str(exp_propzeros))
    
    def test_7(self):
        #randomseed = 1984
        nreps = 75 
        s = 0.6 #[0.0, 0.06, 0.13, 0.20]:
        m, mu, k = 30, 1, 500 
        height, bl = 6, 0.2
        tree = get_balanced_tree(height, bl) # set balanced tree 
        tree = treeswift.read_tree_newick(tree)
        d = (height + 1)* bl # tree.height(weighted=True)

        allreps = dict()
        leaves = [l.label for l in tree.traverse_leaves()]
        sampled_leaves = random.sample(leaves, 10)

        counter = dict() 
        for rep in range(nreps):
            Q = sim_Q(k, m)
            cmtx, char_full = simulate_seqs(tree, Q, mu, silencing_rate=s)
            allreps[rep] = cmtx
            for leaf in sampled_leaves:
                if leaf not in counter:
                    counter[leaf] = []
                counter[leaf].append(len(list(c for c in cmtx[leaf] if c == 0))/k)
        nzeros, nmissing, total = count_all(allreps)
        
        propmissing = float(nmissing)/total
        exp_propmissing = 1 - exp(-d*s) 
        self.assertAlmostEqual(propmissing, exp_propmissing, places=2, msg="SimulatorTest: test_7A failed, silencing rate:" + str(s) + " propmissing: " + str(propmissing) + " exp_propmissing:" + str(exp_propmissing))

        for leaf in sampled_leaves:
            propzeros = sum(counter[leaf])/len(counter[leaf]) #float(nzeros)/(total - nmissing)
            exp_propzeros = exp(-d*(1+s)) 
            self.assertAlmostEqual(propzeros, exp_propzeros, places=2, msg="SimulatorTest: test_7B failed, propzeros: " + str(propzeros) + " exp_propzeros:" + str(exp_propzeros))
    
    def test_8(self):
        #randomseed = 1984
        nreps = 75 
        s = 0.13 #[0.0, 0.06, 0.13, 0.20]:
        m, mu, k = 30, 1, 500 
        height, bl = 6, 0.2
        tree = get_balanced_tree(height, bl) # set balanced tree 
        tree = treeswift.read_tree_newick(tree)
        d = (height + 1)* bl # tree.height(weighted=True)

        allreps = dict()
        leaves = [l.label for l in tree.traverse_leaves()]
        sampled_leaves = random.sample(leaves, 10)

        counter = dict() 
        for rep in range(nreps):
            Q = sim_Q(k, m)
            cmtx, char_full = simulate_seqs(tree, Q, mu, silencing_rate=s)
            allreps[rep] = cmtx
            for leaf in sampled_leaves:
                if leaf not in counter:
                    counter[leaf] = []
                counter[leaf].append(len(list(c for c in cmtx[leaf] if c == 0))/k)
        nzeros, nmissing, total = count_all(allreps)
        
        propmissing = float(nmissing)/total
        exp_propmissing = 1 - exp(-d*s) 
        self.assertAlmostEqual(propmissing, exp_propmissing, places=2, msg="SimulatorTest: test_8A failed, silencing rate:" + str(s) + " propmissing: " + str(propmissing) + " exp_propmissing:" + str(exp_propmissing))

        for leaf in sampled_leaves:
            propzeros = sum(counter[leaf])/len(counter[leaf]) #float(nzeros)/(total - nmissing)
            exp_propzeros = exp(-d*(1+s)) 
            self.assertAlmostEqual(propzeros, exp_propzeros, places=2, msg="SimulatorTest: test_8B failed, propzeros: " + str(propzeros) + " exp_propzeros:" + str(exp_propzeros))
    
    def test_9(self):
        #randomseed = 1984
        nreps = 75 
        s = 0.20 #[0.0, 0.06, 0.13, 0.20]:
        m, mu, k = 30, 1, 500 
        height, bl = 6, 0.2
        tree = get_balanced_tree(height, bl) # set balanced tree 
        tree = treeswift.read_tree_newick(tree)
        d = (height + 1)* bl # tree.height(weighted=True)

        allreps = dict()
        leaves = [l.label for l in tree.traverse_leaves()]
        sampled_leaves = random.sample(leaves, 10)

        counter = dict() 
        for rep in range(nreps):
            Q = sim_Q(k, m)
            cmtx, char_full = simulate_seqs(tree, Q, mu, silencing_rate=s)
            allreps[rep] = cmtx
            for leaf in sampled_leaves:
                if leaf not in counter:
                    counter[leaf] = []
                counter[leaf].append(len(list(c for c in cmtx[leaf] if c == 0))/k)
        nzeros, nmissing, total = count_all(allreps)
        
        propmissing = float(nmissing)/total
        exp_propmissing = 1 - exp(-d*s) 
        self.assertAlmostEqual(propmissing, exp_propmissing, places=2, msg="SimulatorTest: test_9A failed, silencing rate:" + str(s) + " propmissing: " + str(propmissing) + " exp_propmissing:" + str(exp_propmissing))

        for leaf in sampled_leaves:
            propzeros = sum(counter[leaf])/len(counter[leaf]) #float(nzeros)/(total - nmissing)
            exp_propzeros = exp(-d*(1+s)) 
            self.assertAlmostEqual(propzeros, exp_propzeros, places=2, msg="SimulatorTest: test_9B failed, propzeros: " + str(propzeros) + " exp_propzeros:" + str(exp_propzeros))
    
    def test_10(self):
        # test that different replicates give different cmtxs
        nreps = 10 
        m, mu, k = 30, 1, 500 
        tree = treeswift.read_tree_newick("unit_tests/test_data/test_simulator/test.tre") #set unbalanced cass tree 

        allreps = dict()
        for rep in range(nreps):
            Q = sim_Q(k, m)
            allreps[rep], _ = simulate_seqs(tree, Q, mu)
        for x, y in list(combinations(range(nreps), 2)):
            self.assertFalse(allreps[x] == allreps[y], msg="SimulatorTest: test_10 failed, produces equal character matrices over different replicates.")

    def test_11(self):
        # write chi squared test on tree with 
        nreps = 10
        m, mu, k, s = 30, 1, 500, 0.05
        tree = treeswift.read_tree_newick("unit_tests/test_data/test_simulator/test.tre") #set unbalanced cass tree 
        allreps = dict()
        Q = sim_Q(k, m)
        for rep in range(nreps):
            allreps[rep], _ = simulate_seqs(tree, Q, mu, silencing_rate=s)
        for c in list(range(1,m)) + ["?"]:
            for node_label in allreps[rep]:
                self.assertTrue(calc_expected(node_label, 3.98, k, c, Q, allreps, s), msg="SimulatorTest: test_11 failed, unexpected character distribution for node_label: " + str(node_label) + " and char: " + str(c) ) 

    # Test both types of missing data and full data matrix
    def test_12(self):
        """
        # write chi squared test on tree with 
        nreps = 100
        m, mu, k, s, drop = 30, 1, 500, 0.13353, 0.14285
        height, bl = 6, 0.2
        tree = get_balanced_tree(height, bl) # set balanced tree 
        tree = treeswift.read_tree_newick(tree)
        d = (height + 1)* bl # tree.height(weighted=True)
        
        leaves = [l.label for l in tree.traverse_leaves()]
        sampled_leaves = random.sample(leaves, 10)

        counterzero, countermissing, allreps, allseqs = dict(), dict(), dict(), dict()
        Q = sim_Q(k, m)
        for rep in range(nreps):
            cmtx, all_cmtx = simulate_seqs(tree, Q, mu, silencing_rate=s, dropout_rate=drop)
            allreps[rep], allseqs[rep] = cmtx, all_cmtx
            for leaf in sampled_leaves:
                if leaf not in counterzero:
                    counterzero[leaf] = []
                    countermissing[leaf] = []
                #countermissing[leaf] = list(c for c in cmtx[leaf] if c == '?')
                counterzero[leaf].append(len(list(c for c in cmtx[leaf] if c == 0))/k) #/(k - sum(countermissing[leaf])))
        
        for leaf in sampled_leaves:
            # check zero proportions
            propzeros = sum(counterzero[leaf])/len(counterzero[leaf])
            
            exp_propzeros = exp(-d*(1+s)) #* (1 - drop)
            exp_propzeros += (1 - exp_propzeros) * drop 
            self.assertAlmostEqual(propzeros, exp_propzeros, places=2, msg="SimulatorTest: test_12 failed, propzeros: " + str(propzeros) + " exp_propzeros:" + str(exp_propzeros))
        """

    def test_13(self):
        pass 
"""
        # write chi squared test on tree with 
        nreps = 100
        m, mu, k, s, drop = 30, 1, 500, 0.13353, 0.14285
        height, bl = 6, 0.2
        tree = get_balanced_tree(height, bl) # set balanced tree 
        tree = treeswift.read_tree_newick(tree)
        d = (height + 1)* bl # tree.height(weighted=True)
        
        leaves = [l.label for l in tree.traverse_leaves()]
        sampled_leaves = random.sample(leaves, 10)

        counterzero, countermissing, allreps, allseqs = dict(), dict(), dict(), dict()
        Q = sim_Q(k, m)
        for rep in range(nreps):
            cmtx, all_cmtx = simulate_seqs(tree, Q, mu, silencing_rate=s, dropout_rate=drop)
            allreps[rep], allseqs[rep] = cmtx, all_cmtx
            for leaf in sampled_leaves:
                if leaf not in counterzero:
                    counterzero[leaf] = []
                    countermissing[leaf] = []
                countermissing[leaf].append(len(list(c for c in cmtx[leaf] if c == '?'))/k)
                counterzero[leaf].append(len(list(c for c in cmtx[leaf] if c == 0))/(k - sum(countermissing[leaf])))
        nzeros, nmissing, total = count_all(allreps)
        
        for leaf in sampled_leaves:
            # check missing proportions
            propmissing = sum(countermissing[leaf])/len(countermissing[leaf]) 
            exp_propmissing = (1 - exp(-d*s)) + (exp(-d*s) * drop)
            self.assertAlmostEqual(propmissing, exp_propmissing, places=2, msg="SimulatorTest: test_13 failed" + str(s) + " propmissing: " + str(propmissing) + " exp_propmissing:" + str(exp_propmissing))

    def test_14(self):
        # write chi squared test on tree with 
        nreps = 100
        m, mu, k, s, drop = 30, 1, 500, 0.13353, 0.14285
        height, bl = 6, 0.2
        tree = get_balanced_tree(height, bl) # set balanced tree 
        tree = treeswift.read_tree_newick(tree)
        d = (height + 1)* bl # tree.height(weighted=True)
        
        allreps, allseqs = dict(), dict()
        Q = sim_Q(k, m)
        for rep in range(nreps):
            cmtx, all_cmtx = simulate_seqs(tree, Q, mu, silencing_rate=s, dropout_rate=drop)
            allreps[rep], allseqs[rep] = cmtx, all_cmtx
        
        nsilencing, total = 0, 0
        for rep in allseqs:
            mtx = allseqs[rep]
            for c in mtx:
                seq = mtx[c]
                nsilencing += sum([1 if (ch == 's') else 0 for ch in seq])
                total += len(seq)

        propsilencing = nsilencing/total
        exp_propsilencing = 1 - exp(-d*s) 

        self.assertAlmostEqual(propsilencing, exp_propsilencing, places=2, msg="SimulatorTest: test_14 failed: propsilencing: " + str(propsilencing) + " exp_propsilencing:" + str(exp_propsilencing))


    def test_15(self):
        nreps = 100
        m, mu, k, s, drop = 30, 1, 500, 0.134, 0.143
        height, bl = 6, 0.2
        tree = get_balanced_tree(height, bl) # set balanced tree 
        tree = treeswift.read_tree_newick(tree)
        
        allreps, allseqs = dict(), dict()
        Q = sim_Q(k, m)
        for rep in range(nreps):
            _, allseqs[rep] = simulate_seqs(tree, Q, mu, silencing_rate=s, dropout_rate=drop)
        
        ndropout, total = 0, 0
        for rep in allseqs:
            for c in allseqs[rep]:
                seq = allseqs[rep][c]
                ndropout += sum([1 if (ch == 'd') else 0 for ch in seq])
                total += len(seq)

        propdropout = ndropout/total
        exp_propdropout = 0.25 * 0.5

        self.assertAlmostEqual(propdropout, exp_propdropout, places=2, msg="SimulatorTest: test_15 failed: propdropout: " + str(propdropout) + " exp_propdropout:" + str(exp_propdropout))

"""

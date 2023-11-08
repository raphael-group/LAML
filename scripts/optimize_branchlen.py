import numpy as np
# from scripts.felsenstein import felsenstein 
import dendropy
import itertools
import scipy
from scipy import optimize
import random
import math
from Bio import SeqIO
import itertools
from math import log, exp
import subprocess
from numpy import linalg as LA
import pickle

def sets(seq_a, seq_b):
    # get the msa
    k = len(seq_a)
    
    ## calculate the sets  
    s_0, s_1a, s_1b, s_2, s_3 = set(), set(), set(), set(), set()
    for idx in range(len(seq_a)):
        c_a, c_b = seq_a[idx], seq_b[idx]
        if c_a == c_b:
            if c_a == 0:
                s_0.add(idx)
            else:
                s_2.add(idx)
        elif c_a == 0: # then c_b != 0
            s_1b.add(idx)
        elif c_b == 0: # then c_a != 0
            s_1a.add(idx)
        else:
            s_3.add(idx)
    
    assert len(s_0) + len(s_1a) + len(s_1b) + len(s_2) + len(s_3) == k
    return [s_0, s_1a, s_1b, s_2, s_3]

def get_taxon_name(node):
    return int(node.taxon.__str__().replace("'", ""))

def est_d(seq_a, seq_b):
    # print("est_d")
    # print(seq_a, seq_b)
    
    # leaf nodes a and b
    Z_a = seq_a.count(0)
    Z_b = seq_b.count(0)
    Z_ab = set()
    
    for i in range(len(seq_a)):
        a_i, b_i = seq_a[i], seq_b[i]
        if a_i == 0:
            Z_ab.add(i)
        if b_i == 0:
            Z_ab.add(i)
        if a_i != b_i and a_i != 0 and b_i != 0:
            Z_ab.add(i)
    Z_ab = len(Z_ab)
    
    # print("Z_a, Z_b, Z_ab", Z_a, Z_b, Z_ab)
    
    d_a = - np.log(Z_a/Z_ab)
    d_b = - np.log(Z_b/Z_ab)
    
    
    if d_a == -0.0:
        d_a = 0.0
    if d_b == -0.0:
        d_b = 0.0
    
    # print("d_a, d_b", d_a, d_b)
    return d_a, d_b

# Optimize the branch lengths w/ marginal likelihood, using scipy.optimize

def likelihood(x): 
    global k
    global s_0, s_1a, s_1b, s_2, s_3
    s_0_len, s_1a_len, s_1b_len, s_2_len, s_3_len = len(s_0), len(s_1a), len(s_1b), len(s_2), len(s_3)
    
    d_a, d_b, d_r = x
    q = 0.2
    
    p1 = -(s_1b_len + s_0_len) * d_a + (s_1a_len + s_3_len) * np.log(1 - math.exp(-d_a))
    p2 = -(s_1a_len + s_0_len) * d_b + (s_1b_len + s_3_len) * np.log(1 - math.exp(-d_b)) - (k - s_2_len) * d_r
    p3 = 0.0
    
    for i in range(s_2_len): # assuming that prob to any alpha is the same
        # q_ialpha is the transition rate at site i from 0 -> character at node a at site i
        # iterate over sites in s_2, get q for each alpha
        p3 += np.log(q**2 * (1 - math.exp(-d_a)) * (1 - math.exp(-d_b)) * math.exp(-d_r) + q*(1 - math.exp(-d_r)))
    
    return -(p1 + p2 + p3)


# def optimize_len(alphabet_size, k, a, b):
def optimize_len(alphabet_size, a, b, x0):
    x_star = []
    num_iter = 20    
    for i in range(num_iter):
        if i > 0:
            x0 = np.random.uniform(0, 5.0, 3)
        global s_0, s_1a, s_1b, s_2, s_3
        s_0, s_1a, s_1b, s_2, s_3 = sets(a, b)
        # eps = 1e-10
        dmax = -log(1/k)*2
        dmin = -log(1-1/k)/2
        bound = (dmin,dmax)

        # out = optimize.minimize(likelihood, x0, method="L-BFGS-B", options={'disp': False}, bounds=[bound, bound, bound])
        out = optimize.minimize(likelihood, x0, method="SLSQP", options={'disp': False}, bounds=[bound, bound, bound])

        x_star.append((out['fun'], out['x']))
        # w, v = LA.eig(out.hess_inv.todense())
        # print(w > 0.0)
        
    return x_star


def get_seq(seq_dict, leaf_node):
    return seq_dict[int(leaf_node.taxon.__str__().replace("'", ""))]
               
def get_idx(leaf_node):
    # print(leaf_node)
    return int(leaf_node.taxon.__str__().replace("'", ""))

def my_print_tree(t):
    for n in t.postorder_node_iter():
        n.label = ''
    print(t)

def run_trees(fname, topo, m, num_to_run=2):
    trees = []
    sets_of_four = []
    with open(fname) as r:
        four_leaves = []
        for record in SeqIO.parse(r, "fasta"):
            s = str(record.seq)
            if '_________' in record.seq:
                s = s.replace('_________', '')
                four_leaves.append([int(x) for x in s.split('|')])
                sets_of_four.append(four_leaves)
                four_leaves = []
            else:
                four_leaves.append([int(x) for x in s.split('|')])
    print(len(sets_of_four), "trees read.")

    all_dists = dict()
    all_seqs = dict()
    
    for idx, seqs in enumerate(sets_of_four):
        dists = dict()
        seq_dict = dict()
        num_seqs = len(seqs)

        if idx == num_to_run:
                break

        for i in range(num_seqs):
            seq_dict[i] = seqs[i]
            
        for pair in itertools.combinations(list(range(num_seqs)), 2):
            # if idx % 5 == 0:
                # print(idx)

            a, b = pair
            x0 = np.random.uniform(0, 5.0, 3)
            x_star = optimize_len(m, seqs[a], seqs[b], x0) 
            # print the stddev of x_star
            d_ab = sorted(x_star, key=lambda x: x[0], reverse=True)[0][1]
            d_a, d_b, d_r = d_ab
            dists[(a, b)] = [d_a, d_b, d_r]

        all_dists[idx] = dists
        all_seqs[idx] = seq_dict

        # estimate branch length for every pair of leaves
        for n in topo.postorder_node_iter():
            if n.is_internal(): # must be an internal node
                node_a, node_b = n.child_nodes()
                if node_a.is_leaf(): 
                    idx_a = get_idx(node_a)
                    seq_a = get_seq(seq_dict, node_a)
                if node_b.is_leaf():
                    idx_b = get_idx(node_b)
                    seq_b = get_seq(seq_dict, node_b)
                if node_a.is_leaf() and node_b.is_leaf():
                    if (idx_a, idx_b) in dists:
                        d_a, d_b, d_r = dists[(idx_a, idx_b)]
                        node_a.edge.length = d_a
                        node_b.edge.length = d_b
                        n.edge.length = d_r
                    elif (idx_b, idx_a) in dists:
                        d_b, d_a, d_r = dists[(idx_b, idx_a)]
                        node_a.edge.length = d_a
                        node_b.edge.length = d_b
                        n.edge.length = d_r
                    else:
                        print("dists", dists)
        trees.append(topo)
    return trees


def run_mats(fname, topo, m, num_to_run=2):
    trees = []
    sets_of_four = []
    with open(fname) as r:
        four_leaves = []
        for record in SeqIO.parse(r, "fasta"):
            s = str(record.seq)
            if '_________' in record.seq:
                s = s.replace('_________', '')
                four_leaves.append([int(x) for x in s.split('|')])
                sets_of_four.append(four_leaves)
                four_leaves = []
            else:
                four_leaves.append([int(x) for x in s.split('|')])
    print(len(sets_of_four), "trees read.")

    all_dists = dict()
    all_seqs = dict()
    
    for idx, seqs in enumerate(sets_of_four):
        dists = dict()
        seq_dict = dict()
        num_seqs = len(seqs)

        if idx == num_to_run:
                break

        for i in range(num_seqs):
            seq_dict[i] = seqs[i]
            
        for pair in itertools.combinations(list(range(num_seqs)), 2):

            a, b = pair
            x0 = np.random.uniform(0, 5.0, 3)
            x_star = optimize_len(m, seqs[a], seqs[b], x0) 
            # print the stddev of x_star
            d_ab = sorted(x_star, key=lambda x: x[0], reverse=True)[0][1]
            d_a, d_b, d_r = d_ab
            dists[(a, b)] = [d_a, d_b, d_r]

        all_dists[idx] = dists
        all_seqs[idx] = seq_dict

    return all_dists

def main():
    sets_of_four = []
    with open("MP_inconsistent/seqs_m10_k20.txt") as r:
        four_leaves = []
        for record in SeqIO.parse(r, "fasta"):
            s = str(record.seq)
            if '_________' in record.seq:
                s = s.replace('_________', '')
                four_leaves.append([int(x) for x in s.split('|')])
                sets_of_four.append(four_leaves)
                four_leaves = []
            else:
                four_leaves.append([int(x) for x in s.split('|')])
        
    # true tree
    t = "[&R] ((0:0.0360971597765934,1:3.339535381892265):0.0360971597765934,(2:0.0360971597765934,3:3.339535381892265):0.0360971597765934);"
    true_tree = dendropy.Tree.get(data=t, schema="newick")

    topo = '[&R] ((0,1),(2,3));'
    topo = dendropy.Tree.get(data=topo, schema="newick")
    global k 

    for k in (20, 30, 40, 50, 100, 200, 300, 400, 500, 1000, 5000):
        outdir = "/Users/gillianchu/raphael/repos/scmail/results_estbl/"
        # subprocess.run(["mkdir", "-p", outdir])

        est_dicts = run_mats("MP_inconsistent/seqs_m10_k{0}.txt".format(k), topo, 11, 1000)
        # print(est_dicts)
        print("Number of est trees", len(est_dicts))
        with open(outdir+"estdict_m10_k{0}_estbl.pkl".format(k), 'wb') as f:
            pickle.dump(est_dicts, f)
            # pd.DataFrame(mat).to_csv(outdir+"/mat"+str(idx)+".csv")
            # tree.write(path=outdir+"/tree"+str(idx)+".tre", schema="newick")

    # k=20
    # matdists = run_all("MP_inconsistent/seqs_m10_k20.txt", topo, 11, 1000)
    
    # print("filename: MP_inconsistent/seqs_m10_k5000.txt")
    # for x in out5000:
    #     print(x) 
    # print(true_tree)
    
    # print("filename: MP_inconsistent/seqs_m10_k20.txt")
    # for x in out20:
    #     print(x)
    # print(true_tree)

if __name__ == "__main__":
    main()

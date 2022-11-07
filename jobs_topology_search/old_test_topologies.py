import numpy as np
from scipy import optimize
import math
from math import log, exp
import sys
sys.path.append("/Users/gillianchu/raphael/repos/problin")
from problin_libs.ml import wrapper_felsenstein
from problin_libs.compare_two_trees import compare_trees
from sequence_lib import read_sequences

true_tree = ''
# x_star, f_star = wrapper_felsenstein(T, Q, msa, 3, root_edge_len=0.2)
# print(x_star, f_star)

# enumerate all fifteen topologies for 4 leaves
topologies = ["((a,b),(c,d));","((a,c),(b,d));","((a,d),(b,c));",
             "(a,(b,(c,d)));","(a,(c,(b,d)));","(a,(d,(b,c)));",
             "(b,(a,(c,d)));","(b,(c,(a,d)));","(b,(d,(a,c)));",
             "(c,(a,(b,d)));","(c,(b,(a,d)));","(c,(d,(a,b)));",
             "(d,(a,(b,c)));","(d,(b,(a,c)));","(d,(c,(a,b)));"]
# topologies = ["((a,b),(c,d));","(c,(b,(a,d)));"]

true_topology = '((a,b),(c,d));',

m = 10
for k in [200]: #30,40,50,100,200,300,400,500,1000,5000]:
    with open("ML_felsenstein_results_k" + str(k) + ".txt",'w') as fout:

        correct = 0
        total = 0

        S = read_sequences("/Users/gillianchu/raphael/repos/problin/MP_inconsistent/seqs_m10_k" + str(k) + ".txt")
        n_total = 0
        n_correct = 0
        for D in S[:100]:

            Q = []
            for i in range(k):
                q = {j+1:1/m for j in range(m)}
                q[0] = 0
                Q.append(q)

            top_x_star =  -float("inf")
            top_f_star =  -float("inf")
            top_topology = None
            for t, topology in enumerate(topologies): 
                f_star, est_tree, x_star = wrapper_felsenstein(topology, Q, D, use_log=True, optimize_branchlengths=True) #, init_tree=true_tree)
                print(t, f_star, est_tree, x_star)
                print("EST_TREE", est_tree)
                if f_star > top_f_star:
                    top_f_star = f_star
                    top_x_star = x_star
                    top_topology = topology
            if top_topology == true_topology:
                correct += 1
            total += 1
            fout.write(top_topology + "\t" + str(top_f_star) + "\n", flush=True)
            # print(k, "characters", correct, total, top_topology)
        # print(k, "characters", correct/total, "true topologies selected.")
        # fout.write(str(k) + " " + str(n_correct/n_total) + "\n")

    with open("ML_felsenstein_results_k" + str(k) + ".txt", "r") as r:
        lines = r.readlines()
        for line in lines:
            topology, likelihood = line.split("\t")
            print(compare_trees(true_topology, topology))

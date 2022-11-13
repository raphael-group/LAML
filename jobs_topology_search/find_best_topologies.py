import glob
import numpy as np
from scipy import optimize
import math
from math import log, exp
import sys
from problin_libs.ml import wrapper_felsenstein
from problin_libs.compare_two_trees import compare_trees
from problin_libs.sequence_lib import read_sequences

true_tree = ''

# enumerate all fifteen topologies for 4 leaves
topologies = ["((a,b),(c,d));","((a,c),(b,d));","((a,d),(b,c));",
             "(a,(b,(c,d)));","(a,(c,(b,d)));","(a,(d,(b,c)));",
             "(b,(a,(c,d)));","(b,(c,(a,d)));","(b,(d,(a,c)));",
             "(c,(a,(b,d)));","(c,(b,(a,d)));","(c,(d,(a,b)));",
             "(d,(a,(b,c)));","(d,(b,(a,c)));","(d,(c,(a,b)));"]

true_topology = '((a,b),(c,d));'

m = 10
correct = 0
total = 0



k = sys.argv[1]
reps = sys.argv[2]

for rep in range(int(reps)):
    dirname = "/n/fs/ragr-research/projects/problin/jobs_topology_search/ml_results_k" + str(k) + "_rep" + str(rep)
    # print(dirname)
    top_x_star =  -float("inf")
    top_f_star =  -float("inf")
    top_topology = None
    for idx, filename in enumerate(glob.iglob(f'{dirname}/topo*.txt')):
        # print(filename)
        with open(filename, "r") as r:
            lines = r.readlines()
            if len(lines) == 0:
                continue
            t, f_star, est_tree = lines[0].split('\t')
            x_star = lines[-1].split(',')

            # print(t)
            # print(f_star)
            # print(est_tree)
            # print(x_star)
            # print()

        if float(f_star) > top_f_star:
            top_f_star = float(f_star)
            top_x_star = x_star
            top_topology = t
    print(rep, top_topology, true_topology, top_f_star)
    if top_topology == true_topology:
        correct += 1
    total += 1

print("reps:", reps, "correct:", correct/total, correct, total)

    #with open("ML_felsenstein_results_k" + str(k) + ".txt", "r") as r:
    #    lines = r.readlines()
    #    for line in lines:
    #        topology, likelihood = line.split("\t")
    #        print(compare_trees(true_topology, topology))

import dendropy 
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
true_branches = [0.0360971597765934, 3.339535381892265, 0.0360971597765934, 0.0360971597765934, 3.339535381892265, 0.0360971597765934, 0.0]

m = 10
reps = 1000
ks = [20, 30, 40, 50, 100, 200, 300]
results = dict()
results["correct"] = []
results["branches"] = dict()

fname = "ML_fels_results.csv"
with open(fname, "w+") as f:
    f.write("method\tk\trep\ttopo\tbranch_idx\ttrue_branch\test_branch\n")
    for k in ks:
        correct = 0
        total = 0
        print("k:", k)
        branch_errors = dict()
        files_read = 0
        for rep in range(int(reps)):
            dirname = "/n/fs/ragr-research/projects/problin/jobs_topology_search/ml_results_k" + str(k) + "/ml_results_k" + str(k) + "_rep" + str(rep)
            # print(dirname)
            if rep % 100 == 0:
                print(str(rep) + "/" + str(reps))
            top_x_star =  -float("inf")
            top_f_star =  -float("inf")
            top_topology = None
            for idx, filename in enumerate(glob.iglob(f'{dirname}/topo*.txt')):
                with open(filename, "r") as r:
                    lines = r.readlines()
                    if len(lines) == 0:
                        continue
                    t, f_star, est_tree = lines[0].split('\t')
                    x_star = lines[-1].split(',')
                    nwkt = dendropy.Tree.get(data=est_tree, schema="newick", rooting="force-rooted")

                if float(f_star) > top_f_star:
                    top_f_star = float(f_star)
                    top_x_star = x_star
                    top_topology = t
                    est_branches = []
                    for e in nwkt.postorder_edge_iter():
                        est_branches.append(e.length)

            if top_topology == true_topology:
                correct += 1

            for bidx, estb in enumerate(est_branches):
                if bidx not in branch_errors:
                    branch_errors[bidx] = []
                # method\tk\trep\ttopo\t\tbranch_idx\ttrue_branch\test_branch
                f.write("ML_fels\t" + str(k) + "\t" + str(rep) + "\t" + str(top_topology) + "\t" + str(bidx) + "\t" + str(true_branches[bidx]) + "\t" + str(estb) + "\n")
                f.flush()
                branch_errors[bidx].append(abs(estb - true_branches[bidx]))
            
            total += 1

        print("reps:", reps, "correct:", correct/total, correct, total)
        print("files_read:", files_read/(int(reps)*len(topologies)), files_read, int(reps) * len(topologies))
        for bidx in branch_errors:
            print(bidx, np.mean(branch_errors[bidx]))

            if bidx not in results["branches"]:
                results["branches"][bidx] = []
            results["branches"][bidx].append(np.mean(branch_errors[bidx]))
        results["correct"].append(correct/total)

    with open("ML_fels_results.txt", "w+") as f:
        f.write("method k percent_correct" + str(bidx) + "\n")
        for kidx, v in enumerate(results["correct"]):
            f.write("ML_fels " + str(kidx) + " " + str(results["correct"][kidx]) + "\n")
    for bidx in results["branches"]:
        with open("ML_fels_branch" + str(bidx) + ".txt", "w+") as f:
            f.write("method k branch\n")
            for kidx, er in enumerate(results["branches"][bidx]):
                f.write("ML_fels " + str(ks[kidx]) + " " + str(er) + "\n")
            #with open("ML_felsenstein_results_k" + str(k) + ".txt", "r") as r:
            #    lines = r.readlines()
            #    for line in lines:
            #        topology, likelihood = line.split("\t")
            #        print(compare_trees(true_topology, topology))

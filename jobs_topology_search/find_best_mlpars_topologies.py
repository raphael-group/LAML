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


# k = sys.argv[1]
# reps = sys.argv[2]
reps = 1000

results = dict()
results["correct"] = []
results["branches"] = dict()
ks = [20, 30, 40, 50, 100, 200, 300, 400, 500, 5000]
for k in ks:
    correct = 0
    total = 0
    print("k:", k)
    branch_errors = dict()
    files_read = 0
    for rep in range(int(reps)):
    # for rep in range(100, 100 + int(reps)):
        dirname = "/n/fs/ragr-research/projects/problin/jobs_topology_search/mlpars_results_k" + str(k) + "/rep" + str(rep)
        top_ll =  -float("inf")
        top_topology = None
        if rep % 100 == 0:
            print(str(rep) + "/" + str(reps))
        for idx, filename in enumerate(glob.iglob(f'{dirname}/topo*.txt')):
            with open(filename, "r") as r:
                lines = r.readlines()
                if len(lines) == 0:
                    continue
                files_read += 1
                ll = lines[0]
                T = lines[1]
                branches = lines[2]
            idx = int(filename.split('/')[-1].split('.')[0].replace('topo',''))
            if float(ll) > top_ll:
                top_ll = float(ll)
                top_topology = topologies[idx]
                branches = [float(x) for x in branches.replace('[','').replace(']','').split(',')]
                for bidx, estb in enumerate(branches):
                    if bidx not in branch_errors:
                        branch_errors[bidx] = []
                    branch_errors[bidx].append(abs(abs(estb) - true_branches[bidx]))
                
        if top_topology == true_topology:
            correct += 1
        total += 1

    print("correct:", str(correct/total), correct, total)
    print("files_read:", files_read/(int(reps) * len(topologies)), files_read, int(reps) * len(topologies))

    for bidx in branch_errors:
        print(bidx, np.mean(branch_errors[bidx]))
        if bidx not in results["branches"]:
            results["branches"][bidx] = []
        results["branches"][bidx].append(np.mean(branch_errors[bidx]))
    results["correct"].append(correct/total)

with open("ML_pars_results.txt", "w+") as f:
    f.write("method k percent_correct" + str(bidx) + "\n")
    for kidx, v in enumerate(results["correct"]):
        f.write("ML_pars " + str(ks[kidx]) + " " + str(results["correct"][kidx]) + "\n")

for bidx in results["branches"]:
    with open("ML_pars_branch" + str(bidx) + ".txt", "w+") as f:
        f.write("method k branch\n")
        for kidx, er in enumerate(results["branches"][bidx]):
            f.write("ML_pars " + str(ks[kidx]) + " " + str(er) + "\n")

import os
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

dirname = "/n/fs/ragr-research/projects/problin/jobs_topology_search/run_test"
outfile = "/n/fs/ragr-research/projects/problin/jobs_topology_search/mean_var.txt"

d = dict()
# read all files in the directory
total_num = len(os.listdir(f'{dirname}'))
for idx, filename in enumerate(glob.iglob(f'{dirname}/slurm-*.out')):
# read first line
# read array of branches 
# calculate average variance of the branches for one replicate across the 20 restarts

        if idx % 10000 == 0:
            print(idx, '/', total_num)
        with open(filename, "r") as r: 
            lines = r.readlines()
            if len(lines) < 4: 
                continue
        topo_idx, k, rep_idx = lines[0].split()[-3:]
        # print(topo_idx, k, rep_idx)
        # print(filename)

        if len(lines) == 25:
            a = ''
            for line in lines[3:24]:
                l = line[:-1].rstrip()
                a += l # ','.join(l) # .decode()
                # a += line.replace('\t', '').rstrip()
        elif len(lines) == 23:
            a = ''
            for line in lines[1:22]:
                l = line[:-1].rstrip()
                a += l
        else:
            print(len(lines), "not handled.")
            break
        a += '] '
        # print(a)
        a = a.split('array')[1:]

        branches = [x[1:-3] for x in a] 
        branches = [[xx.strip('[]()').replace(',', '') for xx in x.split()] for x in branches]
        branches = [np.array([float(xx) for xx in x if xx != '']) for x in branches]
        m = np.matrix(np.array(branches), dtype=object)
        if m.shape != (20, 7):
             continue
        # print(m, m.shape)
        if k not in d:
            d[k] = dict()
        if topo_idx not in d[k]:
            d[k][topo_idx] = dict()
        if rep_idx not in d[k][topo_idx]:
            d[k][topo_idx][rep_idx] = [] 
        a = np.mean(m.T, axis=1) # axis of each row
        v = np.var(m.T, axis=1)
        d[k][topo_idx][rep_idx] = [np.array2string(a.flatten()).replace('\n', ''), np.array2string(v.flatten()).replace('\n', '')]
        # print("HERE", d[k][topo_idx][rep_idx])
        #if idx > 10:
        #    break

# output: # k topo rep avg variance
with open(outfile, "w+") as w:
    for k in d:
        for topo_idx in d[k]:
            for rep_idx in d[k][topo_idx]:
                w.write('\t'.join([k, topo_idx, rep_idx, d[k][topo_idx][rep_idx][0], d[k][topo_idx][rep_idx][1]]) + "\n")



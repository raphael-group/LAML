#! /usr/bin/env python
import sys, glob
from problin_libs.sequence_lib import read_sequences
from problin_libs.ml import wrapper_felsenstein
import numpy as np

# dirname="/n/fs/ragr-research/projects/problin/jobs/run_true_branches"
dirname = sys.argv[1]
output = sys.argv[2]

likelihoods = []
trees = []
branch_vectors = []

# get the last line and parse it into (likelihood, tree, branches)
for filename in glob.iglob(f'{directory}/*.logfile'):
    with open(filename) as f:
        for line in f:
            pass
        likelihood, tree, branches = line.split()
        print(likelihood, tree, branches)
        likelihoods.append(likelihood)
        trees.append(tree)
        branch_vectors.append(branches)
    break

runtimes = []
for filename in glob.iglob(f'{directory}/*.logrun'):
    with open(filename) as f:
        for line in f:
            runtimes.append(line) 

# if this works, then write the likelihood to an output file
with open(output + ".likelihood", "w+") as w:
    for likelihood in likelihoods:
        w.write(likelihood)
with open(output + ".tree", "w+") as w:
    for tree in trees:
        w.write(tree)
with open(output + ".branches", "w+") as w:
    for branches in branch_vectors:
        w.write(branches)
with open(output + ".runtime", "w+") as w:
    for sec in runtime:
        w.write(sec)


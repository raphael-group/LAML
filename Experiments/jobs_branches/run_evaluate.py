#! /usr/bin/env python
from os.path import exists
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
order = []
runtimes = []
# get the last line and parse it into (likelihood, tree, branches)
for idx, filename in enumerate(glob.iglob(f'{dirname}/*.logfile')):

    fname, ext = filename.split('.')
    fname2 = fname + ".logrun"
    file_exists = exists(fname2)
    print(filename, fname2)

    if file_exists:
        with open(fname2) as f:
            for line in f:
                runtimes.append(line) 

        with open(filename) as f:
            last_lines = f.readlines()[-2:]
            print(last_lines)

            likelihood, o = last_lines[0].split('[&R]')
            o += last_lines[1]
            tree, branch_vector = o.split('array')

            likelihood = likelihood[1:-3]
            tree = tree[:-5]
            branch_vector = branch_vector[:-2].replace('\n', '').replace('\t', '').replace(',', '')[2:-2]
            print("likelihood", likelihood)
            print("tree", tree)
            print("branch_vector", branch_vector)
            
            likelihoods.append(likelihood)
            trees.append(tree)
            branch_vectors.append(branch_vector)
            order.append(filename)
            print()
    #if idx == 2:
    #    break

print(len(likelihoods), len(trees), len(branch_vectors), len(order), len(runtimes))
# if this works, then write the likelihood to an output file
with open(output + ".likelihood", "w+") as w:
    for i, likelihood in enumerate(likelihoods):
        w.write(str(order[i]) + "\t" + str(likelihood) + "\n" )
with open(output + ".tree", "w+") as w:
    for i, tree in enumerate(trees):
        w.write(str(order[i]) + "\t" + str(tree) + "\n")
with open(output + ".branches", "w+") as w:
    for i, branches in enumerate(branch_vectors):
        w.write(str(order[i]) + f"\t{branches}\n")
with open(output + ".runtime", "w+") as w:
    for i, sec in enumerate(runtimes):
        w.write(str(order[i]) + "\t" + str(sec))



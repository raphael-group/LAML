#! /usr/bin/env python

from problin_libs.sim_lib import *
from problin_libs.sequence_lib import write_sequences
from os import mkdir

h = 8 # tree height --> 256 leaves
b = 0.2 # branch length
nreps = 20

treeStr = get_balanced_tree(h,b)
with open("model_tree.nwk",'w') as fout:
    fout.write(treeStr)
tree = read_tree_newick(treeStr)

K=[10,30,50,70,90]
m=10

for k in K:
    Q = []
    for i in range(k):
        q = {j+1:1/m for j in range(m)}
        Q.append(q)
    for i in range(nreps):
        d = "k"+str(k)+"_rep"+str(i+1).rjust(3,'0')
        mkdir(d)
        char_mtrx = simulate_seqs(tree,Q)
        write_sequences(char_mtrx,k,d+"/characters.txt")

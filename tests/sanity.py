#! /usr/bin/env python

from problin_libs.sequence_lib import read_sequences
from problin_libs.ML_solver import ML_solver
from problin_libs.ml_log import wrapper_felsenstein as wf_log
from problin_libs.ml import wrapper_felsenstein as wf
from sys import argv
from treeswift import *
import random

random.seed(a=1234)

print("---------- TEST 1 -----------")
T = "((a:1,b:1)1:1,c:1)2:1;"
tree = read_tree_newick(T)

msa = dict()
msa['a'] = [1, 1]
msa['b'] = [1, 1]
msa['c'] = [1, 1]

Q = []
for i in range(2):
    m_i = len(set(msa[x][i] for x in msa if msa[x][i] not in [0,"-"]))
    q = {j+1:1/m_i for j in range(m_i)}
    q[0] = 0
    Q.append(q)

mySolver = ML_solver(msa,Q,tree.newick())
print("optimal",mySolver.optimize(initials=1),mySolver.params.phi,mySolver.params.nu)
mySolver.params.tree.write_tree_newick("output_sanitytest1.nwk")


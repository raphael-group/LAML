#! /usr/bin/env python

from problin_libs.sequence_lib import read_sequences
from problin_libs.ML_solver import ML_solver
from problin_libs.ml_log import wrapper_felsenstein as wf_log
from problin_libs.ml import wrapper_felsenstein as wf
from sys import argv
from treeswift import *
import random

random.seed(a=1234)

msa = read_sequences("3724_NT_T1_character_matrix.txt",filetype="charMtrx",delimiter="\t")
with open("3724_NT_T1_tree_subsampled.nwk",'r') as f:
    T = f.read().strip()

tree = read_tree_newick(T)

k = len(msa[next(iter(msa.keys()))])
Q = []
for i in range(k):
    m_i = len(set(msa[x][i] for x in msa if msa[x][i] not in [0,"-"]))
    q = {j+1:1/m_i for j in range(m_i)}
    q[0] = 0
    Q.append(q)

mySolver = ML_solver(msa,Q,tree.newick())
print("optimal",mySolver.optimize(initials=1),mySolver.params.phi,mySolver.params.nu)
mySolver.params.tree.write_tree_newick("3724_NT_T1_tree_subsampled_brlen.nwk")

#mySolver1 = ML_solver(msa,Q,tree.newick())
#print("phi0",mySolver1.optimize(fixed_phi=0),mySolver1.params.phi,mySolver1.params.nu)
#mySolver1.params.tree.write_tree_newick("rand10_brlen_phi0.nwk")

#mySolver2 = ML_solver(msa,Q,tree.newick())
#print("nu0",mySolver2.optimize(fixed_nu=0),mySolver2.params.phi,mySolver2.params.nu)
#mySolver2.params.tree.write_tree_newick("rand10_brlen_nu0.nwk")

#print(mySolver.compute_llh(mySolver.params))
#print(wf_log(tree.newick(), Q, msa, optimize_branchlengths=False))
#print(wf(tree.newick(), Q, msa, optimize_branchlengths=False))

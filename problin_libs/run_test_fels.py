import sys
import time
from treeswift import *
# from math import log,exp
from ml_log import wrapper_felsenstein as wf_log
from ml import wrapper_felsenstein as wf
# from ml3 import wrapper_felsenstein as wf2
from sequence_lib import read_sequences
from distance_based_lib import ML_pairwise_estimate
import numpy as np


k = int(sys.argv[1])
m = int(sys.argv[2])
n = int(sys.argv[3])

Q = []
for i in range(k):
    q = {j+1:1/m for j in range(m)}
    #q[0] = 1/m
    Q.append(q)

print("TEST (binary) n:", n)
# test_idx += 1
m = 10
k = 50
T = read_tree("/n/fs/ragr-research/projects/problin/MP_inconsistent/n" + str(n) + "_m" + str(m) + ".bintre", schema="newick")
T = T.newick()
S = read_sequences("../MP_inconsistent/seqs_n" + str(n) + "_m" + str(m) + "_k" + str(k) + "bintre.txt")
msa = S[0]
start_time = time.time()
print("log", wf_log(T, Q, msa, optimize_branchlengths=True, initials=20))
print("-- %s seconds -- " % (time.time() - start_time))
start_time = time.time()
print("og", wf(T, Q, msa, optimize_branchlengths=True, initials=20))
print("-- %s seconds -- " % (time.time() - start_time))


import dendropy
import time
from treeswift import *
# from math import log,exp
from ml_log import mlpars 
from ml_log import wrapper_felsenstein as wf_log
from ml import wrapper_felsenstein as wf
# from ml3 import wrapper_felsenstein as wf2
from sequence_lib import read_sequences
from distance_based_lib import ML_pairwise_estimate
import numpy as np

#T = "(a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934;"
#msa = {'a':[1],'b':[0],'c':[0],'d':[0]}
# T = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,d:3.339535381892265)f:0.0360971597765934)r:0.0;"
# T = read_tree("/n/fs/ragr-research/projects/problin/trial_simulated_tree.tree", schema="newick")


k= 10 
m = 10 
Q = []
for i in range(k):
    q = {j+1:1/m for j in range(m)}
    #q[0] = 1/m
    Q.append(q)

print("TEST 1")
T = "((a,b)e,(c,d)f)r;"
nwkt = dendropy.Tree.get(data=T, schema="newick", rooting="force-rooted")
print(T)
S = read_sequences("../MP_inconsistent/seqs_m10_k" + str(k) + ".txt")
msa = S[0]
print(msa)
print("mlpars", mlpars(nwkt, Q, msa)) #sankoff(nwkt, Q, msa, 0))



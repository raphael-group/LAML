import os
import numpy as np
from scipy import optimize
import math
from math import log, exp
import sys
from problin_libs.ml_log import mlpars 
from problin_libs.sequence_lib import read_sequences
import dendropy

true_tree = ''

# enumerate all fifteen topologies for 4 leaves
topologies = ["((a,b),(c,d));","((a,c),(b,d));","((a,d),(b,c));",
             "(a,(b,(c,d)));","(a,(c,(b,d)));","(a,(d,(b,c)));",
             "(b,(a,(c,d)));","(b,(c,(a,d)));","(b,(d,(a,c)));",
             "(c,(a,(b,d)));","(c,(b,(a,d)));","(c,(d,(a,b)));",
             "(d,(a,(b,c)));","(d,(b,(a,c)));","(d,(c,(a,b)));"]
# topologies = ["((a,b),(c,d));","(c,(b,(a,d)));"]

true_topology = '((a,b),(c,d));',

m = 10

idx = int(sys.argv[1])
k = int(sys.argv[2])
rep = int(sys.argv[3])

Q = []
for i in range(k):
    q = {j+1:1/m for j in range(m)}
    q[0] = 0
    Q.append(q)

S = read_sequences("/n/fs/ragr-research/projects/problin/MP_inconsistent/seqs_m10_k" + str(k) + ".txt")
topology = topologies[idx]

D = S[rep]
nwkt = dendropy.Tree.get(data=topology, schema="newick", rooting="force-rooted")


dirname = "mlpars_results_k" + str(k) + "_rep" + str(rep)
try:
    os.makedirs(dirname)
except FileExistsError:
    pass
with open(dirname + "/topo" + str(idx) + ".txt",'w') as fout:
    ll = mlpars(nwkt, Q, D)
    fout.write(str(ll) + "\n")

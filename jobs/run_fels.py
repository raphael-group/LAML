#! /usr/bin/env python
import sys
from problin_libs.sequence_lib import read_sequences
from problin_libs.ml import wrapper_felsenstein

# ./run_test_fels.sh ${k} ${m} ${i} ${initials} ${t}

k= int(sys.argv[1]) ##5000
m = int(sys.argv[2]) #10
idx = int(sys.argv[3])
initials = int(sys.argv[4])
tree = str(sys.argv[5])

print("k:", k, "m:", m, "idx:", idx, "initials:", initials, "tree:", tree)
Q = []
for i in range(k):
    q = {j+1:1/m for j in range(m)}
    q[0] = 1/m
    Q.append(q)

#T = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,d:3.339535381892265)f:0.0360971597765934)r:0.0;"

T = tree
S = read_sequences("../MP_inconsistent/seqs_m" + str(m) + "_k" + str(k) + ".txt")
print("len S", len(S))
msa = S[idx]

print(wrapper_felsenstein(T, Q, msa, use_log=True, optimize_branchlengths=True,initials=initials))

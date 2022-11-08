# from math import log,exp
from ml import wrapper_felsenstein as wf
# from ml3 import wrapper_felsenstein as wf2
from sequence_lib import read_sequences
from distance_based_lib import ML_pairwise_estimate
import numpy as np

k= 200
m = 10
Q = []
for i in range(k):
    q = {j+1:1/m for j in range(m)}
    #q[0] = 1/m
    Q.append(q)
#T = "(a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934;"
#msa = {'a':[1],'b':[0],'c':[0],'d':[0]}
T = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,d:3.339535381892265)f:0.0360971597765934)r:0.0;"

#T = "((a,b)e,(c,d)f)r;"
S = read_sequences("../MP_inconsistent/seqs_m10_k" + str(k) + ".txt")
msa = S[0]

print("TEST 1")
# print(wf(T, Q, msa, use_log=True, optimize_branchlengths=True, initials=20))
# print(wf2(T, Q, msa, optimize_branchlengths=True, initials=20))

print(wf(T, Q, msa, optimize_branchlengths=True, initials=20))


print("TEST 2")
T = "((a,b)e,(c,d)f)r;"
# print(wf(T, Q, msa, use_log=True, optimize_branchlengths=True, initials=20))
# print(wf2(T, Q, msa, optimize_branchlengths=True, initials=20))
print(wf(T, Q, msa, optimize_branchlengths=True, initials=20))

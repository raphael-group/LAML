import time
from treeswift import *
# from math import log,exp
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


k= 50 
m = 10 
Q = []
for i in range(k):
    q = {j+1:1/m for j in range(m)}
    #q[0] = 1/m
    Q.append(q)

#print("TEST 1")
#T = "((a,b)e,(c,d)f)r;"
#S = read_sequences("../MP_inconsistent/seqs_m10_k" + str(k) + ".txt")
#msa = S[0]

#print("log", wf_log(T, Q, msa, optimize_branchlengths=True, initials=20))
#print("og", wf(T, Q, msa, optimize_branchlengths=True, initials=20))

k = 50 
m = 10 
Q = []
for i in range(k):
    q = {j+1:1/m for j in range(m)}
    #q[0] = 1/m
    Q.append(q)

T = "((a:1,b:1)1:1,c:1)2:1;"
msa = dict()
msa['a'] = [1]
msa['b'] = [1]
msa['c'] = [1]
print("og", wf(T, Q, msa, use_log=False, optimize_branchlengths=False))
print("log", wf_log(T, Q, msa, optimize_branchlengths=False))

# print("TEST 2")
#T = "((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,d:3.339535381892265)f:0.0360971597765934)r:0.0;"
#S = read_sequences("../MP_inconsistent/seqs_m10_k" + str(k) + ".txt")
#msa = S[0]
#print("log", wf_log(T, Q, msa, optimize_branchlengths=False, initials=20))
#print("og", wf(T, Q, msa, optimize_branchlengths=False, initials=20))
'''
k = 30 
print("TEST 3")
T = read_tree("/n/fs/ragr-research/projects/problin/MP_inconsistent/n10_m10.tre", schema="newick")
T = T.newick()
S = read_sequences("../MP_inconsistent/seqs_n10_m10_k" + str(k) + ".txt")
msa = S[0]
start_time = time.time()
print("log", wf_log(T, Q, msa, optimize_branchlengths=False, initials=20))
print("-- %s seconds -- " % (time.time() - start_time))
start_time = time.time()
print("og", wf(T, Q, msa, optimize_branchlengths=False, initials=20))
print("-- %s seconds -- " % (time.time() - start_time))

print("TEST 4")
start_time = time.time()
print("log", wf_log(T, Q, msa, optimize_branchlengths=True, initials=20))
print("-- %s seconds -- " % (time.time() - start_time))
start_time = time.time()
print("og", wf(T, Q, msa, optimize_branchlengths=True, initials=20))
print("-- %s seconds -- " % (time.time() - start_time))

def write_out(outfile, ll, t, bl, tfile, sfile, cmd):
    with open(outfile, "w+") as w:
        w.write(tfile)
        w.write(sfile)
        w.write(cmd)
        w.write(ll)
        w.write(t)
        w.write(bl)
test_idx = 5
# for n in [25, 30, 50, 70, 100]: 
for n in [70]: #50, 70, 100]: 
    print("TEST " + str(test_idx), "n:", n)
    # test_idx += 1
    m = 10
    k = 50
    outfile = "n" + str(n) + "_m" + str(m) + "_k" + str(k) + "optbl.loglikelihood"

    T = read_tree("/n/fs/ragr-research/projects/problin/MP_inconsistent/n" + str(n) + "_m" + str(m) + ".tre", schema="newick")
    T = T.newick()
    S = read_sequences("../MP_inconsistent/seqs_n" + str(n) + "_m" + str(m) + "_k" + str(k) + ".txt")
    msa = S[0]
    start_time = time.time()
    ll, t, bl = wf_log(T, Q, msa, optimize_branchlengths=True, initials=20)

    print("-- %s seconds -- " % (time.time() - start_time))
    print("log", ll, t, bl)
    tfile = "/n/fs/ragr-research/projects/problin/MP_inconsistent/n" + str(n) + "_m" + str(m) + ".tre"
    sfile = "../MP_inconsistent/seqs_n" + str(n) + "_m" + str(m) + "_k" + str(k) + ".txt"
    cmd = "wf_log(T, Q, msa, optimize_branchlengths=True, initials=20))"
    write_out(outfile, ll, t, bl, tfile, sfile, cmd)
    
    start_time = time.time()
    l, t, bl = wf(T, Q, msa, optimize_branchlengths=True, initials=20)
    print("-- %s seconds -- " % (time.time() - start_time))
    
    outfile = "n" + str(n) + "_m" + str(m) + "_k" + str(k) + "optbl.likelihood"
    
    print("og", l, t, bl)
    tfile = "/n/fs/ragr-research/projects/problin/MP_inconsistent/n" + str(n) + "_m" + str(m) + ".tre"
    sfile = "../MP_inconsistent/seqs_n" + str(n) + "_m" + str(m) + "_k" + str(k) + ".txt"
    cmd = "wf(T, Q, msa, optimize_branchlengths=True, initials=20))"
    write_out(outfile, ll, t, bl, tfile, sfile, cmd)

test_idx = 6
for n in [30, 50, 70, 100]: 
    print("TEST " + str(test_idx), "(binary) n:", n)
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
'''

# T = "((a,b)e,(c,d)f)r;"
## print(wf(T, Q, msa, use_log=True, optimize_branchlengths=True, initials=20))
## print(wf2(T, Q, msa, optimize_branchlengths=True, initials=20))
# print(wf(T, Q, msa, optimize_branchlengths=True, initials=20))'''

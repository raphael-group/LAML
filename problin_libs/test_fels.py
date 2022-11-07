from ml import wrapper_felsenstein
from sequence_lib import read_sequences
from distance_based_lib import ML_pairwise_estimate

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
#for msa in S:
print(wrapper_felsenstein(T, Q, msa, use_log=True, optimize_branchlengths=True,initials=1))

#print(wrapper_felsenstein(T, Q, msa, use_log=True, optimize_branchlengths=True,initials=20))
# three-way
#t = ML_pairwise_estimate(msa['a'],msa['b'],Q,x0=[0.0360971597765934,3.339535381892265,0.0360971597765934],do_optimize=False)
#t = ML_pairwise_estimate(msa['a'],msa['b'],Q,do_optimize=True)
#t = ML_pairwise_estimate(msa['a'],msa['b'],Q,x0=[0.69317159,0.69357885,0.34657359],do_optimize=False)
#print(t)

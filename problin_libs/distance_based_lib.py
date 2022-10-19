#! /usr/bin/env python
from math import log,exp
from sequence_lib import read_sequences
from scipy import optimize
from random import random

q = 0.1

def simple_pairwise_estimate(a,b):
    za = 0
    zb = 0
    z = 0
    for ca,cb in zip(a,b):
        za += (ca == 0)	
        zb += (cb == 0)
        z += (ca == 0 or cb == 0 or ca != cb)	
    da = -log(za/z)	
    db = -log(zb/z)
    return da,db

def ML_pairwise_estimate(a,b,initials=20):
    def __sets__(seq_a, seq_b):
        # get the msa
        k = len(seq_a)
        
        ## calculate the sets  
        s_0, s_1a, s_1b, s_2, s_3 = set(), set(), set(), set(), set()
        for idx in range(len(seq_a)):
            c_a, c_b = seq_a[idx], seq_b[idx]
            if c_a == c_b:
                if c_a == 0:
                    s_0.add(idx)
                else:
                    s_2.add(idx)
            elif c_a == 0: # then c_b != 0
                s_1b.add(idx)
            elif c_b == 0: # then c_a != 0
                s_1a.add(idx)
            else:
                s_3.add(idx)
        
        assert len(s_0) + len(s_1a) + len(s_1b) + len(s_2) + len(s_3) == k
        return [s_0, s_1a, s_1b, s_2, s_3,k]

    def l(x): # negative log-likelihood
        d_a, d_b, d_r = x
        
        p1 = -(s_1b_len + s_0_len) * d_a + (s_1a_len + s_3_len) * log(1 - exp(-d_a))
        p2 = -(s_1a_len + s_0_len) * d_b + (s_1b_len + s_3_len) * log(1 - exp(-d_b)) - (k - s_2_len) * d_r
        p3 = 0.0
        
        for i in range(s_2_len): # assuming that prob to any alpha is the same
            # q_ialpha is the transition rate at site i from 0 -> character at node a at site i
            # iterate over sites in s_2, get q for each alpha
            p3 += log(q**2 * (1 - exp(-d_a)) * (1 - exp(-d_b)) * exp(-d_r) + q*(1 - exp(-d_r)))
        
        return -(p1 + p2 + p3)
    
    x_star = []
    s_0, s_1a, s_1b, s_2, s_3, k = __sets__(a, b)
    s_0_len,s_1a_len,s_1b_len,s_2_len,s_3_len = len(s_0), len(s_1a), len(s_1b), len(s_2), len(s_3)
    dmax = -log(1/k)*2
    dmin = -log(1-1/k)/2
    bound = (dmin,dmax)
    x_star = None
    f_star = float("inf")
    for i in range(initials):
        x0 = (random()*5,random()*5,random()*5)
        out = optimize.minimize(l, x0, method="SLSQP", options={'disp': False,'maxiter':1000}, bounds=[bound,bound,bound])
        if out.success and out.fun < f_star:
            x_star = out.x
            f_star = out.fun 
    return x_star,f_star

def triplet_estimate_old(d_ab,d_ac,d_bc):
# Note: d_xy = (d(x,LCA(x,y)),d(y,LCA(x,y)))
    delta_a = abs(d_ab[0]-d_ac[0])
    delta_b = abs(d_ab[1]-d_bc[0])
    delta_c = abs(d_ac[1]-d_bc[1])
    delta = min(delta_a,delta_b,delta_c)
    if delta == delta_a:
        return 1 # a|bc
    if delta == delta_b:
        return 2 # b|ac
    else:
        return 3 # c|ab        

def triplet_estimate(dr_ab,dr_ac,dr_bc):
# test the topology of the triplet (a,b,c)
    dr_max = max(dr_ab,dr_bc,dr_ac)
    if dr_max == dr_bc:
        return 1 # a|bc
    if dr_max == dr_ac:
        return 2 # b|ac
    return 3 # c|ab        

'''
S = read_sequences("../MP_inconsistent/seqs_m10_k40.txt")
count = 0
i = 0
for D in S:
    a = D['a']
    b = D['b']
    c = D['c']
    d = D['d']

    dr_ab = ML_pairwise_estimate(a,b)[0][2]
    dr_ac = ML_pairwise_estimate(a,c)[0][2]
    dr_bc = ML_pairwise_estimate(b,c)[0][2]
    t_abc = triplet_estimate(dr_ab,dr_ac,dr_bc)

    dr_ad = ML_pairwise_estimate(a,d)[0][2]
    dr_cd = ML_pairwise_estimate(c,d)[0][2]
    t_acd = triplet_estimate(dr_ac,dr_ad,dr_cd)
    
    i += 1
    count += (t_abc == 3 and t_acd == 1)
    print(i,count)
'''    

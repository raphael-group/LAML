#! /usr/bin/env python
from math import log
from sequence_lib import read_sequences

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

def triplet_estimate(d_ab,d_ac,d_bc):
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

S = read_sequences("seqs_k1000.txt")
for D in S:
    a = D['a']
    b = D['b']
    c = D['c']

    d_ab = simple_pairwise_estimate(a,b)
    d_ac = simple_pairwise_estimate(a,c)
    d_bc = simple_pairwise_estimate(b,c)
    #print(d_ab,d_ac,d_bc)
    print(triplet_estimate(d_ab,d_ac,d_bc))

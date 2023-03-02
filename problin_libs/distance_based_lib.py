#! /usr/bin/env python
from math import log,exp
from problin_libs.sequence_lib import read_sequences, read_Q
from scipy import optimize
from random import random,seed
from treeswift import *

def ML_pairwise_estimate(a,b,Q,initials=10,do_optimize=True,x0=None):
# Note: optimize = False must be used with x0 not None to compute the likelihood
    def __sets__(seq_a, seq_b):
        # get the msa
        k = len(seq_a)
        s0 = s1a = s1b = s3 = 0
        m0a = m0b = m1a = m1b = 0
        S2 = []
        for idx in range(k):
            c_a, c_b = seq_a[idx], seq_b[idx]         
            s0 += int(c_a == c_b == 0)
            s1a += int(c_a > 0 and c_b == 0)
            s1b += int(c_b > 0 and c_a == 0)
            s3 += int(c_a > 0 and c_b > 0 and c_a != c_b)
            m0a += int(c_a == 0 and c_b == -1)
            m0b += int(c_b == 0 and c_a == -1)
            m1a += int(c_a > 0 and c_b == -1)
            m1b += int(c_b > 0 and c_a == -1)
            if c_a > 0 and c_a == c_b:
                S2.append(idx)
        return s0,s1a,s1b,S2,s3,m0a,m0b,m1a,m1b

    def l(x): # negative log-likelihood
        d_a, d_b, d_r = x
        
        p1 = -(s1b + s0) * d_a + (s1a + s3 + m0a) * log(1 - exp(-d_a))
        p2 = -(s1a + s0) * d_b + (s1b + s3 + m0b) * log(1 - exp(-d_b)) - (s0+s1a+s1b+s3+m0a+m0b) * d_r
        p3 = 0.0        
        for i in S2: # assuming that prob to any alpha is the same
            # q is the transition rate at site i from 0 -> character at node a at site i
            # iterate over sites in s_2, get q for each 
            q = Q[i][a[i]]
            p3 += log(q**2 * (1 - exp(-d_a)) * (1 - exp(-d_b)) * exp(-d_r) + q*(1 - exp(-d_r)))        
        p4 = m1a*log(1-exp(-d_a-d_r))
        p5 = m1b*log(1-exp(-d_b-d_r))
        return -(p1 + p2 + p3 + p4 + p5)
    
    s0,s1a,s1b,S2,s3,m0a,m0b,m1a,m1b = __sets__(a,b)
    s2 = len(S2)    
    if do_optimize:
        k = len(a)
        dmax = -log(1/k)*2
        dmin = -log(1-1/k)/2
        bound = (dmin,dmax)
        x_star = None
        f_star = float("inf")
        if x0 is not None:
            out = optimize.minimize(l, x0, method="SLSQP", options={'disp': False,'maxiter':1000}, bounds=[bound,bound,bound])
            if out.success and out.fun < f_star:
                x_star = out.x
                f_star = out.fun
        else:        
            for i in range(initials):
                x0 = (random()*5+dmin,random()*5+dmin,random()*5+dmin)
                out = optimize.minimize(l, x0, method="SLSQP", options={'disp': False,'maxiter':1000}, bounds=[bound,bound,bound])
                if out.success and out.fun < f_star:
                    x_star = out.x
                    f_star = out.fun
    else:
        x_star = x0
        f_star = l(x_star)
                    
    return x_star,f_star

def triplet_estimate(dr_ab,dr_ac,dr_bc):
# test the topology of the triplet (a,b,c)
    dr_max = max(dr_ab,dr_bc,dr_ac)
    if dr_max == dr_bc:
        return 1 # a|bc
    if dr_max == dr_ac:
        return 2 # b|ac
    return 3 # c|ab       

def greedy_triplet(sequences,Q):
# sequences is a dictionary mapping a name to a sequence 
    def __get_key__(a,b):
        return (a,b) if a<=b else (b,a)
    def __compute_Dr__():
        seq_names = list(sequences.keys())
        N = len(seq_names)
        Dr = {}
        for i in range(N-1):
            a = seq_names[i]
            sa = sequences[a]
            for j in range(i+1,N):
                b = seq_names[j]
                sb = sequences[b]
                key_ab = __get_key__(a,b)
                Dr[key_ab] = ML_pairwise_estimate(sa,sb,Q)[0][2]
        return Dr

    def __add_one_seq__(c,sc,root_node,Dr):
    # find a branch to add c to the tree by checking the triplets
    # return the node below the branch
        if root_node.is_leaf():
            return root_node
        left_child, right_child = root_node.children # assuming a perfect binary tree
        left_side = list(x.label for x in left_child.traverse_leaves())
        right_side = list(x.label for x in right_child.traverse_leaves())
        
        a = left_side[0]
        dr_ac = Dr[__get_key__(a,c)]
        for x in left_side[1:]:
            dr_xc = Dr[__get_key__(x,c)]
            if dr_xc > dr_ac:
                a = x
                dr_ac = dr_xc
        
        b = right_side[0]
        dr_bc = Dr[__get_key__(b,c)]
        for x in right_side[1:]:
            dr_xc = Dr[__get_key__(x,c)]
            if dr_xc > dr_bc:
                b = x
                dr_bc = dr_xc
        dr_ab = Dr[__get_key__(a,b)]
        trpl_abc = triplet_estimate(dr_ab,dr_ac,dr_bc)
        
        if trpl_abc == 1: # a|bc --> recurse on the right side of the tree
            return __add_one_seq__(c,sc,right_child,Dr)
        elif trpl_abc == 2: # b|ac --> recurse on the left side of the tree
            return __add_one_seq__(c,sc,left_child,Dr)
        else: # c|ab --> c is an outgroup of root_node
            return root_node
    def __add_node__(c,T,v):
        # add c on the branch above the specified v
        w = Node()
        w.label = c
        u = Node()
        if v.is_root():
            T.root = u
        else:    
            p = v.parent
            p.remove_child(v)
            p.add_child(u)
        u.add_child(v)
        u.add_child(w)
                    
    Dr = __compute_Dr__() # mapping a pair of sequence names to their dr
    seq_names = list(sequences.keys())
    a,b = seq_names[:2]
    sa,sb = sequences[a],sequences[b]
    tree = read_tree_newick("("+a+","+b+");")
    
    for c in seq_names[2:]:
        sc = sequences[c]
        v = __add_one_seq__(c,sc,tree.root,Dr)
        __add_node__(c,tree,v)
    return tree.newick()

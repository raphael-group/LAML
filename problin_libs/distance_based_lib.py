#! /usr/bin/env python
from math import log,exp
from sequence_lib import read_sequences
from scipy import optimize
from random import random,seed
from treeswift import *

q = 0.1

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

def triplet_estimate(dr_ab,dr_ac,dr_bc):
# test the topology of the triplet (a,b,c)
    dr_max = max(dr_ab,dr_bc,dr_ac)
    if dr_max == dr_bc:
        return 1 # a|bc
    if dr_max == dr_ac:
        return 2 # b|ac
    return 3 # c|ab       

def greedy_triplet(sequences):
# sequences is a dictionary mapping a name to a sequence 
    def __get_key__(a,b):
        return (a,b) if a<=b else (b,a)
    def __add_one_seq__(c,sc,root_node,Dr,sequences):
    # find a branch to add c to the tree by greedily checking the triplets
    # return the node below the branch
        if root_node.is_leaf():
            return root_node
        left_child, right_child = root_node.children # assuming a perfect binary tree
        a = list(x.label for x in left_child.traverse_leaves())[0]
        b = list(x.label for x in right_child.traverse_leaves())[0]
        sa = sequences[a]
        sb = sequences[b]
        key_ab = __get_key__(a,b)
        key_ac = __get_key__(a,c)
        key_bc = __get_key__(b,c)
        if key_ab in Dr:
            dr_ab = Dr[key_ab] 
        else:    
            dr_ab = ML_pairwise_estimate(sa,sb)[0][2]    
            Dr[key_ab] = dr_ab
        if key_ac in Dr:
            dr_ac = Dr[key_ac] 
        else:    
            dr_ac = ML_pairwise_estimate(sa,sc)[0][2]    
            Dr[key_ac] = dr_ac
        if key_bc in Dr:
            dr_bc = Dr[key_bc] 
        else:    
            dr_bc = ML_pairwise_estimate(sb,sc)[0][2]    
            Dr[key_bc] = dr_bc
        trpl_abc = triplet_estimate(dr_ab,dr_ac,dr_bc)
        
        if trpl_abc == 1: # a|bc --> recurse on the right side of the tree
            return __add_one_seq__(c,sc,right_child,Dr,sequences)
        elif trpl_abc == 2: # b|ac --> recurse on the left side of the tree
            return __add_one_seq__(c,sc,left_child,Dr,sequences)
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
                    
    Dr = {} # mapping a pair of sequence names to their dr
    seq_names = list(sequences.keys())
    a,b = seq_names[:2]
    sa,sb = sequences[a],sequences[b]
    tree = read_tree_newick("("+a+","+b+");")
    
    for c in seq_names[2:]:
        sc = sequences[c]
        v = __add_one_seq__(c,sc,tree.root,Dr,sequences)
        __add_node__(c,tree,v)
    return tree.newick()

'''
S = read_sequences("../MP_inconsistent/seqs_m10_k1000.txt")
i = 0
for D in S:
    print(greedy_triplet(D))
    i += 1
    if i >= 100:
        break
'''        

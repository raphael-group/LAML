#! /usr/bin/env python

from treeswift import *
from math import *
from random import random
from scipy.optimize import fsolve

def f(p,q):
    return 2*q*p**4 - q**2*p**2 + 2*(1-q**2)*p + 2*q**2-2

def root_f(q):
    def g(p):
        return f(p,q)
    sol = fsolve(g,1)
    return sol[-1] 

def simulate_tree(q):
    p = (1+root_f(q))/2
    d_b = d_d = -log(1-p)
    d_a = d_c = d_e = d_f = -log(p)

    E = {'a':d_a,'b':d_b,'c':d_c,'d':d_d,'e':d_e,'f':d_f,'r':None}
    T = read_tree_newick("((a,b)e,(c,d)f)r;")

    for node in T.traverse_preorder():
        node.edge_length = E[node.label]
    return T

def simulate_seqs(tree,k,m):
    tree.root.seq = [0]*k
    for node in tree.traverse_preorder():
        if node.is_root():
            continue
        p = 1-exp(-node.edge_length)
        pnode = node.parent
        seq = []
        for c in pnode.seq:
            if c != 0:
                seq += [c]
            else:
                r = random()
                if r < p: # change state
                    seq += [int(random()*m) + 1] 
                else: # keep 0
                    seq += [0]  
        node.seq = seq

def compute_MP(T,seqs):
    score = 0
    for node in T.traverse_postorder():
        if node.is_leaf():
            node.seq = seqs[node.label]
        else:
            a,b = node.children
            s1 = a.seq
            s2 = b.seq
            s = [x if x==y and x!=0 else 0 for x,y in zip(s1,s2)]     
            node.seq = s
            curr_scores = [int(x!=z) + int(y!=z) for (x,y,z) in zip(s1,s2,s)]
            curr_score = sum(curr_scores)
            score += curr_score
    return score


treeList = ["((a,b),(c,d));","((a,c),(b,d));","((a,d),(b,c));",
            "(a,(b,(c,d)));","(a,(c,(b,d)));","(a,(d,(b,c)));",
            "(b,(a,(c,d)));","(b,(c,(a,d)));","(b,(d,(a,c)));",
            "(c,(a,(b,d)));","(c,(b,(a,d)));","(c,(d,(a,b)));",
            "(d,(a,(b,c)));","(d,(b,(a,c)));","(d,(c,(a,b)));"]
reps = 1000

m = 10
q = 1/m # collision probability
T = simulate_tree(q)
T.write_tree_newick("m" + str(m) + ".tre")

for k in [20,30,40,50,100,200,300,400,500,1000,5000]:
    count = 0
    outfile = "seqs_m" + str(m) + "_k" + str(k) + ".txt"
    with open(outfile,'w') as fout:
        for i in range(reps):
            simulate_seqs(T,k,m)
            seqs = {}
            for node in T.traverse_leaves():
                fout.write(">" + node.label + "\n")
                fout.write("|".join([str(x) for x in node.seq])+"\n")
                seqs[node.label] = node.seq
            fout.write("_________\n")    
            MP_scores = []    
            for treeStr in treeList:
                T1 = read_tree_newick(treeStr)   
                MP_scores.append(compute_MP(T1,seqs))
            count += min(MP_scores[1:]) > MP_scores[0]
    print(k,count/reps)    

import numpy as np
import dendropy
import math
from random import random,seed
from math import log,exp
from scipy import optimize


def prob_same(node_likelihood, char, site, curr_node, use_log):
    # print("prob same")
    c = curr_node.child_nodes()[0]
    tp = math.exp(-get_branchlen(c))
    return tp * node_likelihood[c][site][char]

def prob_change(msa, q_dict, node_likelihood, site, curr_node, use_log):
    all_child_prob = 1.0
    for c in curr_node.child_nodes():
        char = get_char(msa, c, site)
        if char != 0: 
            q_ialpha = q_dict[site][char] 
            tp = q_ialpha * (1 - math.exp(-get_branchlen(c)))
            all_child_prob *= tp * node_likelihood[c][site][char]
        else:
            q_ialpha = q_dict[site][0]
            tp = math.exp(-get_branchlen(c))
            all_child_prob *= tp * node_likelihood[c][site][0]
    return all_child_prob



def likelihood_under_n(node_likelihood, n, site, msa, q_dict, use_log):
    child_states = set()
    # print(n)
    if n not in node_likelihood:
        node_likelihood[n] = dict()
        node_likelihood[n][site] = dict()
        
    child_states = []
    for child in n.child_nodes():
        if child.is_leaf():
            child_states.append(get_char(msa, child, site))
        else:
            for x in node_likelihood[child][site]:
                state_prob = node_likelihood[child][site][x]
                if state_prob > 0.0:
                    child_states.append(x)
                    
    parent_poss_states = dict()
    if 0 in child_states: # probability 0 -> 0
        if len(set(child_states)) == 1: # both children are state 0 
            parent_poss_states[0] = prob_same(node_likelihood, 0, site, n, use_log) ** 2
        else: 
            for c in child_states: # probability c -> c != 0
                parent_poss_states[c] = 0.0
            # probability 0 -> c (alpha)
            parent_poss_states[0] = prob_change(msa, q_dict, node_likelihood, site, n, use_log)  
    else:
        if len(set(child_states)) == 1: # both children are same nonzero state
            c = child_states[0]
            parent_poss_states[c] = 1.0 
            parent_poss_states[0] = prob_change(msa, q_dict, node_likelihood, site, n, use_log)
        else:
            parent_poss_states[0] = 1.0

    for x in parent_poss_states.keys():
        node_likelihood[n][site][x] = parent_poss_states[x]

    return node_likelihood

def get_branchlen(child_node):
    if child_node.edge_length is None:
        print(child_node.child_nodes())
    return child_node.edge_length

def get_char(msa, leaf_node, site):
    return msa[int(leaf_node.taxon.__str__().replace("'", ""))-1][site]

def felsenstein(T, Q, msa, root_edge_len=0.2, use_log=False):
    numsites = len(msa[0])

    alphabet = dict()
    for site in range(numsites):
        alphabet[site] = Q[site].keys()

    print("q_dict", Q)
    nwkt = dendropy.Tree.get(data=T, schema="newick")
    print(nwkt)

    for n in nwkt.leaf_node_iter():
        print(n.taxon, ''.join([str(get_char(msa, n, s)) for s in range(numsites)]))

    node_likelihood = dict()

    ## CALCULATE THE LIKELIHOOD
    for n in nwkt.postorder_node_iter():
        # print("node:", n)
        if n.taxon is not None: # must be a leaf node, set up 
            node_likelihood[n] = dict()
            for site in range(numsites):
                node_likelihood[n][site] = dict()
                for char in alphabet[site]:
                    node_likelihood[n][site][char] = 0.0
                char_state = get_char(msa, n, site)
                node_likelihood[n][site][char_state] = 1.0
            
        elif n.taxon is None: # must be an internal node
            for site in range(numsites):
                if n not in node_likelihood.keys():
                    node_likelihood[n] = dict()
                node_likelihood[n][site] = dict()
                for char in alphabet[site]:
                    node_likelihood[n][site][char] = 0.0
                
                node_likelihood = likelihood_under_n(node_likelihood, n, site, msa, Q, use_log)

    tree_likelihood = 1.0
    for site in range(numsites):
        for rootchar in node_likelihood[n][site].keys():
            prob_rootchar = node_likelihood[n][site][rootchar]
            print(rootchar, prob_rootchar)
            if prob_rootchar > 0.0: 
                q_ialpha = Q[site][rootchar]
                if rootchar == 0:
                    if use_log:
                        tree_likelihood += (-root_edge_len) + np.log(prob_rootchar) # + np.log(q_ialpha) 
                    else:
                        tree_likelihood *= (math.exp(-root_edge_len)) * prob_rootchar # * q_ialpha 
                    
                else:
                    if use_log:
                        tree_likelihood += np.log((1 - math.exp(-root_edge_len))) + np.log(q_ialpha) + np.log(prob_rootchar)
                    else:
                        tree_likelihood *= ((1 - math.exp(-root_edge_len)) * q_ialpha * prob_rootchar)
    
    return tree_likelihood


def wrapper_felsenstein(T, Q, msa, k, root_edge_len=0.2, use_log=False, initials=20):
    numsites = len(msa[0])

    alphabet = dict()
    for site in range(numsites):
        alphabet[site] = Q[site].keys()

    print("q_dict", Q)
    nwkt = dendropy.Tree.get(data=T, schema="newick")
    num_edges = len(list(nwkt.postorder_edge_iter()))
    print(nwkt)

    for n in nwkt.leaf_node_iter():
        print(n.taxon, ''.join([str(get_char(msa, n, s)) for s in range(numsites)]))

    def felsenstein(x):
        # x is a vector containing all branch lengths
        # map branch lengths to tree edges
        for i, e in enumerate(nwkt.postorder_edge_iter()):
            e.length = x[i]

        node_likelihood = dict()

        for n in nwkt.postorder_node_iter():
            if n.taxon is not None: # must be a leaf node, set up 
                node_likelihood[n] = dict()
                for site in range(numsites):
                    node_likelihood[n][site] = dict()
                    for char in alphabet[site]:
                        node_likelihood[n][site][char] = 0.0
                    char_state = get_char(msa, n, site)
                    node_likelihood[n][site][char_state] = 1.0
                
            elif n.taxon is None: # must be an internal node
                for site in range(numsites):
                    if n not in node_likelihood:
                        node_likelihood[n] = dict()
                    node_likelihood[n][site] = dict()
                    for char in alphabet[site]:
                        node_likelihood[n][site][char] = 0.0
                    
                    node_likelihood = likelihood_under_n(node_likelihood, n, site, msa, Q, use_log)

        tree_likelihood = 1.0
        for site in range(numsites):
            for rootchar in node_likelihood[n][site].keys():
                prob_rootchar = node_likelihood[n][site][rootchar]
                # print(rootchar, prob_rootchar)
                if prob_rootchar > 0.0: 
                    q_ialpha = Q[site][rootchar]
                    if rootchar == 0:
                        if use_log:
                            tree_likelihood += (-root_edge_len) + np.log(prob_rootchar) # + np.log(q_ialpha) 
                        else:
                            tree_likelihood *= (math.exp(-root_edge_len)) * prob_rootchar # * q_ialpha 
                        
                    else:
                        if use_log:
                            tree_likelihood += np.log((1 - math.exp(-root_edge_len))) + np.log(q_ialpha) + np.log(prob_rootchar)
                        else:
                            tree_likelihood *= ((1 - math.exp(-root_edge_len)) * q_ialpha * prob_rootchar)
        return tree_likelihood

    x_star = []
    dmax = -log(1/k)*2
    dmin = -log(1-1/k)/2
    bound = (dmin, dmax)
    x_star = None
    f_star = float("inf")

    for i in range(initials):
        x0 = [random()] * num_edges
        print(x0)
        out = optimize.minimize(felsenstein, x0, method="SLSQP", options={'disp':True,'maxiter':1000}, bounds=[bound]*num_edges)
        if out.success and out.fun < f_star:
            x_star = out.x
            f_star = out.fun

    return x_star, f_star
    
from problin_libs.sequence_lib import read_sequences
from treeswift import *
k=30
m=10
Q = []
for i in range(k):
    q = {j+1:1/m for j in range(m)}
    q[0] = 0
    Q.append(q)

S = read_sequences("../MP_inconsistent/seqs_m10_k" + str(k) + ".txt")
D = S[1]
print(D)

lb2num = {'a':1,'b':2,'c':3,'d':4}
tree = read_tree_newick("../MP_inconsistent/m10.tre")
for node in tree.traverse_leaves():
    node.label = lb2num[node.label]
T = tree.newick()

msa = []
for x in ['a','b','c','d']:
    seq = [y for y in D[x]]
    msa.append(seq)
msa = np.array(msa)
print(msa)

print("likelihood",felsenstein(T, Q, msa, root_edge_len=0, use_log=True))

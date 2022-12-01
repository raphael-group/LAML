import numpy as np
import dendropy
import math
from random import random,seed
from math import log,exp
from scipy import optimize
from problin_libs.sequence_lib import read_sequences

def prob_same(node_likelihood, char, site, curr_node):
    c1,c2 = curr_node.child_nodes()
    tp1, tp2 = -get_branchlen(c1), -get_branchlen(c2)
    # tp1,tp2 = exp(-get_branchlen(c1)),exp(-get_branchlen(c2))
    # print(node_likelihood[c1][site][char], node_likelihood[c2][site][char])
    if node_likelihood[c1][site][char] == 0 or node_likelihood[c2][site][char] == 0:
        return 0
    return exp(tp1 + log(node_likelihood[c1][site][char]) + tp2 + log(node_likelihood[c2][site][char]))
    # GCLOG: return tp1*node_likelihood[c1][site][char]*tp2*node_likelihood[c2][site][char]

def num_children(n):
    if n.is_leaf():
        return 1
    else:
        s = 0
        for c in n.child_nodes():
            s += num_children(c) 
        return s 

def prob_change(msa, q_dict, node_likelihood, site, curr_node, child_states):
    all_child_prob = 0.0
    if len(curr_node.child_nodes()) == 1:
        num_nodes = 0
        print(num_children(curr_node))
        print(curr_node)
    c1, c2 = curr_node.child_nodes()
    l1 = get_branchlen(c1)
    l2 = get_branchlen(c2)
    for char1 in node_likelihood[c1][site]:
        if node_likelihood[c1][site][char1] == 0:
            continue
        p1 = log(node_likelihood[c1][site][char1])
        # GCLOG: p1 = node_likelihood[c1][site][char1]
        p1 += -l1 if char1==0 else log(1-exp(-l1)) + log(q_dict[site][char1])
        # GCLOG: p1 *= exp(-l1) if char1==0 else (1-exp(-l1))*q_dict[site][char1]
        for char2 in node_likelihood[c2][site]:
            if node_likelihood[c2][site][char2] == 0:
                continue
            p2 = log(node_likelihood[c2][site][char2])
            # GCLOG: p2 = node_likelihood[c2][site][char2]
            p2 += -l2 if char2==0 else log(1-exp(-l2)) + log(q_dict[site][char2])
            # GCLOG: p2 *= exp(-l2) if char2==0 else (1-exp(-l2))*q_dict[site][char2]
            all_child_prob += exp(p1 + p2)
            # GCLOG: all_child_prob += p1*p2
    return all_child_prob        
    
def likelihood_under_n(node_likelihood, n, site, msa, q_dict, is_root):
    child_states = set()
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
    if 0 in set(child_states): # probability 0 -> 0
        if len(set(child_states)) == 1: # both children are state 0 
            tmp = prob_same(node_likelihood, 0, site, n)
            parent_poss_states[0] = tmp
        else: 
            parent_poss_states[0] = prob_change(msa, q_dict, node_likelihood, site, n, child_states)  
    else:
        if len(set(child_states)) == 1: # both children are same nonzero state
            c = child_states[0]
            parent_poss_states[c] = 1.0 
        parent_poss_states[0] = prob_change(msa, q_dict, node_likelihood, site, n, child_states)

    for x in parent_poss_states.keys():
        node_likelihood[n][site][x] = parent_poss_states[x]

    return node_likelihood

def get_branchlen(child_node):
    if child_node.edge_length is None:
        print(child_node.child_nodes())
    return child_node.edge_length

def get_char(msa, leaf_node, site):
    return msa[leaf_node.taxon.label][site]

def wrapper_felsenstein(T, Q, msa, initials=20, optimize_branchlengths=False, init_tree=None, print_all=False):
    numsites = len(msa[next(iter(msa.keys()))])
    alphabet = dict()
    for site in range(numsites):
        alphabet[site] = Q[site].keys()

    nwkt = dendropy.Tree.get(data=T, schema="newick", rooting="force-rooted")
    num_edges = len(list(nwkt.postorder_edge_iter()))
    
    def felsenstein(x): 
        for i, e in enumerate(nwkt.postorder_edge_iter()): # visit the descendents before visiting edge
            e.length = x[i]
        # check what the root edge length is. if none, set to 0.0
        root_edge_len = get_branchlen(nwkt.seed_node)
        if root_edge_len == None:
            root_edge_len = 0.0

        alphabet = dict()
        for site in range(numsites):
            alphabet[site] = Q[site].keys()

        node_likelihood = dict()

        for n in nwkt.postorder_node_iter():
            if n.is_leaf(): 
                node_likelihood[n] = dict()
                for site in range(numsites):
                    node_likelihood[n][site] = dict()
                    char_state = get_char(msa, n, site)
                    node_likelihood[n][site][char_state] = 1.0                
            elif n.is_internal(): 
                for site in range(numsites):
                    if n not in node_likelihood.keys():
                        node_likelihood[n] = dict()
                    node_likelihood[n][site] = dict()
                    node_likelihood = likelihood_under_n(node_likelihood, n, site, msa, Q, nwkt.seed_node is n)
        
        tree_likelihood = 0.0
        for site in range(numsites):
            site_likelihood = 0.0 
            for rootchar in node_likelihood[n][site].keys():
                prob_rootchar = node_likelihood[n][site][rootchar]
                if rootchar == 0:
                    if prob_rootchar == 0:
                        site_likelihood += 0
                    else:
                        site_likelihood += exp( log(math.exp(-root_edge_len)) + log(prob_rootchar) )
                    # GCLOG: site_likelihood += (math.exp(-root_edge_len)) * prob_rootchar # * q_ialpha 
                else:
                    q_ialpha = Q[site][rootchar]
                    if prob_rootchar == 0:
                        site_likelihood += 0
                    else:
                        site_likelihood += exp(log(1 - math.exp(-root_edge_len)) + log(q_ialpha) + log(prob_rootchar) )
                    # GCLOG: site_likelihood += ((1 - math.exp(-root_edge_len)) * q_ialpha * prob_rootchar)
            if site_likelihood > 0:
                tree_likelihood += log(site_likelihood)
            else:
                tree_likelihood += 0
        return -tree_likelihood

    if optimize_branchlengths: 
        x_star = []
        dmax = -log(1/numsites)*2
        dmin = -log(1-1/numsites)/2
        bound = (dmin, dmax)
        # print("bound", bound)

        x_star = None
        f_star = float("inf")

        x0 = []
        if init_tree: 
            init_tree = dendropy.Tree.get(data=init_tree, schema="newick", rooting="force-rooted")
            for i, e in enumerate(init_tree.postorder_edge_iter()): # visit the descendents before visiting edge
                x0.append(e.length)
            # print("initial likelihood", -felsenstein(x0))
            out = optimize.minimize(felsenstein, x0, method="SLSQP", options={'disp':False,'maxiter':1000}, bounds=[bound]*num_edges)
            x_star = out.x
            f_star = out.fun
        else: 
            all_results_x = []
            all_results_f = []
            for i in range(initials):
                x0 = [random() * (dmax - dmin) + dmin] * num_edges
                out = optimize.minimize(felsenstein, x0, method="SLSQP", options={'disp':False,'maxiter':1000}, bounds=[bound]*num_edges)
                all_results_x.append(out.x)
                all_results_f.append(out.fun)
                if out.success and out.fun < f_star:
                    x_star = out.x
                    f_star = out.fun
        if print_all:
            print(all_results_x)
            print(all_results_f)
        for i, e in enumerate(nwkt.postorder_edge_iter()):
            e.length = x_star[i]
        return -f_star, nwkt.as_string("newick"), x_star

    else:
        x0 = []
        # same way that we put it into the tree
        for i, e in enumerate(nwkt.postorder_edge_iter()):
            x0.append(e.length)

        return -felsenstein(x0), nwkt.as_string("newick"), x0

### Defining Max Parsimony Likelihood

def sankoff(T, C, msa, site):
    def get_state(n, site):
        if n.is_leaf():
            return msa[n.taxon.label][site]

    # post order traversal of T
    nodedict = dict()
    for n in T.postorder_node_iter():
        if n.is_leaf():
            nodedict[n] = dict()
            # print("n.taxon.label", n.taxon.label)
            state = get_state(n, site)
            nodedict[n][site] = dict()
            nodedict[n][site][state] = 0
        elif n.is_internal():
            cost = 0 
            possible_states = set()
            child_states = set()
            nodedict[n] = dict()

            # fill in possible states for parent
            for c in n.child_nodes():
                if c.is_leaf():
                    state = get_state(c, site)
                    child_states.add(state)
                else:
                    for state in nodedict[c].keys():
                        child_states.add(state)
            if len(child_states) == 1:
                if 0 not in child_states:
                    child_states.add(0)
                    possible_states = child_states
                else:
                    possible_states = {0}
   

            for c in n.child_nodes():
                min_state, min_val = None, float('inf')
                for s_p in possible_states:    

                    for s_c in nodedict[c].keys():
                        # check costs
                        if s_c == 0 and s_p != 0:
                            v = float('inf')
                        else: 
                            if s_p == s_c and s_p == 0:
                                v = nodedict[c][site][s_c]
                            else:
                                v = nodedict[c][site][s_c] + C[s_p][s_c]
                        if v < min_val:
                            min_val = v
                            min_state = s_c
                    nodedict[c][site] = dict()
                    nodedict[c][site][s_p] = min_val 
            

    return nodedict 

def num_zeros(l):
    return np.count_zeros(l)

def pars_likelihood(T, labels):
    # on each branch (i, j)
    # p_j = 1 - z_j, z_i where these are the numbers of zeros
    # log likelihood is sum of all log(p_j)
    ll = 0
    for e in T.postorder_edge_iter():
        i, j = e.rootedge, e.head_node
        a, b = labels[i], labels[j]
        # count number of zeros
        z_i, z_j = num_zeros(a), num_zeros(b)
        p_j = 1 - z_j, z_i
        ll += log(p_j)
    return ll


def mlpars(T, Q):
    # sankoff
    # calculate likelihood
    pass






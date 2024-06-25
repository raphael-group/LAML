import numpy as np
import dendropy
import math
from random import random,seed
from math import log,exp
from scipy import optimize
from scmail_libs.sequence_lib import read_sequences

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
            state = get_state(n, site)
            nodedict[n][site] = dict()
            nodedict[n][site][state] = 0
            print(n, nodedict)

        else:
            cost = 0 
            possible_states = set()
            child_states = set()
            nodedict[n] = dict()
            nodedict[n][site] = dict()

            # fill in possible states for parent
            for c in n.child_nodes():
                if c.is_leaf():
                    state = get_state(c, site)
                    child_states.add(state)
                else:
                    for state in nodedict[c].keys():
                        child_states.add(state)
            
            # print("child_states", child_states)
            if len(child_states) == 1 and 0 not in child_states:
                child_states.add(0)
                possible_states = child_states
            else:
                possible_states = {0}
            # print("possible_states", possible_states)

            for c in n.child_nodes():
                min_state, min_val = None, float('inf')

                for s_p in possible_states:    

                    for s_c in nodedict[c][site].keys():
                        # check costs
                        if s_c == 0 and s_p != 0:
                            v = float('inf')
                        else:
                            # print("s_p", s_p, "s_c", s_c)
                            if s_p == s_c and s_p == 0:
                                 # print(nodedict, c, site, s_c)
                                 v = nodedict[c][site][s_c]
                            else:
                                 v = nodedict[c][site][s_c] + C[s_p][s_c]
                        if v < min_val:
                            min_val = v
                            min_state = s_c
                    nodedict[n][site] = dict()
                    nodedict[n][site][s_p] = min_val 
            print(n, nodedict)
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






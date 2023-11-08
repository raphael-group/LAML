import numpy as np
import dendropy
import math
from treeswift import *
from math import *
from random import random
from sys import argv


# def prob_same(nodedict, node_likelihood, char, site, curr_node):
#     all_child_prob = 1.0
#     for c in curr_node.child_nodes():
#         char_state_prob = 0.0
#         for alpha in nodedict[c][site]:
#             tp = math.exp(-get_branchlen(c))
#             char_state_prob += tp * node_likelihood[c][site][char]
#         all_child_prob *= char_state_prob
#     return all_child_prob

def prob_same(node_likelihood, char, site, curr_node, use_log):
    # print("prob same")
    c = curr_node.child_nodes()[0]
    # if use_log:
    #     return -get_branchlen(c) + node_likelihood[c][site][char]
    # else:
    tp = math.exp(-get_branchlen(c))
    return tp * node_likelihood[c][site][char]

def prob_change(msa, q_dict, node_likelihood, site, curr_node, use_log):
    # print("prob change")
    # if use_log:
    #     all_child_prob = 0.0
    # else:
    all_child_prob = 1.0
    for c in curr_node.child_nodes():
        char = get_char(msa, c, site)
        if char != 0: 
            q_ialpha = q_dict[site][char] 
            # print("q_ialpha HERE", q_ialpha)
            # if use_log:
            #     tp = np.log(q_ialpha) + np.log(1 - math.exp(-get_branchlen(c)))    
            #     all_child_prob += tp + node_likelihood[c][site][char]
            # else:
            tp = q_ialpha * (1 - math.exp(-get_branchlen(c)))
            all_child_prob *= tp * node_likelihood[c][site][char]
            # print("tp:", tp)
            # print("nl", node_likelihood[c][site][char], "char", char)
        else:
            #char = get_char(msa, c, site)
            #print("other zero child", char, "should be 0")
            q_ialpha = q_dict[site][0]
            # if use_log:
            #   all_child_prob += -get_branchlen(c) + node_likelihood[c][site][0]
            # else:
            tp = math.exp(-get_branchlen(c))
            all_child_prob *= tp * node_likelihood[c][site][0]
    return all_child_prob

# def prob_change(q_dict, nodedict, node_likelihood, site, char, curr_node):
#     all_child_prob = 1.0
#     for c in curr_node.child_nodes():
#         char_state_prob = 0.0
#         for alpha in nodedict[c][site]:
#             q_ialpha = q_dict[site][alpha] 
#             tp = q_ialpha * (1 - math.exp(get_branchlen(c)))
#             char_state_prob += tp * node_likelihood[c][site][char]
#         all_child_prob *= char_state_prob
#     return all_child_prob


def likelihood_under_n(node_likelihood, n, site, msa, q_dict, use_log):
    # n is an internal node
    child_states = set()
        
    if n not in node_likelihood:
        node_likelihood[n] = dict()
        node_likelihood[n][site] = dict()
        
    # identify all child states. 
    # this constrains n's possible states.
    child_states = []
    for child in n.child_nodes():
        if child.is_leaf():
            child_states.append(get_char(msa, child, site))
        else:
            for x in node_likelihood[child][site]:
                state_prob = node_likelihood[child][site][x]
                # if not use_log and state_prob > 0.0:
                if state_prob > 0.0:
                    child_states.append(x)
                # elif use_log:
                #     child_states.append(x)
                    
    parent_poss_states = dict()

    print("Child states:", child_states)
    if 0 in child_states: # probability 0 -> 0
        if len(set(child_states)) == 1: # both children are state 0 
            # print("both children are 0")
            # if use_log:
            #     parent_poss_states[0] = prob_same(node_likelihood, 0, site, n, use_log) * 2 
            # else:
                # parent_poss_states[0] = prob_same(node_likelihood, 0, site, n, use_log) ** 2 
            parent_poss_states[0] = prob_same(node_likelihood, 0, site, n, use_log) ** 2
            # print("Probability of 0->0 is:", prob_same(node_likelihood, 0, site, n))
        else: # one child has state 0, other has non-zero state
            # print("one child has state 0, other is nonzero")
            # if not use_log:
            for c in child_states: # probability c -> c != 0
                parent_poss_states[c] = 0.0
            # probability 0 -> c (alpha)
            parent_poss_states[0] = prob_change(msa, q_dict, node_likelihood, site, n, use_log)  
            # print("Probability of 0->", child_states, "is:", prob_change(msa, q_dict, node_likelihood, site, n))
    else:
        if len(set(child_states)) == 1: # both children are same nonzero state
            # print("both children are alpha", site)
            c = child_states[0]
            # if use_log:
            #     parent_poss_states[c] = 0.0
            # else:
            parent_poss_states[c] = 1.0 
            parent_poss_states[0] = prob_change(msa, q_dict, node_likelihood, site, n, use_log)
            # print("Probability of 0->", c, "is:", prob_change(msa, q_dict, node_likelihood, site, n))
        else:
            # print("children are alpha, beta, both nonzero")
            # if use_log:
            #     parent_poss_states[0] = 0.0
            # else:
            parent_poss_states[0] = 1.0
            # print("Probability of ALPHA->(ALPHA, BETA) is: 1.0")
    # print("parent possible states", "site", site)
    # print(parent_poss_states)

    for x in parent_poss_states.keys():
        # nodedict[n][site][x] = parent_poss_states[x]
        node_likelihood[n][site][x] = parent_poss_states[x]
    
    # print("node likelihood")
    # print(node_likelihood)

    return node_likelihood

def get_branchlen(child_node):
    if child_node.edge_length is None:
        print(child_node.child_nodes())
    # print("Getting branch length of:", child_node.edge_length)
    return child_node.edge_length

def get_char(msa, leaf_node, site):
    return msa[int(leaf_node.taxon.__str__().replace("'", ""))-1][site]

# def felsenstein(T, Q, msa, ordering, root_edge_len=0.2, use_log=False):
def felsenstein(T, Q, msa, root_edge_len=0.2, use_log=False):
    # takes in tree with branch lengths as input
    # output a likelihood

    ## HOUSEKEEPING 
    numsites = len(msa[0])

    # q_dict = dict()
    # for site in range(numsites):
    #     q_dict[site] = dict()
    #     # get alphabet 
    #     for char in np.unique(msa.T[site]):
    #         q_dict[site][char] = Q[site][char]
    alphabet = dict()
    for site in range(numsites):
        alphabet[site] = Q[site].keys()

    print("q_dict", Q)

    nwkt = dendropy.Tree.get(data=T, schema="newick")
    print(nwkt)

    for n in nwkt.leaf_node_iter():
        print(n.taxon, ''.join([str(get_char(msa, n, s)) for s in range(numsites)]))

    # nodedict = dict()
    node_likelihood = dict()

    ## CALCULATE THE LIKELIHOOD
    for n in nwkt.postorder_node_iter():
        # print("node:", n)
        if n.taxon is not None: # must be a leaf node, set up 
            # nodedict[n] = dict()
            node_likelihood[n] = dict()
            for site in range(numsites):
                # char_state = get_char(msa, n, site)
                # nodedict[n][site] = dict()
                # nodedict[n][site][char_state] = 1.0
                node_likelihood[n][site] = dict()
                for char in alphabet[site]:
                    node_likelihood[n][site][char] = 0.0
                char_state = get_char(msa, n, site)
                node_likelihood[n][site][char_state] = 1.0
            
        elif n.taxon is None: # must be an internal node
            for site in range(numsites):
                # print("site:", site)
                if n not in node_likelihood:
                    # nodedict[n] = dict()
                    node_likelihood[n] = dict()
                node_likelihood[n][site] = dict()
                
                # nodedict[n][site] = dict()
                # if use_log: 
                #   for char in alphabet[site]:
                #       node_likelihood[n][site][char] = 1.0
                for char in alphabet[site]:
                    node_likelihood[n][site][char] = 0.0
                
                node_likelihood = likelihood_under_n(node_likelihood, n, site, msa, Q, use_log)
                # node_likelihood = likelihood_under_n(nodedict, node_likelihood, n, site, msa, q_dict)
    print(node_likelihood)
    # last n is the provided root node 
    # print("Calculating likelihood according to a root node r*")
    # SETTING UP r*, say r* -> r is dist 0.2
    tree_likelihood = 1.0
    # print("Root characters and probabilities:")
    for site in range(numsites):
        # under node_likelihood, calculate the prob    
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

def main():
    # t = '[&R] ((1:0.5,2:0.5):0.60,(3:0.5,4:0.5):0.60);'
    # t_A = '[&R] (1:1.10,3:1.10,2:1.10,4:1.10);' # star tree
    # t_B = '[&R] ((1:0.96,3:0.96):0.14,2:1.10,4:1.10);' # non-star tree
    # t_C = '[&R] (1:1.167,2:1.167,3:1.167:4:1.167);' # non-star tree
    t = '[&R] (1:0.96,3:0.96);'
    # increase k
    # k = 1000
    # mu = 0.025
    # nstate=9
    # # borrowed from https://github.com/uym2/scmail/blob/main/scripts/sim_test.py
    # tree = read_tree_newick(t_B)
    # tree.root.seq = '0'*k
    # msa = []

    # for node in tree.traverse_preorder():
    #     if node.is_root():
    #         continue
    #     p = 1-exp(-mu*node.edge_length)
    #     pnode = node.parent
    #     seq = ''
    #     for c in pnode.seq:
    #         if c != '0':
    #             seq += c
    #         else:
    #             r = random()
    #             if r < p: # change state
    #                 seq += str(int(random()*nstate) + 1) 
    #             else: # add '0'
    #                 seq += '0'  
    #     node.seq = seq
    #     msa.append(np.array(list(node.seq)))
    # msa = np.array(msa)
    # q = [[mu] * nstate] * k
    # # need to fix this 
    # q_ordering = [[i for i in range(nstate)] * k]

    # evolve down a tree -> look at uyen's code to do this
    msa = np.array([[1,0,2], # species 1
                    [1,0,0], # species 2
                    [0,1,2], # species 3
                    [0,1,0]] # species 4
                )
    
    # prob of transitioning 0 -> 0, 0 -> 1, 0 -> 2
    q = [{0: 0.2, 1:0.3, 2:0.5}, # site 1
         {0: 0.2, 1:0.3, 2:0.5}, # site 2
         {0: 0.2, 1:0.3, 2:0.5}, # site 3
        ]

    likelihood = felsenstein(t, q, msa, 0.14, use_log=False)
    print("Tree Likelihood:", likelihood)

if __name__ == "__main__":
    main()

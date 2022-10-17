import numpy as np
import dendropy
import math

import math

def prob_same(nodedict, node_likelihood, site, curr_node):
    all_child_prob = 1.0
    for c in curr_node.child_nodes():
        char_state_prob = 0.0
        for alpha in nodedict[c][site]:
            tp = math.exp(-get_branchlen(c))
            char_state_prob += tp * node_likelihood[c][site]
        all_child_prob *= char_state_prob
    return all_child_prob

def prob_change(q_dict, nodedict, node_likelihood, site, curr_state, curr_node):
    all_child_prob = 1.0
    for c in curr_node.child_nodes():
        char_state_prob = 0.0
        for alpha in nodedict[c][site]:
            q_ialpha = q_dict[site][alpha] 
            tp = q_ialpha * (1 - math.exp(get_branchlen(c)))
            char_state_prob += tp * node_likelihood[c][site]
        all_child_prob *= char_state_prob
    return all_child_prob


def likelihood_under_n(nodedict, node_likelihood, n, site, msa, q_dict):
    # n is an internal node
    child_states = set()
        
    if n not in nodedict:
        nodedict[n] = dict()
        nodedict[n][site] = dict()
        
    # identify all child states. 
    # this constrains n's possible states.
    child_states = set()
    for child in n.child_nodes():
        if child.is_leaf():
            child_states.add(get_char(msa, child, site))
        else:
            for x in nodedict[child][site]:
                state_prob = nodedict[child][site][x]
                if state_prob > 0.0:
                    child_states.add(x)
                    
    parent_poss_states = dict()
    if len(child_states) == 1:
        if 0 in child_states: # probability 0 -> 0
            parent_poss_states[0] = prob_same(nodedict, node_likelihood, site, n) 
        else:
            for c in child_states: # probability c -> c != 0
                parent_poss_states[c] = 1.0 
            # probability 0 -> c (alpha)
            parent_poss_states[0] = prob_change(q_dict, nodedict, node_likelihood, site, 0, n)  
    else:
        # probability 0 -> 1 and 0 -> 2 or
        # probability 0 -> 0 and 0 -> 1 WLOG
        parent_poss_states[0] = 1.0
    for x in parent_poss_states.keys():
        nodedict[n][site][x] = parent_poss_states[x]
        node_likelihood[n][site] *= parent_poss_states[x]
    
    return nodedict, node_likelihood

def get_branchlen(child_node):
    if child_node.edge_length is None:
        print(child_node.child_nodes())
    return child_node.edge_length

def get_char(msa, leaf_node, site):
    return msa[int(leaf_node.taxon.__str__().replace("'", ""))-1][site]

def felsenstein(T, Q, msa):
    # takes in tree with branch lengths as input
    # output a likelihood

    ## HOUSEKEEPING 
    numsites = len(msa[0])
        
    char_probs = dict()
    for site in range(numsites):
        char_probs[site] = dict()
        chars, counts = np.unique(msa.T[site], return_counts=True)
        for i in range(len(chars)):
            char_probs[site][chars[i]] = counts[i]/len(msa.T[site])

    q_dict = dict()
    for site in range(numsites):
        q_dict[site] = dict()
        # get alphabet 
        for idx, char in enumerate(np.unique(msa.T[site])):
            q_dict[site][char] = Q[site][idx]

    nwkt = dendropy.Tree.get(data=T, schema="newick")
    print(nwkt)

    nodedict = dict()
    node_likelihood = dict()

    ## CALCULATE THE LIKELIHOOD
    for n in nwkt.postorder_node_iter():
        # print("node:", n)
        if n.taxon is not None: # must be a leaf node, set up 
            nodedict[n] = dict()
            node_likelihood[n] = dict()
            for site in range(numsites):
                char_state = get_char(msa, n, site)
                nodedict[n][site] = dict()
                nodedict[n][site][char_state] = 1.0
                node_likelihood[n][site] = 1.0
            
        elif n.taxon is None: # must be an internal node
            for site in range(numsites):
                # print("site:", site)
                if n not in nodedict:
                    nodedict[n] = dict()
                    node_likelihood[n] = dict()
                
                nodedict[n][site] = dict()
                node_likelihood[n][site] = 1.0
                
                nodedict, node_likelihood = likelihood_under_n(nodedict, node_likelihood, n, site, msa, q_dict)

    # last n is the provided root node 
    # print("Calculating likelihood according to a root node r*")
    # SETTING UP r*, say r* -> r is dist 0.2
    root_edge_len = 0.2
    tree_likelihood = 1.0
    for site in range(numsites):
        # under node_likelihood, calculate the prob    
        for rootchar in nodedict[n][site].keys():
            prob_rootchar = nodedict[n][site][rootchar]
            if prob_rootchar > 0.0: 
                if rootchar == 0:
                    tree_likelihood *= (math.exp(-root_edge_len)) * node_likelihood[n][site]
                else:
                    q_ialpha = q_dict[site][rootchar]
                    tree_likelihood *= (1 - math.exp(-root_edge_len)) * q_ialpha * node_likelihood[n][site]
    
    print("Tree Likelihood:", tree_likelihood)
    return tree_likelihood

def main():
    t = '[&R] ((1:0.5,2:0.5):0.60,(3:0.5,4:0.5):0.60);'
    msa = np.array([[1,0,2], # species 1
                    [1,0,0], # species 2
                    [0,1,2], # species 3
                    [0,1,0]] # species 4
                )
    q = [[0.2, 0.3, 0.5], # site 1
         [0.2, 0.3, 0.5], # site 2
         [0.2, 0.3, 0.5]  # site 3
        ]
    likelihood = felsenstein(t, q, msa)

if __name__ == "__main__":
    main()
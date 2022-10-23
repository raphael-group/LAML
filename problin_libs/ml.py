import numpy as np
import dendropy
import math



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



def likelihood_under_n(nodedict, node_likelihood, n, site, msa, q_dict, use_log):
    child_states = set()
        
    if n not in nodedict:
        nodedict[n] = dict()
        nodedict[n][site] = dict()
        
    child_states = []
    for child in n.child_nodes():
        if child.is_leaf():
            child_states.append(get_char(msa, child, site))
        else:
            for x in nodedict[child][site]:
                state_prob = nodedict[child][site][x]
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

def felsenstein(T, Q, msa, ordering, root_edge_len=0.2, use_log=False):
    numsites = len(msa[0])
    def char_to_index(c, site, ordering):
        return ordering[site].index(c)

    q_dict = dict()
    for site in range(numsites):
        q_dict[site] = dict()
        for char in np.unique(msa.T[site]):
            q_dict[site][char] = Q[site][char_to_index(char, site, ordering)]
    alphabet = dict()
    for site in range(numsites):
        alphabet[site] = q_dict[site].keys()

    print("q_dict", q_dict)
    nwkt = dendropy.Tree.get(data=T, schema="newick")
    print(nwkt)

    for n in nwkt.leaf_node_iter():
        print(n.taxon, ''.join([str(get_char(msa, n, s)) for s in range(numsites)]))

    nodedict = dict()
    node_likelihood = dict()

    ## CALCULATE THE LIKELIHOOD
    for n in nwkt.postorder_node_iter():
        # print("node:", n)
        if n.taxon is not None: # must be a leaf node, set up 
            nodedict[n] = dict()
            node_likelihood[n] = dict()
            for site in range(numsites):
                node_likelihood[n][site] = dict()
                for char in alphabet[site]:
                    node_likelihood[n][site][char] = 0.0
                char_state = get_char(msa, n, site)
                node_likelihood[n][site][char_state] = 1.0
            
        elif n.taxon is None: # must be an internal node
            for site in range(numsites):
                if n not in nodedict:
                    node_likelihood[n] = dict()
                node_likelihood[n][site] = dict()
                
                for char in alphabet[site]:
                    node_likelihood[n][site][char] = 1.0
                
                node_likelihood = likelihood_under_n(nodedict, node_likelihood, n, site, msa, q_dict, use_log)

    tree_likelihood = 1.0
    for site in range(numsites):
        for rootchar in node_likelihood[n][site].keys():
            prob_rootchar = node_likelihood[n][site][rootchar]
            print(rootchar, prob_rootchar)
            if prob_rootchar > 0.0: 
                q_ialpha = q_dict[site][rootchar]
                if rootchar == 0:
                    if use_log:
                        tree_likelihood += (-root_edge_len) + np.log(q_ialpha) + np.log(prob_rootchar)
                    else:
                        tree_likelihood *= (math.exp(-root_edge_len)) * q_ialpha * prob_rootchar
                    
                else:
                    if use_log:
                        tree_likelihood += np.log((1 - math.exp(-root_edge_len))) + np.log(q_ialpha) + np.log(prob_rootchar)
                    else:
                        tree_likelihood *= ((1 - math.exp(-root_edge_len)) * q_ialpha * prob_rootchar)
    
    return tree_likelihood


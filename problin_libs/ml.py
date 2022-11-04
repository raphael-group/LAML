import numpy as np
import dendropy
import math
from random import random,seed
from math import log,exp
from scipy import optimize
from problin_libs.sequence_lib import read_sequences

def prob_same(node_likelihood, char, site, curr_node, use_log):
    # print("prob same")
    #c = curr_node.child_nodes()[0]
    #tp = math.exp(-get_branchlen(c))
    #return tp * node_likelihood[c][site][char]
    c1,c2 = curr_node.child_nodes()
    tp1,tp2 = exp(-get_branchlen(c1)),exp(-get_branchlen(c2))
    return tp1*node_likelihood[c1][site][char]*tp2*node_likelihood[c2][site][char]

def prob_change(msa, q_dict, node_likelihood, site, curr_node, child_states, use_log):
    all_child_prob = 0.0
    c1,c2 = curr_node.child_nodes()
    l1 = get_branchlen(c1)
    l2 = get_branchlen(c2)
    for char1 in node_likelihood[c1][site]:
        p1 = node_likelihood[c1][site][char1]
        p1 *= exp(-l1) if char1==0 else (1-exp(-l1))*q_dict[site][char1]
        for char2 in node_likelihood[c2][site]:
            p2 = node_likelihood[c2][site][char2]
            p2 *= exp(-l2) if char2==0 else (1-exp(-l2))*q_dict[site][char2]
            all_child_prob += p1*p2
    return all_child_prob        
    '''
    all_child_prob = 1.0
    for c in curr_node.child_nodes():
        # print("[prob_change]: call to get_char()", c)
        if c.is_leaf():
            char = get_char(msa, c, site)
            if char != 0: 
                # print(q_dict[site], site, char)
                q_ialpha = q_dict[site][char] 
                tp = q_ialpha * (1 - math.exp(-get_branchlen(c)))
                all_child_prob *= tp * node_likelihood[c][site][char]
            else:
                q_ialpha = q_dict[site][0]
                tp = math.exp(-get_branchlen(c))
                all_child_prob *= tp * node_likelihood[c][site][0]
        else:
            for char in node_likelihood[c][site].keys():
                # print("[prob_change]: internal child char", char)
                if char != 0: 
                    # print(q_dict[site], site, char)
                    q_ialpha = q_dict[site][char] 
                    if node_likelihood[c][site][char] > 0:
                        tp = q_ialpha * (1 - math.exp(-get_branchlen(c)))
                        # print("[prob_change]: tp", tp)
                        all_child_prob += tp * node_likelihood[c][site][char] # GC: Edited this
                else:
                    q_ialpha = q_dict[site][0]
                    if node_likelihood[c][site][0] > 0:
                        tp = math.exp(-get_branchlen(c))
                        all_child_prob += tp * node_likelihood[c][site][0] # GC: Edited this
                # print("[prob_change]:", all_child_prob)
    return all_child_prob
    '''
    
def likelihood_under_n(node_likelihood, n, site, msa, q_dict, is_root, use_log):
    # print("[likelihood_under_n]", n, "site", site)
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

    # print("[likelihood_under_n]: child_states", child_states)
    parent_poss_states = dict()
    if 0 in set(child_states): # probability 0 -> 0
        if len(set(child_states)) == 1: # both children are state 0 
            tmp = prob_same(node_likelihood, 0, site, n, use_log)
            #if is_root:
            parent_poss_states[0] = tmp
            #else:
            #    parent_poss_states[0] = tmp **2
        else: 
            for c in child_states: # probability c -> c != 0
                parent_poss_states[c] = 0.0
            # probability 0 -> c (alpha)
            parent_poss_states[0] = prob_change(msa, q_dict, node_likelihood, site, n, child_states, use_log)  
    else:
        if len(set(child_states)) == 1: # both children are same nonzero state
            c = child_states[0]
            parent_poss_states[c] = 1.0 
        parent_poss_states[0] = prob_change(msa, q_dict, node_likelihood, site, n, child_states, use_log)
        #else:
        #    parent_poss_states[0] = 1.0

    for x in parent_poss_states.keys():
        node_likelihood[n][site][x] = parent_poss_states[x]

    return node_likelihood

def get_branchlen(child_node):
    if child_node.edge_length is None:
        print(child_node.child_nodes())
    return child_node.edge_length

def get_char(msa, leaf_node, site):
    # print(msa)
    return msa[leaf_node.taxon.label][site]

'''
def felsenstein(T, Q, msa, use_log=False): #, root_edge_len=0.0):
    # print("MSA", msa)
    numsites = len(msa[next(iter(msa.keys()))])
    # numsites = len(msa.key[0])

    alphabet = dict()
    for site in range(numsites):
        alphabet[site] = Q[site].keys()

    # print("q_dict", Q)
    nwkt = dendropy.Tree.get(data=T, schema="newick", rooting="force-rooted")
    # print(nwkt)

    # for n in nwkt.leaf_node_iter():
    #     print(n.taxon, ''.join([str(get_char(msa, n, s)) for s in range(numsites)]))

    node_likelihood = dict()

    ## CALCULATE THE LIKELIHOOD
    for n in nwkt.postorder_node_iter():
        # print("node:", n)
        if n.is_leaf(): # n.taxon is not None: # must be a leaf node, set up 
            node_likelihood[n] = dict()
            for site in range(numsites):
                node_likelihood[n][site] = dict()
                for char in alphabet[site]:
                    node_likelihood[n][site][char] = 0.0
                # print("[felsenstein]: in site", site)
                char_state = get_char(msa, n, site)
                node_likelihood[n][site][char_state] = 1.0
            
        elif n.is_internal(): # n.taxon is None: # must be an internal node
            for site in range(numsites):
                if n not in node_likelihood.keys():
                    node_likelihood[n] = dict()
                node_likelihood[n][site] = dict()
                for char in alphabet[site]:
                    node_likelihood[n][site][char] = 0.0

                node_likelihood = likelihood_under_n(node_likelihood, n, site, msa, Q, nwkt.seed_node is n, use_log)
    # print(node_likelihood)
    if use_log:
        # print("Using log.")
        tree_likelihood = 0.0
    else:
        # print("NOT using log.")
        tree_likelihood = 1.0

    if nwkt.is_rooted:
        # print("Tree provided was rooted.")
        # print("Tree is rooted, here is likelihood at root.", node_likelihood[n])
        for site in range(numsites):
            if use_log:
                tree_likelihood += np.log(node_likelihood[n][site][0])
            else:
                tree_likelihood *= node_likelihood[n][site][0]
        
    elif not nwkt.is_rooted:
        # print("Tree provided was NOT rooted.")
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
'''

def wrapper_felsenstein(T, Q, msa, use_log=True, initials=20, optimize_branchlengths=False, init_tree=None):
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
        # print("branch length of seed node", get_branchlen(nwkt.seed_node))
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
                    for char in alphabet[site]:
                        node_likelihood[n][site][char] = 0.0
                    char_state = get_char(msa, n, site)
                    node_likelihood[n][site][char_state] = 1.0                
            elif n.is_internal(): 
                for site in range(numsites):
                    if n not in node_likelihood.keys():
                        node_likelihood[n] = dict()
                    node_likelihood[n][site] = dict()
                    for char in alphabet[site]:
                        node_likelihood[n][site][char] = 0.0
                    node_likelihood = likelihood_under_n(node_likelihood, n, site, msa, Q, nwkt.seed_node is n, use_log)
        
        #print(node_likelihood[n.child_nodes()[0]])
        #print(node_likelihood[n.child_nodes()[1]])
        #print(node_likelihood[n])
        if use_log:
            tree_likelihood = 0.0
        else:
            tree_likelihood = 1.0
        for site in range(numsites):
            site_likelihood = 0.0 if use_log else 1.0
            for rootchar in node_likelihood[n][site].keys():
                prob_rootchar = node_likelihood[n][site][rootchar]
                if rootchar == 0:
                    site_likelihood += (math.exp(-root_edge_len)) * prob_rootchar # * q_ialpha 
                else:
                    q_ialpha = Q[site][rootchar]
                    site_likelihood += ((1 - math.exp(-root_edge_len)) * q_ialpha * prob_rootchar)
            if use_log:
                tree_likelihood += log(site_likelihood)
            else:
                tree_likelihood *= site_likelihood        
        '''    
        # assuming that n is the seed_node (forced or otherwise)
        for site in range(numsites):
            for rootchar in node_likelihood[n][site].keys():
                prob_rootchar = node_likelihood[n][site][rootchar]
                if prob_rootchar > 0.0: 
                    q_ialpha = Q[site][rootchar]
                    if rootchar == 0:
                        if use_log:  # standard log is base e
                            tree_likelihood += (-root_edge_len) + np.log(prob_rootchar) # + np.log(q_ialpha) 
                        else:
                            tree_likelihood *= (math.exp(-root_edge_len)) * prob_rootchar # * q_ialpha 
                    else:
                        if use_log:
                            tree_likelihood += np.log((1 - math.exp(-root_edge_len))) + np.log(q_ialpha) + np.log(prob_rootchar)
                        else:
                            tree_likelihood *= ((1 - math.exp(-root_edge_len)) * q_ialpha * prob_rootchar)
        '''
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

            for i in range(initials):
                x0 = [random() * (dmax - dmin) + dmin] * num_edges
                out = optimize.minimize(felsenstein, x0, method="SLSQP", options={'disp':False,'maxiter':1000}, bounds=[bound]*num_edges)
                if out.success and out.fun < f_star:
                    x_star = out.x
                    f_star = out.fun
                # print(i, [out.fun, out.x], end='')      

        # print("optimal likelihood",-felsenstein(x_star),-f_star)
        for i, e in enumerate(nwkt.postorder_edge_iter()):
            e.length = x_star[i]
        return -f_star, nwkt.as_string("newick"), x_star

    else:
        x0 = []
        # same way that we put it into the tree
        for i, e in enumerate(nwkt.postorder_edge_iter()):
            x0.append(e.length)
        return -felsenstein(x0), nwkt.as_string("newick"), x0

def sets(seq_a, seq_b):
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
    return [s_0, s_1a, s_1b, s_2, s_3]

def get_taxon_name(node):
    return int(node.taxon.__str__().replace("'", ""))

def est_d(seq_a, seq_b):
    # print("est_d")
    # print(seq_a, seq_b)
    
    # leaf nodes a and b
    Z_a = seq_a.count(0)
    Z_b = seq_b.count(0)
    Z_ab = set()
    
    for i in range(len(seq_a)):
        a_i, b_i = seq_a[i], seq_b[i]
        if a_i == 0:
            Z_ab.add(i)
        if b_i == 0:
            Z_ab.add(i)
        if a_i != b_i and a_i != 0 and b_i != 0:
            Z_ab.add(i)
    Z_ab = len(Z_ab)
    
    # print("Z_a, Z_b, Z_ab", Z_a, Z_b, Z_ab)
    
    d_a = - np.log(Z_a/Z_ab)
    d_b = - np.log(Z_b/Z_ab)
    
    
    if d_a == -0.0:
        d_a = 0.0
    if d_b == -0.0:
        d_b = 0.0
    
    # print("d_a, d_b", d_a, d_b)
    return d_a, d_b


# def optimize_len(alphabet_size, k, a, b):
def optimize_len(k, a, b, x0):
    x_star = []
    num_iter = 20    
    dmax = -log(1/k)*2
    dmin = -log(1-1/k)/2
    bound = (dmin,dmax)

    def branch_likelihood(x): 
        s_0_len, s_1a_len, s_1b_len, s_2_len, s_3_len = len(s_0), len(s_1a), len(s_1b), len(s_2), len(s_3)
        
        d_a, d_b, d_r = x
        q = 0.2
        
        p1 = -(s_1b_len + s_0_len) * d_a + (s_1a_len + s_3_len) * np.log(1 - math.exp(-d_a))
        p2 = -(s_1a_len + s_0_len) * d_b + (s_1b_len + s_3_len) * np.log(1 - math.exp(-d_b)) - (k - s_2_len) * d_r
        p3 = 0.0
        
        for i in range(s_2_len): # assuming that prob to any alpha is the same
            # q_ialpha is the transition rate at site i from 0 -> character at node a at site i
            # iterate over sites in s_2, get q for each alpha
            p3 += np.log(q**2 * (1 - math.exp(-d_a)) * (1 - math.exp(-d_b)) * math.exp(-d_r) + q*(1 - math.exp(-d_r)))
        
        return -(p1 + p2 + p3)

    for i in range(num_iter):
        if i > 0:
            x0 = np.random.uniform(dmin, dmax, 3)
        s_0, s_1a, s_1b, s_2, s_3 = sets(a, b)
        # eps = 1e-10


        # out = optimize.minimize(likelihood, x0, method="L-BFGS-B", options={'disp': False}, bounds=[bound, bound, bound])
        out = optimize.minimize(branch_likelihood, x0, method="SLSQP", options={'disp': False}, bounds=[bound, bound, bound])

        x_star.append((out['fun'], out['x']))
        # w, v = LA.eig(out.hess_inv.todense())
        # print(w > 0.0)
        
    return x_star


def get_seq(seq_dict, leaf_node):
    return seq_dict[int(leaf_node.taxon.__str__().replace("'", ""))]
               
def get_idx(leaf_node):
    # print(leaf_node)
    return int(leaf_node.taxon.__str__().replace("'", ""))

def my_print_tree(t):
    for n in t.postorder_node_iter():
        n.label = ''
    print(t)

def run_trees(fname, topo, m, num_to_run=2):
    trees = []
    sets_of_four = []
    with open(fname) as r:
        four_leaves = []
        for record in SeqIO.parse(r, "fasta"):
            s = str(record.seq)
            if '_________' in record.seq:
                s = s.replace('_________', '')
                four_leaves.append([int(x) for x in s.split('|')])
                sets_of_four.append(four_leaves)
                four_leaves = []
            else:
                four_leaves.append([int(x) for x in s.split('|')])
    print(len(sets_of_four), "trees read.")

    all_dists = dict()
    all_seqs = dict()
    
    for idx, seqs in enumerate(sets_of_four):
        dists = dict()
        seq_dict = dict()
        num_seqs = len(seqs)

        if idx == num_to_run:
                break

        for i in range(num_seqs):
            seq_dict[i] = seqs[i]
            
        for pair in itertools.combinations(list(range(num_seqs)), 2):
            # if idx % 5 == 0:
                # print(idx)

            a, b = pair
            x0 = np.random.uniform(0, 5.0, 3)
            x_star = optimize_len(m, seqs[a], seqs[b], x0) 
            # print the stddev of x_star
            d_ab = sorted(x_star, key=lambda x: x[0], reverse=True)[0][1]
            d_a, d_b, d_r = d_ab
            dists[(a, b)] = [d_a, d_b, d_r]

        all_dists[idx] = dists
        all_seqs[idx] = seq_dict

        # estimate branch length for every pair of leaves
        for n in topo.postorder_node_iter():
            if n.is_internal(): # must be an internal node
                node_a, node_b = n.child_nodes()
                if node_a.is_leaf(): 
                    idx_a = get_idx(node_a)
                    seq_a = get_seq(seq_dict, node_a)
                if node_b.is_leaf():
                    idx_b = get_idx(node_b)
                    seq_b = get_seq(seq_dict, node_b)
                if node_a.is_leaf() and node_b.is_leaf():
                    if (idx_a, idx_b) in dists:
                        d_a, d_b, d_r = dists[(idx_a, idx_b)]
                        node_a.edge.length = d_a
                        node_b.edge.length = d_b
                        n.edge.length = d_r
                    elif (idx_b, idx_a) in dists:
                        d_b, d_a, d_r = dists[(idx_b, idx_a)]
                        node_a.edge.length = d_a
                        node_b.edge.length = d_b
                        n.edge.length = d_r
                    else:
                        print("dists", dists)
        trees.append(topo)
    return trees


def run_mats(fname, topo, m, num_to_run=2):
    trees = []
    sets_of_four = []
    with open(fname) as r:
        four_leaves = []
        for record in SeqIO.parse(r, "fasta"):
            s = str(record.seq)
            if '_________' in record.seq:
                s = s.replace('_________', '')
                four_leaves.append([int(x) for x in s.split('|')])
                sets_of_four.append(four_leaves)
                four_leaves = []
            else:
                four_leaves.append([int(x) for x in s.split('|')])
    print(len(sets_of_four), "trees read.")

    all_dists = dict()
    all_seqs = dict()
    
    for idx, seqs in enumerate(sets_of_four):
        dists = dict()
        seq_dict = dict()
        num_seqs = len(seqs)

        if idx == num_to_run:
                break

        for i in range(num_seqs):
            seq_dict[i] = seqs[i]
            
        for pair in itertools.combinations(list(range(num_seqs)), 2):

            a, b = pair
            x0 = np.random.uniform(0, 5.0, 3)
            x_star = optimize_len(m, seqs[a], seqs[b], x0) 
            # print the stddev of x_star
            d_ab = sorted(x_star, key=lambda x: x[0], reverse=True)[0][1]
            d_a, d_b, d_r = d_ab
            dists[(a, b)] = [d_a, d_b, d_r]

        all_dists[idx] = dists
        all_seqs[idx] = seq_dict

    return all_dists


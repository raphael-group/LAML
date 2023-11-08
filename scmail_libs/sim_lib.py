import os
from treeswift import *
from math import *
import random
from scmail_libs.sequence_lib import load_pickle
from random import lognormvariate, randint
#import networkx as nx
#import cassiopeia as cass
import pickle

"""
def simTree_lnorm(nLeaves,scale,std,randseed=None):
    # simulate tree using Cassiopeia where branch lengths
    # follow a lognormal distribution
    if randseed is None:
        randseed = randint(0,1024)        
    log_sd = sqrt(log((1+sqrt(1+4*(std*std)))/2))
    bd_sim = cass.sim.BirthDeathFitnessSimulator(
        birth_waiting_distribution = lambda scale:lognormvariate(log(scale),log_sd),
        initial_birth_scale = scale,
        num_extant = nLeaves,
        random_seed = randseed
    )
    ground_truth_tree = bd_sim.simulate_tree()
    nwstr = ground_truth_tree.get_newick(record_branch_lengths=True)
    return nwstr
"""

def get_balanced_tree(tree_height,branch_length,num_nodes=None):
# create a fully balanced tree with height = `tree_height`
# each branch length = `branch_length`
    root = Node("n0",branch_length)
    #root = Node("n0",branch_length)
    root.h = 0
    node_list = [root] # serve as a stack
    idx = 1
    while node_list:
        pnode = node_list.pop()
        h = pnode.h
        if h < tree_height:
            #cnode1 = Node(str(idx),branch_length)
            cnode1 = Node("n"+str(idx),branch_length)
            cnode1.h = h+1
            node_list.append(cnode1)
            pnode.add_child(cnode1)
            if num_nodes:
                if num_nodes < idx:
                    break
            idx += 1
            
            cnode2 = Node("n"+str(idx),branch_length)
            #cnode2 = Node(str(idx),branch_length)
            cnode2.h = h+1

            node_list.append(cnode2)
            pnode.add_child(cnode2)
            if num_nodes:
                if num_nodes < idx:
                    break
            
            idx += 1
    tree = Tree()
    tree.root = root
    return tree.newick() 

def simulate_seqs(tree,Q, mu=1.0, silencing_rate=0, dropout_rate=0, s=None ):
    if s is not None:
       random.seed(s)
    k = len(Q)
    z_seq = [0]*k
    for node in tree.traverse_preorder():
        d = node.edge_length * mu if node.edge_length is not None else 0
        #p = 1-exp(-d) # prob of nonzero mutation
        p = exp(-d) 
        p_seq = node.parent.seq if not node.is_root() else z_seq
        seq = []
        # simulate the child sequence
        for i, c in enumerate(p_seq):
            alphabet = list(Q[i].keys())
            m = len(alphabet)
            w = [0]*(m+2)
            nu = silencing_rate
            w[-1] = 1-p**nu if c != 's' else 1
            if c == 0:
                w[0] = p**(nu+1)
                for j,alpha in enumerate(alphabet):
                    w[j+1] = Q[i][alpha]*p**nu*(1-p)
            elif c != 's':
                j = alphabet.index(c)
                w[j+1] = p**nu
            nc = random.choices([0]+alphabet+['s'],weights=w)
            # simulate dropout
            if node.is_leaf() and nc[0] != 's' and random.random() < dropout_rate:
                nc = ['d']
            seq += nc 
        # set the sequence
        node.seq = seq    

    # full (internal included) character matrix    
    char_full = {}
    for node in tree.traverse_preorder():
        char_full[node.label] = node.seq
    # leaf (observed) character matrix    
    char_mtrx = {}
    for node in tree.traverse_leaves():
        char_mtrx[node.label] = [c if c not in ['s','d'] else '?' for c in node.seq]
    return char_mtrx,char_full         

def simulate_dropout(C, d):
    n_cmtx = dict()
    for nlabel in C:
        seq = C[nlabel]
        ns = [0]*len(seq) 
        for i,c in enumerate(seq):
            if c == -1:
                ns[i] = '?'
            else:    
                r = random.random()
                ns[i] = '?' if r<d else c
        n_cmtx[nlabel] = ns
    return n_cmtx

def sim_Q(k, m, prior_outfile=""):
    Q = []
    for i in range(k):
        q = {j+1:1/m for j in range(m)}
        Q.append(q)
    if prior_outfile != "":
        with open(prior_outfile, "w") as fout:
            for i in range(k):
                for x in Q[i]:
                    fout.write(str(i) + "," + str(x) + "," + str(Q[i][x]) + "\n") 
    return Q

def concat_Q(d):
    all_priors = []
    for f in os.listdir(d):
        #p = load_pickle(d+"/"+f)
        fstream = open(d+"/"+f,'rb')
        p = pickle.load(fstream)
        fstream.close()
        for k in p:
            Sum = sum(p[k].values())
            p[k] = {x:p[k][x]/Sum for x in p[k]}
            all_priors.append(p[k])
    return all_priors

def sample_Q(k, all_priors,s=None):
    if s is not None:
        random.seed(s)
    newQ = dict()
    for i in range(k):
        p = random.choice(all_priors)
        newQ[i] = p
    '''    
    if prior_outfile != "":
        with open(prior_outfile, "w") as fout:
            for i in range(k):
                for x in Q[i]:
                    fout.write(str(i) + " " + str(x) + " " + str(Q[i][x]) + "\n") '''
    return newQ

if __name__=="__main__":
    treeStr = get_balanced_tree(2,1.0)
    tree = read_tree_newick(treeStr)
    k=30
    m=5
    Q = []
    for i in range(k):
        q = {j+1:1/m for j in range(m)}
        Q.append(q)
    char_mtrx = simulate_seqs(tree,Q)
    print(treeStr)
    print(char_mtrx)

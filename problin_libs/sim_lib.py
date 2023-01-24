from treeswift import *
from math import *
import random

def get_balanced_tree(tree_height,branch_length):
# create a fully balanced tree with height = `tree_height`
# each branch length = `branch_length`
    root = Node("n0",branch_length)
    root.h = 0
    node_list = [root] # serve as a stack
    idx = 1
    while node_list:
        pnode = node_list.pop()
        h = pnode.h
        if h < tree_height:
            cnode1 = Node("n"+str(idx),branch_length)
            cnode1.h = h+1
            cnode2 = Node("n"+str(idx+1),branch_length)
            cnode2.h = h+1
            pnode.add_child(cnode1)
            pnode.add_child(cnode2)
            node_list.append(cnode1)
            node_list.append(cnode2)
            idx += 2
    tree = Tree()
    tree.root = root
    return tree.newick() 

def simulate_seqs(tree,Q, mu=1.0, with_heritable=False, silencing_rate=1e-4, s=1984):
    random.seed(s)
    k = len(Q)
    z_seq = [0]*k
    for node in tree.traverse_preorder():
        d = node.edge_length * mu if node.edge_length is not None else 0
        p = 1-exp(-d)
        p_seq = node.parent.seq if not node.is_root() else z_seq
        seq = []
        for i,c in enumerate(p_seq):
            if c != 0:
                nc = [c]
                r = random.random()
                if r < p and with_heritable and c != -1:
                    r = random.random()
                    if r < silencing_rate:
                        nc = [-1]
                seq += nc
            else:
                r = random.random()
                if r < p: # change state
                    seq += random.choices(list(Q[i].keys()),weights=list(Q[i].values()),k=1)
                else: # keep 0
                    seq += [0]  
        node.seq = seq
    char_mtrx = {}
    for node in tree.traverse_leaves():
        char_mtrx[node.label] = node.seq
    return char_mtrx          

def simulate_dropout(C, d, s):
    random.seed(s)
    n_cmtx = dict()
    for nlabel in C:
        seq = C[nlabel]
        ns = [] 
        for c in seq:
            r = random.random()
            if r < d:
                ns += [-1]
            else:
                ns += [c]
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
                    fout.write(str(i) + " " + str(x) + " " + str(Q[i][x]) + "\n")
    return Q

def setup(tree, m, mu, k):
    Q = sim_Q(k, m)
    cmtx = simulate_seqs(tree, Q, mu)
    return cmtx

def count_missing(mtx):
    nzeros, nmissing, total = 0, 0, 0
    for c in mtx:
        seq = mtx[c]
        nzeros += sum([1 if ch == 0 else 0 for ch in seq])
        nmissing += sum([1 if ch == -1 else 0 for ch in seq])
        total += len(seq)

    # print(nzeros, nmissing, total)
    return nzeros, nmissing, total

def count_all(repdict):
    all_zero, all_missing, all_total = 0, 0, 0
    for i in repdict:
        mtx = repdict[i]
        nz, nm, total = count_missing(mtx)
        all_zero += nz
        all_missing += nm
        all_total += total
    return all_zero, all_missing, all_total
    
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

from treeswift import *
from math import *
from random import random, choices

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

def simulate_seqs(tree,Q):
    k = len(Q)
    z_seq = [0]*k
    for node in tree.traverse_preorder():
        d = node.edge_length if node.edge_length is not None else 0
        p = 1-exp(-d)
        p_seq = node.parent.seq if not node.is_root() else z_seq
        seq = []
        for i,c in enumerate(p_seq):
            if c != 0:
                seq += [c]
            else:
                r = random()
                if r < p: # change state
                    seq += choices(list(Q[i].keys()),weights=list(Q[i].values()),k=1)
                else: # keep 0
                    seq += [0]  
        node.seq = seq
    char_mtrx = {}
    for node in tree.traverse_leaves():
        char_mtrx[node.label] = node.seq
    return char_mtrx          

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

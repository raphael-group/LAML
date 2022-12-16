#! /usr/bin/env python
from math import *

def MLP_label(T,seqs):
    k = len(list(seqs.values())[0])
    d_max = -log(1/k)*2
    d_min = -log(1-1/k)/2
    #d_min = 0.005025167926750725
    #d_max = 9.210340371976182 
    for node in T.traverse_postorder():
        if node.is_leaf():
            node.seq = seqs[node.label]
            node.z = len([x for x in node.seq if x==0])
        else:
            a,b = node.children
            s1 = a.seq
            s2 = b.seq
            s = [x if x==y and x!=0 else 0 for x,y in zip(s1,s2)]     
            node.seq = s
            node.z = len([x for x in s if x==0])
            if node.z == 0:
                a.edge_length = b.edge_length = d_min # unidentifiable branch lengths
                continue
            p_a = 1-a.z/node.z
            if p_a == 0:
                d_a = d_min
            elif p_a == 1:
                d_a = d_max
            else:
                d_a = -log(1-p_a) 
            p_b = 1-b.z/node.z
            if p_b == 0:
                d_b = d_min
            elif p_b == 1:
                d_b = d_max
            else:
                d_b = -log(1-p_b) 
            a.edge_length = d_a
            b.edge_length = d_b
            #a.edge_length = p_a
            #b.edge_length = p_b
    # special care for the root
    p_0 = 1-T.root.z/k
    if p_0 == 0:
        d_0 = d_min
    elif p_0 == 1:
        d_0 = d_max
    else:
        d_0 = -log(1-p_0) 
    T.root.edge_length = d_0
    #T.root.edge_length = p_0

if __name__ == "__main__":
    from sys import argv
    from problin_libs.sequence_lib import *
    from treeswift import *

    tree = read_tree_newick(argv[1])
    char_mtrx = read_sequences(argv[2],filetype="charMtrx") 
    MLP_label(tree,char_mtrx)
    tree.root.h = 1 # the root node has a branch above
    for node in tree.traverse_preorder():
        if not node.is_root():
            node.h = node.parent.h + 1
        if node.edge_length is not None:
            print(node.label + " " + str(node.h) + " " + str(node.edge_length))

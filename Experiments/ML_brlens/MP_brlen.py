#! /usr/bin/env python

def MP_label(T,seqs):
    MP_cost = 0
    for node in T.traverse_postorder():
        if node.is_leaf():
            node.seq = seqs[node.label]
        else:
            a,b = node.children
            s1 = a.seq
            s2 = b.seq
            s = [x if x==y and x!=0 else 0 for x,y in zip(s1,s2)]     
            node.seq = s
            a.edge_length = sum([int(x!=z) for (x,z) in zip(s1,s)])
            b.edge_length = sum([int(x!=z) for (x,z) in zip(s2,s)])
            MP_cost += (a.edge_length + b.edge_length)
    # special care for the root
    T.root.edge_length = len([x for x in T.root.seq if x != 0])
    MP_cost += T.root.edge_length
    return MP_cost

if __name__ == "__main__":
    from sys import argv
    from problin_libs.sequence_lib import *
    from treeswift import *

    tree = read_tree_newick(argv[1])
    char_mtrx = read_sequences(argv[2],filetype="charMtrx") 
    MP_label(tree,char_mtrx)
    tree.root.h = 1 # the root node has a branch above
    for node in tree.traverse_preorder():
        if not node.is_root():
            node.h = node.parent.h + 1
        print(node.label + " " + str(node.h) + " " + str(node.edge_length))

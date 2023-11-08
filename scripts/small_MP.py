#! /usr/bin/env python
from treeswift import *
from scmail_libs.sequence_lib import * 
from math import *

def small_MP(nwkStr,seqs,Q=None):
    T = read_tree_newick(nwkStr)
    if len(T.root.children) > 1:
        new_root = Node()
        new_root.add_child(T.root)
        T.root = new_root 
    MP_cost = 0
    all_labels = set()
    idx = 0
    for node in T.traverse_postorder():
        if node.is_leaf():
            k = len(seqs[node.label])
            node.seq = seqs[node.label]
        elif node.is_root():
            node.seq = [0]*k    
        else:
            C = node.children
            s0 = C[0].seq
            for cnode in C[1:]:
                s1 = cnode.seq
                s = []
                for x,y in zip(s1,s0):
                    if x == y or y == '?':
                        z = x
                    elif x == '?':
                        z = y
                    else:
                        z = 0
                    s.append(z)
                s0 = s
            node.seq = s0
    lb2seqs = {}
    for node in T.traverse_preorder():
        if node.label is None or node.label in all_labels:
            node.label = "I" + str(idx)
            idx += 1
        lb2seqs[node.label] = node.seq        
        for cnode in node.children:
            cnode.edge_length = 0
            for i,(x,z) in enumerate(zip(node.seq,cnode.seq)):
                if z == '?':
                    z = x
                    cnode.seq[i] = z
                if z != 0:    
                    cost = (x != z)*(-log(Q[i][z]) if Q is not None else 1)    
                    cnode.edge_length += cost
            MP_cost += cnode.edge_length
            
    return T.newick(),lb2seqs,MP_cost

if __name__ == "__main__":
    from sys import argv
    #msa_file = "s0d100p01_character_matrix.csv"
    msa_file = argv[1]
    tree_file = argv[2]
    out_file = argv[3]

    msa, site_names = read_sequences(msa_file,filetype="charMtrx",delimiter=",",masked_symbol='?')
    #prior_file = "prior_k30_r01.csv"
    #Q = read_priors(prior_file, site_names)
    #tree_file = "s0d100p01_tree.nwk"
    with open(tree_file,'r') as fin:
        treeStr = fin.read().strip()
        
    annTree,lb2seqs,MP_cost = small_MP(treeStr,msa,Q=None)
    #with open("s0d100p01_MP.txt",'w') as fout:
    with open(out_file,'w') as fout:
        fout.write(annTree + "\n")
        for lb in lb2seqs:
            fout.write(lb + "," + ",".join([str(x) for x in lb2seqs[lb]]) + "\n")
    print("MP cost:" + str(MP_cost))        

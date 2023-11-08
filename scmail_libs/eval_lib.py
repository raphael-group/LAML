#! /usr/bin/env python
from scmail_libs.sequence_lib import read_sequences,read_charMtrx
from scmail_libs.lca_lib import find_LCAs
from scmail_libs.EM_solver import EM_solver
from math import *
from treeswift import *

def score_char(trueChar,estCharProbs,silencing_symbol='-1'):
    if trueChar == 'd':
        score = 1
        if '-1' in estCharProbs:
            score -= estCharProbs['-1']
    else:
        if trueChar == 's':
            trueChar = silencing_symbol    
        score = 0
        if trueChar in estCharProbs:
            score = estCharProbs[trueChar]
    return score    

def get_charProbs(tokens,soft_assignment=True,masked_symbol='?'):
    estCharProbs = {}
    if len(tokens) == 1:
        if tokens[0] == masked_symbol:
            return None
        estCharProbs[tokens[0]] = 1
    else:
        pmax = 0
        best_c = None
        for token in tokens:
            c,p = token.split(":")   
            p = float(p) 
            if p > pmax:
                pmax = p
                best_c = c
            estCharProbs[c] = float(p)
        if not soft_assignment:
            estCharProbs = {}
            estCharProbs[best_c] = 1
    return estCharProbs

def score_seq(trueCharList,estCharList,soft_assignment=True,silencing_symbol='-1',masked_symbol='?'):
    score = 0
    k = 0
    for (x,y) in zip(trueCharList,estCharList):
        tokens = y.split("/")
        estCharProbs = get_charProbs(tokens,soft_assignment=soft_assignment,masked_symbol=masked_symbol)   
        if estCharProbs is not None:
            score += score_char(x,estCharProbs,silencing_symbol=silencing_symbol)        
            k += 1
    return score/k

def allelic_coupling(char_mtrx,cells,masked_symbol='?'):
    N = len(cells)
    AC_answer = {} 
    for i in range(N-1):
        for j in range(i+1,N):
            seq_i = char_mtrx[cells[i]]
            seq_j = char_mtrx[cells[j]]
            seq_ij = []
            d_ij = 0
            for (x,y) in zip(seq_i,seq_j):
                if x == y or y == masked_symbol:
                    z = x
                elif x == masked_symbol:
                    z = y
                else:
                    z = '0'
                seq_ij.append(z)
                if x != masked_symbol and str(x) != str(z):
                    d_ij += 1
                if y != masked_symbol and str(y) != str(z):
                    d_ij += 1
            if cells[i] < cells[j]:
                c1,c2 = cells[i],cells[j]
            else:    
                c1,c2 = cells[j],cells[i]
            AC_answer[(c1,c2)] = (seq_ij,d_ij)
    return AC_answer            

def tree_coupling(tree,cells,charMtrx):
    # assume that all nodes in tree are present in the charMtrx
    selected_leaves = []
    cells_set = set(cells)
    for node in tree.traverse_leaves():
        if node.label in cells_set:    
            selected_leaves.append(node)
    N = len(selected_leaves)
    answer = {}
    D = tree.distance_matrix()
    selected_pairs = []
    for i in range(N-1):
        for j in range(i+1,N):
            selected_pairs.append((selected_leaves[i],selected_leaves[j]))
    selected_lcas = find_LCAs(tree,[(x.label,y.label) for (x,y) in selected_pairs])

    for (leaf_i,leaf_j),lca_ij in zip(selected_pairs,selected_lcas):
        d_ij = D[leaf_i][leaf_j]
        s_ij = charMtrx[lca_ij.label]
        if leaf_i.label < leaf_j.label:
            c1,c2 = leaf_i.label,leaf_j.label
        else:    
            c1,c2 = leaf_j.label,leaf_i.label
        answer[(c1,c2)] = (s_ij,d_ij)
    return answer            

def read_annotation(infile):
    fin = open(infile,'r')
    charMtrx,_ = read_charMtrx(fin,replace_mchar=None,convert_to_int=False)
    fin.close()

    fin = open(infile,'r')
    tree = read_tree_newick(fin.readline())
    fin.close()
   
    # place sequences onto the tree
    for node in tree.traverse_postorder():
        node.seq = charMtrx[node.label]
    
    return charMtrx,tree

def read_groundtruth(infile,restricted_leafset=None):
    fin = open(infile,'r')
    charMtrx,_ = read_charMtrx(fin,replace_mchar=None,convert_to_int=False,stop_key="Evolve tree")
    k = len(charMtrx['I0'])
    fin.close()

    fin = open(infile,'r')
    for line in fin:
        if line.startswith("Evolve tree:"):
            treeStr = line.split("Evolve tree:")[1]
            tree = read_tree_newick(treeStr)
    fin.close()

    if restricted_leafset is not None:
        tree = tree.extract_tree_with(restricted_leafset)

    # HACKING: the groundtruth does not give the actual state for every dropout,
    # but only puts 'd'. Below is an attempt to recover it
    for node in tree.traverse_preorder():
        node.seq = [x for x in charMtrx[node.label]]
        if node.is_leaf():
            p0 = exp(-node.edge_length)
            for i in range(len(node.seq)):
                if node.seq[i] =='d':
                    node.seq[i] = '0' if (node.parent.seq[i] == '0' and p0>0.5) else node.parent.seq[i]

    # convert edge lengths to the number of mutations                    
    for node in tree.traverse_preorder():
        for c in node.children:
                d = k-score_seq(c.seq,node.seq,soft_assignment=False,silencing_symbol='s')*k
                c.edge_length = d 
    return charMtrx,tree

def convert_to_hard_assignment(seq):
    new_seq = []
    for x in seq:
        tokens = x.split("/")
        charProbs = get_charProbs(tokens,soft_assignment=False,masked_symbol='?')
        y = list(charProbs.keys())[0]
        new_seq.append(y)
    return new_seq    
        
def count_mutation(ann_tree,silence_symbol='-1'):
    # count the total number of mutations on the tree for each site
    # ann_tree is a treeswift object; here assume the tree has been "annotated",
    # so that every node has attribute node.seq
    k = len(ann_tree.root.seq)
    N_mus = [0]*k
    for node in ann_tree.traverse_postorder():
        if node.is_root():
            continue
        seq = convert_to_hard_assignment(node.seq)
        pseq = convert_to_hard_assignment(node.parent.seq)
        for i,(x,y) in enumerate(zip(pseq,seq)):
            N_mus[i] += (x!=y and y!=silence_symbol)
    return N_mus 

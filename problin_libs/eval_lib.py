#! /usr/bin/env python
from problin_libs.sequence_lib import read_sequences,read_charMtrx

def score_char(trueChar,estCharProbs):
    if trueChar == 'd':
        score = 1
        if '-1' in estCharProbs:
            score -= estCharProbs['-1']
    else:
        if trueChar == 's':
            trueChar = '-1'    
        score = 0
        if trueChar in estCharProbs:
            score = estCharProbs[trueChar]
    return score    

def get_charProbs(tokens,soft_assignment=True):
    estCharProbs = {}
    if len(tokens) == 1:
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

def score_seq(trueCharList,estCharList,soft_assignment=True):
    score = 0
    for (x,y) in zip(trueCharList,estCharList):
        tokens = y.split("/")
        estCharProbs = get_charProbs(tokens,soft_assignment=soft_assignment)   
        score += score_char(x,estCharProbs)        
    return score

def allelic_coupling(char_mtrx,cells):
    N = len(cells)
    AC_answer = {}
    for i in range(N-1):
        for j in range(i+1,N):
            seq_i = char_mtrx[cells[i]]
            seq_j = char_mtrx[cells[j]]
            seq_ij = []
            d_ij = 0
            for (x,y) in zip(seq_i,seq_j):
                if x == y or y == '?':
                    z = x
                elif x == '?':
                    z = y
                else:
                    z = 0
                seq_ij.append(z)
                if x != '?' and str(x) != str(z):
                    d_ij += 1
                if y != '?' and str(y) != str(z):
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
    for i in range(N-1):
        for j in range(i+1,N):
            leaf_i = selected_leaves[i]
            leaf_j = selected_leaves[j]
            #d_ij = tree.distance_between(leaf_i,leaf_j)
            d_ij = D[leaf_i][leaf_j]
            #lca_ij = ...
            lca_ij = leaf_i
            s_ij = charMtrx[lca_ij.label]
            if leaf_i.label < leaf_j.label:
                c1,c2 = leaf_i.label,leaf_j.label
            else:    
                c1,c2 = leaf_j.label,leaf_i.label
            answer[(c1,c2)] = (s_ij,d_ij)
    return answer            

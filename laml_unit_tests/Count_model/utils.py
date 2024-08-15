from random import random
from laml_libs.Count_model.PMM_base import PMM_model, Alphabet,AlleleTable

def countgen(alphabet,chosen_state,silencing=False,maxcount=1000):
    # generates the UMI counts table for a cell, providing a mapping from cassette state to count
    M = len(alphabet)
    counts = [0]*M
    for i in range(M):
        counts[i] = int(random()*maxcount)
    m = max(counts) + int(maxcount/M)
    C = {}    
    if silencing:
        C[tuple([-1]*len(alphabet[0]))] = 0
    for i,a in enumerate(alphabet):    
        # indicates dropout
        if chosen_state == ('?',):
            C[a] = 0
        else:
            C[a] = counts[i]
            if a == chosen_state:
                C[a] = m
    return C       

def charMtrx_2_alleleTable(charMtrx,alphabet):
    K = alphabet.K
    J = alphabet.J
    data_struct = {}
    for cell_name in charMtrx:
        counts = [{}]*K
        for k in range(K):
            counts[k] = countgen(alphabet.get_cassette_alphabet(k),tuple([charMtrx[cell_name][k]]))
        data_struct[cell_name] = counts
    allele_table = AlleleTable(K,J,data_struct,alphabet)
    return allele_table

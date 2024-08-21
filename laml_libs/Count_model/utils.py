from math import *
from .Alphabet import Alphabet 
from .AlleleTable import AlleleTable
from random import random

def log_sum_exp(numlist):
    # using log-trick to compute log(sum(exp(x) for x in numlist))
    # mitigate the problem of underflow
    maxx = max(numlist)
    result = maxx + log(sum([exp(x-maxx) for x in numlist]))
    return result

def pseudo_log(x):
    return log(x) if x>0 else min_llh

def countgen(alphabet,chosen_state,silencing=False,maxcount=1000):
    # generates the UMI counts table for a cell, providing a mapping from cassette state to count
    M = len(alphabet)
    J = len(chosen_state)
    missing_state = tuple(['?']*J)
    counts = [0]*M
    for i in range(M):
        counts[i] = int(random()*maxcount)
    m = max(counts) + int(maxcount/M)
    C = {}    
    if silencing:
        C[tuple([-1]*len(alphabet[0]))] = 0
    for i,a in enumerate(alphabet):    
        # indicates dropout
        if chosen_state == missing_state:
            C[a] = 0
        else:
            C[a] = counts[i]
            if a == chosen_state:
                C[a] = m
    return C       

def charMtrx_2_alleleTable(charMtrx):
    # charMtrx is an instance of CharMtrx
    # alphabet is an instance of Alphabet
    # charMtrx and alphabet must have the same K and J
    alphabet = charMtrx.alphabet
    K = charMtrx.K
    J = charMtrx.J
    data_struct = {}
    for cell_name in charMtrx.get_cell_names():
        counts = [{}]*K
        for k in range(K):
            chosen_state = charMtrx.get(cell_name,k)
            counts[k] = countgen(alphabet.get_cassette_alphabet(k),chosen_state)
        data_struct[cell_name] = counts
    allele_table = AlleleTable(K,J,data_struct,alphabet)
    return allele_table

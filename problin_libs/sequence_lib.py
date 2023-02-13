#! /usr/bin/env python
from statistics import mean
import pandas as pd

recognized_missing = set(['-', '?', '-1'])

def write_sequences(char_mtrx,nsites,outFile,delimiter=","):
    with open(outFile,'w') as fout:
        # header
        fout.write("cell")
        for i in range(nsites):
            fout.write(delimiter+"r" + str(i+1))
        fout.write("\n")
        # actual data
        for cell in char_mtrx:
            fout.write(cell)
            for x in char_mtrx[cell]:
                fout.write(delimiter+str(x))
            fout.write("\n")


def read_sequences(inFile,filetype="charMtrx",delimiter=",",masked_symbol=None):
    with open(inFile,'r') as fin:
        if filetype == "fasta":
            print("Warning: Reading " + str(inFile) + " as fasta file. Processing missing data in these files is not yet implemented.")
            return read_fasta(fin)
        elif filetype == "charMtrx":

            return read_charMtrx(fin,delimiter=delimiter,masked_symbol=masked_symbol)

def read_fasta(fin):    
    S = [] # will be a list of dictionaries
    D = {}
    for line in fin:
        if line.startswith(">"):
            name = line.strip()[1:]
        elif line.startswith("_"):
            S.append(D)
            D = {}
        else:
            seq = [int(x) for x in line.strip().split("|")]
            D[name] = seq       
    return S

def check_missing(seen_missing, x):
    # returns whether character x is a missing character
    if x in seen_missing:
        return True
    elif x in recognized_missing:
        return True
    elif x.isalpha(): # check alphanumeric 
        return True 
    elif int(x) < 0: # check positivity
        return True 
    else:
        return False

def read_charMtrx(fin,delimiter=",",masked_symbol=None):
    D = {}
    site_names = fin.readline().strip().split(delimiter)[1:]
    

    if masked_symbol != '-':
        seen_missing = set([masked_symbol])
    else: 
        seen_missing = set([])

    for line in fin:
        line_split = line.strip().split(delimiter)
        name = line_split[0]
        # check if any are missing characters or nonnegative
        seq = []
        for x in line_split[1:]:
            if check_missing(seen_missing, x):
                seen_missing.add(x)
                seq.append('?')
            else:
                seq.append(int(x))
        #seq = [int(x) if x != masked_symbol else "?" for x in line_split[1:]]
        D[name] = seq
    if len(seen_missing) > 1:
        print("Warning: Found " + str(seen_missing) + " characters and treated them as missing.")
    elif masked_symbol == None:
        print("Warning: Reading sequences, detected " + str(seen_missing) + " as the missing character(s). We recommend explicitly providing the missing character.")
    return D, site_names    

def read_Q(inFile):
    with open(inFile,'r') as fin:
        fin.readline() # skip the header
        Q = {}
        for line in fin:
            char,state,prob = line.strip().split(',')
            if not char in Q:
                Q[char] = {int(state):float(prob)}
            else:
                Q[char][int(state)] = float(prob)
    return [Q[q] for q in sorted(Q)]

#from treeswift import *

def extract_brlens(tfile, ofile):

    t = read_newick_tree(tfile)
    t.root.h = 1

    with open(ofile, "w+") as w:
        for nidx, node in enumerate(t.traverse_preorder()):
            if not node.is_root():
                node.h = node.parent.h + 1
                s = nidx + " " + str(node.h) + " " + str(node.edge_length)
                w.write(s)

def alphabet_size(mtx):
    df = pd.DataFrame.from_dict(mtx, orient='index')
    unique_series = df.nunique()
    return max(unique_series), min(unique_series), mean(unique_series)


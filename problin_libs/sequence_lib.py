#! /usr/bin/env python

def read_sequences(inFile,filetype="fasta"):
    with open(inFile,'r') as fin:
        if filetype == "fasta":
            return read_fasta(fin)
        elif filetype == "charMtrx":
            return read_charMtrx(fin)

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

def read_charMtrx(fin):    
    D = {}
    fin.readline() # skip the header
    for line in fin:
        line_split = line.strip().split(",")
        name = line_split[0]
        seq = [int(x) for x in line_split[1:]]
        D[name] = seq
    return D    

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

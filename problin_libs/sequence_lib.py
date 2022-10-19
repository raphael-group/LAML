#! /usr/bin/env python

def read_sequences(inFile):
    S = [] # will be a list of dictionaries
    D = {}
    with open(inFile,'r') as fin:
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

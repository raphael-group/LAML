#! /usr/bin/env python

seq_count = {}

with open("3724_NT_T1_character_matrix.txt","r") as fin:
    fin.readline() # skip header
    for line in fin:
        split = line.strip().split()
        cellname = split[0]
        seq = "|".join(split[1:])
        if seq in seq_count:
            seq_count[seq] += 1
        else:
            seq_count[seq] = 1

for seq in seq_count:
    print(seq,seq_count[seq])                

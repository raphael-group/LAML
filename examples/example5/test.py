#! /usr/bin/env python
import re

infile = "TLS_097_unfiltered_cm.txt"
alleles = set([])
a1 = set([])
a2 = set([])
a3 = set([])

with open(infile,'r') as fin:
    fin.readline() # ignore the header
    for line in fin:
        tokens = line.strip().split(",{")
        a = set(re.findall(r"\(.*?\)",tokens[1][:-1]))
        for x in a:
            x1,x2,x3 = x[1:-1].split(",")
            a1.add(x1)
            a2.add(x2)
            a3.add(x3)
        alleles = alleles.union(a)
        #allele = tokens[1]#[1:-1]#.split(",")[0]
        #print(allele)
#print(sorted(a1))
print(len(alleles),len(a1),len(a2),len(a3))

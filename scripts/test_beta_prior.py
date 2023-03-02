#! /usr/bin/env python

from problin_libs.ML_solver import beta_prior
from problin_libs.sequence_lib import read_sequences
from sys import argv

msa,_ = read_sequences(argv[1],masked_symbol="-1",filetype="charMtrx")
total = 0
missing = 0

for x in msa:
    total += len(msa[x])
    missing += len([y for y in msa[x] if y=='?'])
print(missing,total,missing/total)    

beta_prior(msa)


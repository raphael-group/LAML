#! /usr/bin/env python
import pickle
from problin_libs.sequence_lib import read_sequences
from problin_libs.ML_solver import ML_solver
from treeswift import *
import random
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-t","--topology",required=True,help="Input tree topology in newick format. Branch lengths will be ignored.")
parser.add_argument("-c","--characters",required=True,help="The input character matrix. Must have header.")
parser.add_argument("-p","--priors",required=False, default="uniform", help="The input prior matrix Q. Default: 'uniform'.")
parser.add_argument("--delimiter",required=False,default="tab",help="The delimiter of the input character matrix. Can be one of {'comma','tab','whitespace'} .Default: 'tab'.")
parser.add_argument("--nInitials",type=int,required=False,default=20,help="The number of initial points. Default: 20.")
parser.add_argument("-o","--output",required=True,help="The output file.")

args = vars(parser.parse_args())

delim_map = {'tab':'\t','comma':',','whitespace':' '}
delimiter = delim_map[args["delimiter"]]
msa = read_sequences(args["characters"],filetype="charMtrx",delimiter=delimiter)
with open(args["topology"],'r') as f:
    treeStr = f.read().strip()

k = len(msa[next(iter(msa.keys()))])

if args["priors"] == "uniform":
    # use the uniform Q matrix
    Q = []
    for i in range(k):
        m_i = len(set(msa[x][i] for x in msa if msa[x][i] not in [0,"?"]))
        q = {j+1:1/m_i for j in range(m_i)}
        q[0] = 0
        Q.append(q)
else:
    # read in the Q matrix
    priors = open(args["priors"], "rb")
    priors = pickle.load(priors)
    # key is site
    Q = []
    for i in sorted(priors.keys()):
        q = [priors[i][x] for x in sorted(priors[i])] 
        Q.append(q)
    # print(Q)


mySolver = ML_solver(msa,Q,treeStr)
optimal_llh = mySolver.optimize(initials=args["nInitials"])
with open(args["output"],'w') as fout:
    fout.write("Optimal tree: " +  mySolver.params.tree.newick() + "\n")
    fout.write("Optimal llh: " +  str(optimal_llh) + "\n")
    fout.write("Optimal dropout rate: " + str(mySolver.params.phi) + "\n")
    fout.write("Optimal silencing rate: " + str(mySolver.params.nu) + "\n")

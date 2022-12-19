#! /usr/bin/env python
import os.path
import pickle
from problin_libs.sequence_lib import read_sequences
from problin_libs.ML_solver import ML_solver
from treeswift import *
import random
import argparse

eps = 1e-10

parser = argparse.ArgumentParser()

parser.add_argument("-t","--topology",required=True,help="Input tree topology in newick format. Branch lengths will be ignored.")
parser.add_argument("-c","--characters",required=True,help="The input character matrix. Must have header.")
parser.add_argument("-r","--rep",required=False,help="The rep index of the input character matrix.") 
parser.add_argument("--noSilence",action='store_true',help="Assume there is no gene silencing, but allow missing data by dropout in sc-sequencing.")
parser.add_argument("--noDropout",action='store_true',help="Assume there is no sc-sequencing dropout, but allow missing data by gene silencing.")
parser.add_argument("-p","--priors",required=False, default="uniform", help="The input prior matrix Q. Default: if not specified, use a uniform prior.")
parser.add_argument("--delimiter",required=False,default="tab",help="The delimiter of the input character matrix. Can be one of {'comma','tab','whitespace'} .Default: 'tab'.")
parser.add_argument("--nInitials",type=int,required=False,default=20,help="The number of initial points. Default: 20.")
parser.add_argument("-o","--output",required=True,help="The output file.")

args = vars(parser.parse_args())

delim_map = {'tab':'\t','comma':',','whitespace':' '}
delimiter = delim_map[args["delimiter"]]
msa, site_names = read_sequences(args["characters"],filetype="charMtrx",delimiter=delimiter)
if args["rep"]:
    msa = read_sequences(args["characters"],filetype="fasta",delimiter=delimiter)
    print("Using rep:", args["rep"])
    msa = msa[int(args["rep"])]

if os.path.isfile(args["topology"]):
    with open(args["topology"],'r') as f:
        treeStr = f.read().strip()
else:
    treeStr = args["topology"]

k = len(msa[next(iter(msa.keys()))])
fixed_phi = eps if args["noDropout"] else None
fixed_nu = eps if args["noSilence"] else None

if args["priors"] == "uniform":
    # use the uniform Q matrix
    Q = []
    for i in range(k):
        M_i = set(msa[x][i] for x in msa if msa[x][i] not in [0,"?"])
        m_i = len(M_i)
        q = {x:1/m_i for x in M_i}
        q[0] = 0
        Q.append(q)
else:
    # read in the Q matrix
    file_extension = args["priors"].strip().split(".")[-1]
    if file_extension == "pkl": # pickled file
        infile = open(args["priors"], "rb")
        priors = pickle.load(infile)
        infile.close()
        Q = []
        for i in sorted(priors.keys()):
            q = {int(x):priors[i][x] for x in priors[i]}
            q[0] = 0
            Q.append(q)
    else:
        Q = [{0:0} for i in range(k)]
        with open(args["priors"],'r') as fin:
            for line in fin:
                site_idx,char_state,prob = line.strip().split()
                site_idx = int(site_idx)
                char_state = int(char_state)
                prob = float(prob)
                Q[site_idx][char_state] = prob

mySolver = ML_solver(msa,Q,treeStr)
optimal_llh = mySolver.optimize(initials=args["nInitials"],fixed_phi=fixed_phi,fixed_nu=fixed_nu)
with open(args["output"],'w') as fout:
    fout.write("Optimal tree: " +  mySolver.params.tree.newick() + "\n")
    fout.write("Optimal negative-llh: " +  str(optimal_llh) + "\n")
    fout.write("Optimal dropout rate: " + str(mySolver.params.phi) + "\n")
    fout.write("Optimal silencing rate: " + str(mySolver.params.nu) + "\n")

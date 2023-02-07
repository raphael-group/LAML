#! /usr/bin/env python
import pickle
from problin_libs.sequence_lib import read_sequences
from problin_libs.ML_solver import ML_solver
from problin_libs.EM_solver import EM_solver
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
parser.add_argument("-b","--betaPrior",required=False, default=(1,1), help="The beta prior of the dropout rate phi. Default: (1,1) --> uniform prior. Use 'auto' to let the solver automatically estimate the beta prior.")
parser.add_argument("--solver",required=False,default="Generic",help="Specify a solver. Options are 'Generic' or 'EM'. Caution: at current stage, EM only works with flag --noSilence. Default: 'Generic'")
parser.add_argument("--delimiter",required=False,default="tab",help="The delimiter of the input character matrix. Can be one of {'comma','tab','whitespace'} .Default: 'tab'.")
parser.add_argument("--nInitials",type=int,required=False,default=20,help="The number of initial points. Default: 20.")
parser.add_argument("-m","--maskedchar",required=True,help="Masked character.")
parser.add_argument("-o","--output",required=True,help="The output file.")
parser.add_argument("-v","--verbose",required=False,help="Print EM updates.",default=False)

args = vars(parser.parse_args())

delim_map = {'tab':'\t','comma':',','whitespace':' '}
delimiter = delim_map[args["delimiter"]]
msa, site_names = read_sequences(args["characters"],filetype="charMtrx",delimiter=delimiter,masked_symbol=args["maskedchar"])
if args["rep"]:
    print("Using rep:", args["rep"])
    msa = msa[int(args["rep"])]
with open(args["topology"],'r') as f:
    treeStr = f.read().strip()

k = len(msa[next(iter(msa.keys()))])
fixed_phi = eps if args["noDropout"] else None
fixed_nu = eps if args["noSilence"] else None
beta_prior = args["betaPrior"] 

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
    elif file_extension == "csv":
        Q = [{0:0} for i in range(k)]
        with open(args["priors"],'r') as fin:
            lines = fin.readlines()
            for line in lines[1:]:
                site_idx,char_state,prob = line.strip().split(',')
                site_idx = int(site_idx[1:])
                char_state = int(char_state)
                prob = float(prob)
                Q[site_idx][char_state] = prob
    else:
        Q = [{0:0} for i in range(k)]
        with open(args["priors"],'r') as fin:
            for line in fin:
                site_idx,char_state,prob = line.strip().split()
                site_idx = int(site_idx)
                char_state = int(char_state)
                prob = float(prob)
                Q[site_idx][char_state] = prob

selected_solver = ML_solver
em_selected = False
if args["solver"].lower() == "em": 
    if not args["noSilence"]:
        print("WARNING: EM algorithm has not been implemented for data with gene silencing. Program automatically switches to generic solver")
    else:
        selected_solver = EM_solver   
        em_selected = True
if em_selected:
    print("Optimization by EM algorithm") 
else:    
    print("Optimization by Generic solver")        
    
mySolver = selected_solver(msa,Q,treeStr,beta_prior=beta_prior)
optimal_llh = mySolver.optimize(initials=args["nInitials"],fixed_phi=fixed_phi,fixed_nu=fixed_nu,verbose=False)
with open(args["output"],'w') as fout:
    fout.write("Optimal tree: " +  mySolver.params.tree.newick() + "\n")
    fout.write("Optimal negative-llh: " +  str(optimal_llh) + "\n")
    fout.write("Optimal dropout rate: " + str(mySolver.params.phi) + "\n")
    fout.write("Optimal silencing rate: " + str(mySolver.params.nu) + "\n")

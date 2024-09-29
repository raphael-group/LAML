#! /usr/bin/env python
import os
import pickle
import laml_libs as laml
from laml_libs.IO_handler.sequence_lib import read_sequences, read_priors
from laml_libs.PMM_original.EM_solver import EM_solver
from math import *
from treeswift import *
import random
import argparse
import timeit
from sys import argv,exit,stdout,setrecursionlimit
import sys
from copy import deepcopy

setrecursionlimit(5000)

class Logger(object):
    def __init__(self, output_prefix):
        self.terminal = sys.stdout
        self.log = open(output_prefix + ".log", "a")
    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
    def flush(self):
        pass

def main():
    parser = argparse.ArgumentParser()
    otherOptions = parser._action_groups.pop()

    requiredNamed = parser.add_argument_group('required arguments')
    inputOptions= parser.add_argument_group('input options')
    outputOptions = parser.add_argument_group('output options')
    numericalOptions = parser.add_argument_group('numerical optimization options')
    parser._action_groups.append(otherOptions)

    # input arguments
    requiredNamed.add_argument("-t","--topology",required=True,help="[REQUIRED] The input tree topology in newick format. Branch lengths will be ignored.") 
    requiredNamed.add_argument("-c","--characters",required=True,help="[REQUIRED] The input character matrix. Must have header.")
    requiredNamed.add_argument("-n","--newSeqs",required=True,help="[REQUIRED] The new sequence to be placed onto the tree.")

    inputOptions.add_argument("-p","--priors",required=False, default="uniform", help="The input prior matrix Q. Default: if not specified, use a uniform prior.")
    inputOptions.add_argument("--delimiter",required=False,default="comma",help="The delimiter of the input character matrix. Can be one of {'comma','tab','whitespace'} .Default: 'comma'.")
    inputOptions.add_argument("-m","--missing_data",required=False,default="?",help="Missing data character. Default: if not specified, assumes '?'.")
    
    # output arguments
    outputOptions.add_argument("-o","--output",required=False,help="Output prefix. Default: LAML_output")
    outputOptions.add_argument("-v","--verbose",required=False,action='store_true',help="Show verbose messages.")
   
    # Numerical Optimization Arguments
    numericalOptions.add_argument("--randseeds",required=False,help="Random seeds for branch length optimization. Can be a single interger number or a list of intergers whose length is equal to the number of initial points (see --nInitials).")

    if len(argv) == 1:
        parser.print_help()
        exit(0)

    args = vars(parser.parse_args())
    if args["output"]:
        prefix = args["output"]
    else:
        prefix = "LAML_output"
    sys.stdout = Logger(prefix)

    lic_file = os.path.join(os.path.expanduser("~"), 'mosek/mosek.lic')
    if 'MOSEKLM_LICENSE_FILE' not in os.environ and not os.path.isfile(lic_file):
        print("MOSEK license not found in environment variables. Please set the MOSEK license!")
        exit(0)

    if not os.path.isfile(args["characters"]) or not os.path.isfile(args["topology"]):
        print("Input files not found.")
        exit(0)
    
    print("Launching " + laml.PROGRAM_NAME + " version " + laml.PROGRAM_VERSION)
    print(laml.PROGRAM_NAME + " was called as follows: " + " ".join(argv))
    
    # preprocessing: read and analyze input
    delim_map = {'tab':'\t','comma':',','whitespace':' '}
    delimiter = delim_map[args["delimiter"]]
    msa, site_names = read_sequences(args["characters"],filetype="charMtrx",delimiter=delimiter,masked_symbol=args["missing_data"])

    with open(args["topology"],'r') as f:
        input_trees = []
        for line in f:
            input_trees.append(line.strip())

    k = len(msa[next(iter(msa.keys()))])
    if args["randseeds"] is None:
        random_seeds = None
    else:
        random_seeds = [int(x) for x in args["randseeds"].strip().split()]
        if args["nInitials"] != 1 and len(random_seeds) == 1:
            random_seeds = random_seeds[0]

    if args["priors"] == "uniform":
        print("No prior file detected, using uniform prior probabilities for each alphabet on each site.")
        # use the uniform Q matrix
        Q = []
        for i in range(k):
            M_i = set(msa[x][i] for x in msa if msa[x][i] not in [0,"?"])
            #check if column has only zeros and missing data
            if len(M_i) == 0: 
                # add pseudo mutated state
                m_i = 1
                q = {"1":1.0}
            else:
                m_i = len(M_i)
                q = {x:1/m_i for x in M_i}
            q[0] = 0
            Q.append(q)
    else:
        Q = read_priors(args["priors"],msa,site_names=site_names)

    # main tasks        
    data = {'charMtrx':msa} 
    prior = {'Q':Q} 
    params = {'nu': laml.eps,'phi': laml.eps}  
    
    mySolver = EM_solver(input_trees,data,prior,params)
    start_time = timeit.default_timer()
    nllh,status = mySolver.optimize(initials=1,verbose=args["verbose"],random_seeds=random_seeds,ultra_constr=True)
    stop_time = timeit.default_timer()
    print("Optimization runtime (s):", stop_time - start_time)
    opt_params = {'nu':mySolver.params.nu,'phi':mySolver.params.phi}
    backbone_trees = mySolver.trees

    # read the new sequence
    newSeqs,_ = read_sequences(args["newSeqs"],filetype="charMtrx",delimiter=delimiter,masked_symbol=args["missing_data"])
    # add the new sequences to msa
    for new_cell in newSeqs:
        msa[new_cell] = newSeqs[new_cell]
    data = {'charMtrx':msa}
    
    # compute the tree height
    h = 0
    for tree in backbone_trees:
        tree.root.d2root = tree.root.edge_length
        for v in tree.traverse_preorder():
            if not v.is_root():
                v.d2root = v.parent.d2root + v.edge_length
            h = max(h,v.d2root)       
    
    placement_result = []
    for new_cell in newSeqs:
        best_nllh = float("inf")
        start_time = timeit.default_timer()
        for tree in backbone_trees:
            for v in tree.traverse_preorder():
                v_el = v.edge_length
                # create a new set of trees
                v_new = Node()
                v_new.label = new_cell
                u_new = Node()
                u_new.add_child(v_new)
                if v.is_root():
                    u_new.add_child(v)
                    tree.root = u_new
                else:    
                    u = v.parent
                    u.remove_child(v)
                    u.add_child(u_new)
                    u_new.add_child(v)
                v.edge_length = v_el/2
                u_new.edge_length = v_el/2
                v_new.edge_length = h - v.d2root + v_el/2
                new_trees = [t.newick() for t in backbone_trees]
                mySolver = EM_solver(new_trees,data,prior,opt_params)    
                curr_nllh = mySolver.negative_llh()
                if curr_nllh < best_nllh:
                    best_nllh = curr_nllh
                    best_trees = [t.newick() for t in mySolver.trees]
                # turn back
                u_new.remove_child(v)    
                if u_new.is_root():
                    tree.root = v
                else:
                    u.remove_child(u_new)
                    u.add_child(v)
                v.edge_length = v_el  
        stop_time = timeit.default_timer()
        print("Placement runtime (s):", stop_time - start_time)
        placement_result.append((new_cell,best_nllh,best_trees))

    with open(prefix+"_placement.txt",'w') as fout:
        for new_cell,best_nllh,best_trees in placement_result:
            fout.write("new cell: " + new_cell + "\n")
            fout.write("nllh: " + str(best_nllh) + "\n") 
            fout.write("augmented tree(s):\n")
            for tree_str in best_trees:
                fout.write(tree_str + "\n")   

if __name__ == "__main__":
    main()

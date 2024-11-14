#! /usr/bin/env python
import json
import os
import pickle
import laml_libs as laml
from laml_libs.IO_handler.sequence_lib import read_sequences, read_priors
from laml_libs.IO_handler.DLT_parser import DLT_parser
from laml_libs.Count_model.Alphabet import Alphabet
from laml_libs.Count_model.PMMN_model import PMMN_model
from laml_libs.Count_model.PMMC_model import PMMC_model 
from laml_libs.Count_model.CharMtrx import CharMtrx
from laml_libs.TopoSearch.Topology_search import Topology_search
from laml_libs.TopoSearch.Topology_search_parallel import Topology_search_parallel as Topology_search_parallel
from laml_libs.TopoSearch.Topology_search import Topology_search as Topology_search_sequential
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
        
def scale_trees_by_lambda(solver,Lambda): 
    for tree in solver.trees:
        for node in tree.traverse_preorder():
            if node.edge_length is not None:
                node.edge_length *= Lambda

def main():
    parser = argparse.ArgumentParser()
    otherOptions = parser._action_groups.pop()

    requiredNamed = parser.add_argument_group('required arguments')
    inputOptions= parser.add_argument_group('input options')
    outputOptions = parser.add_argument_group('output options')
    numericalOptions = parser.add_argument_group('numerical optimization options')
    topologySearchOptions = parser.add_argument_group('topology search options')
    parser._action_groups.append(otherOptions)

    # input arguments
    requiredNamed.add_argument("-t","--topology",required=True,help="[REQUIRED] The input tree topology in newick format. Branch lengths will be ignored.") 
    requiredNamed.add_argument("-c","--characters",required=True,help="[REQUIRED] The input character matrix. Must have header.")
    requiredNamed.add_argument("-M","--readout_model",required=True,help="[REQUIRED] The readout model.")

    inputOptions.add_argument("-p","--priors",required=False, default="uniform", help="The input prior matrix Q. Default: if not specified, use a uniform prior.")
    inputOptions.add_argument("--delimiter",required=False,default="comma",help="The delimiter of the input character matrix. Can be one of {'comma','tab','whitespace'} .Default: 'comma'.")
    inputOptions.add_argument("-m","--missing_data",required=False,default="?",help="Missing data character. Default: if not specified, assumes '?'.")
    inputOptions.add_argument("-s","--max_states_per_cassette",required=False,default=None,help="Max alphabet size per cassette. Default: if not specified, uses all of them.")
   
    # output arguments
    outputOptions.add_argument("-o","--output",required=False,help="Output prefix. Default: LAML_output")
    outputOptions.add_argument("-v","--verbose",required=False,action='store_true',help="Show verbose messages.")
   
    # Numerical Optimization Arguments
    numericalOptions.add_argument("--solver",required=False,default="EM",help="Specify a solver. Options are 'Scipy' or 'EM'. Default: EM")
    numericalOptions.add_argument("-L","--compute_llh",action='store_true',help="Compute log-likelihood of the input tree(s) in -t and the input params in -P (lambda,phi,nu,rho). Will NOT optimize branch lengths, lambda, phi, nu, or rho. The input tree(s) in -t MUST have branch lengths. This option has higher priority than --topology_search and --resolve_search.")
    numericalOptions.add_argument("--timescale",required=False,default=1.0,help="Experiment time. Scales the output tree height to this value. Default: 1.0.")
    numericalOptions.add_argument("--fixedBrlens",action='store_true',help="Keep the branch lengths of the input trees (specified in -t) fixed and optimize other numerical parameters. Preferably the genome edit rate (lambda) should be specified via -P; if lambda is not given, assume lambda=1. This option does NOT work with --topology_search.")
    numericalOptions.add_argument("--noSilence",action='store_true',help="Assume there is no gene silencing, but allow missing data by dropout in single cell sequencing.")
    numericalOptions.add_argument("--noDropout",action='store_true',help="Assume there is no sc-sequencing dropout, but allow missing data by gene silencing.")
    numericalOptions.add_argument("--noError",action='store_true',help="Assume there is no error in the input character matrix.")
    numericalOptions.add_argument("-P","--fixedParams",required=False,help="Fixed the params (lambda,phi,nu,rho) to the specified values. This option has higher priority than --noSilence, --noDropout, --noError, and --timescale. To use with -L, ALL of the four params (lambda,phi,nu,rho) must be specified")
    numericalOptions.add_argument("--nInitials",type=int,required=False,default=20,help="The number of initial points. Default: 20.")
    numericalOptions.add_argument("--randseeds",required=False,help="Random seeds for branch length optimization. Can be a single interger number or a list of intergers whose length is equal to the number of initial points (see --nInitials).")

    # Topology Search Arguments
    topologySearchOptions.add_argument("--topology_search",action='store_true', required=False,help="Perform topology search using NNI operations. Always returns a fully resolved (i.e. binary) tree.")
    topologySearchOptions.add_argument("--resolve_search",action='store_true', required=False,help="Resolve polytomies by performing topology search ONLY on branches with polytomies. This option has higher priority than --topology_search.")
    topologySearchOptions.add_argument("--keep_polytomies",action='store_true', required=False,help="Keep polytomies while performing topology search. This option only works with --topology_search.")
    topologySearchOptions.add_argument("--randomreps", required=False, default=1, type=int, help="Number of replicates to run for the random strategy of topology search.")
    topologySearchOptions.add_argument("--maxIters", required=False, default=500, type=int, help="Maximum number of iterations to run topology search.")
    topologySearchOptions.add_argument("--parallel", required=False,action='store_true', help="Turn on parallel version of topology search.")

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
    start_time = timeit.default_timer()
    
    # preprocessing: read and analyze input
    delim_map = {'tab':'\t','comma':',','whitespace':' '}
    delimiter = delim_map[args["delimiter"]]


    with open(args["topology"],'r') as f:
        input_trees = []
        for line in f:
            input_trees.append(line.strip())

    fixed_params = {}
    if args["noDropout"]:
        fixed_params['phi'] = 0
    if args["noSilence"]:
        fixed_params['nu'] = 0
    if args["noError"]:
        fixed_params['rho'] = 1
    if args["fixedBrlens"]:
        fixed_params["lambda"] = 1 # default value; will be overrided by args["fixedParams"] if lambda is specified there

    if args["fixedParams"]: # will override "noDropout", "noSilence", "noError" if 'phi', 'nu', or 'rho' presents in fixedParams`
        tokens = args["fixedParams"].strip().split()
        for token in tokens:
            key,value = token.strip().split("=")
            fixed_params[key] = float(value)

    if args["randseeds"] is None:
        random_seeds = None
    else:
        random_seeds = [int(x) for x in args["randseeds"].strip().split()]
        if args["nInitials"] != 1 and len(random_seeds) == 1:
            random_seeds = random_seeds[0]

    priorfile = None if args["priors"] == "uniform" else args["priors"]
    parser = DLT_parser(datafile=args["characters"], priorfile=priorfile, max_allele_per_cassette=args["max_states_per_cassette"])
    selected_model = PMMN_model if args["readout_model"] == "PMMN" else PMMC_model
    data = {'DLT_data': parser.DLT_data} 
    prior = {'Q': parser.priors}  

    # main tasks        
    #ini_phi = 0 if fixed_phi is None else fixed_phi
    #ini_nu = 0 if fixed_nu is None else fixed_nu
    #ini_params = {'phi':0.1,'nu':0,'mu':1,'rho':1}
    Topology_search = Topology_search_sequential if not args["parallel"] else Topology_search_parallel
    #myTopoSearch = Topology_search(input_trees, selected_model, data=data, prior=prior, params=ini_params)
    myTopoSearch = Topology_search(input_trees, selected_model, data=data, prior=prior, params=fixed_params)

    print("Setup finished. Beginning computation...")
    if args["compute_llh"]:
        print("Computing the joint likelihood of the input trees and specified parameters without any optimization")
        mySolver = myTopoSearch.get_solver()
        # [hacking] rescale the input branch lengths by the specified lambda
        fixed_lambda = fixed_params['lambda'] #if 'lambda' in fixed_params else 1
        scale_trees_by_lambda(mySolver,fixed_lambda) 
        #for tree in mySolver.trees:
        #    for node in tree.traverse_preorder():
        #        if node.edge_length is not None:
        #            node.edge_length *= fixed_lambda
        nllh = mySolver.negative_llh()
        #opt_trees = myTopoSearch.treeTopoList
        #opt_params = myTopoSearch.params
        print("Tree negative log-likelihood: " + str(nllh))
        print("Tree log-likelihood: " + str(-nllh))
    else:
        # setup the strategy
        my_strategy = deepcopy(laml.DEFAULT_STRATEGY)
        # enforce ultrametric or not?
        my_strategy['ultra_constr'] = True
        # resolve polytomies or not?
        resolve_polytomies = not args["keep_polytomies"]
        # only resolve polytomies or do full search?
        my_strategy['resolve_search_only'] = args["resolve_search"]
        # fixed params
        my_strategy['fixed_params'] = deepcopy(fixed_params)

        # full search or local search to only resolve polytomies? 
        if not args["resolve_search"] and not args["topology_search"]:
            mySolver = myTopoSearch.get_solver()
            fixed_brlen = "All" if args["fixedBrlens"] else None
            if args["fixedBrlens"]:
                print("Fixing branch lengths to the given values and optimizing other numerical parameters without topology search")
                # [hacking] rescale the input branch lengths by the specified lambda
                fixed_lambda = fixed_params['lambda'] #if 'lambda' in fixed_params else 1
                scale_trees_by_lambda(mySolver,fixed_lambda) 
            else:    
                print("Optimizing branch lengths and other numerical parameters without topology search")
            nllh = mySolver.optimize('EM',initials=args["nInitials"],fixed_params=fixed_params,verbose=args["verbose"],random_seeds=random_seeds,ultra_constr=True,fixed_brlen=fixed_brlen)
            myTopoSearch.update_from_solver(mySolver)
            opt_trees = myTopoSearch.treeTopoList
            opt_params = myTopoSearch.params
        else:
            if args["resolve_search"]:
                if not resolve_polytomies:
                    print("WARNING: --resolve_search was specified with --keep_polytomies. Program will only optimize numerical parameters WITHOUT any topology search.")
                else:    
                    print("Starting local topology search to resolve polytomies")
                    if not myTopoSearch.has_polytomies:
                        print("No polytomy detected. Program will only optimize numerical parameters WITHOUT any topology search.")
            else:
                if myTopoSearch.has_polytomy:
                    print("Detected polytomies in the input trees.")
                    if not resolve_polytomies:
                        print("Flag --keep_polytomies is on. All polytomies will be kept.")                 
                    else:
                        print("Flag --keep_polytomies is off. All polytomies will be resolved.") 
                else:
                    print("No polytomy detected.")         
                print("Starting topology search")

            if args["parallel"]:
                print("Running topology search in parallel...")
            else:
                print("Running topology search sequentially...")
            checkpoint_file = f"{prefix}_ckpt.txt"
            opt_trees,max_score,opt_params = myTopoSearch.search(resolve_polytomies=resolve_polytomies,maxiter=args["maxIters"], verbose=args["verbose"], strategy=my_strategy, nreps=args['randomreps'],checkpoint_file=checkpoint_file) 
            nllh = -max_score        
    
    # post-processing: analyze results and output 
    if not args["compute_llh"]:
        out_tree = prefix + "_trees.nwk"
        out_params = prefix + "_params.txt"
        
        with open(out_tree,'w') as fout:
            for tstr in opt_trees:
                tree = read_tree_newick(tstr)
                # get the height of the tree
                tree_height = tree.height(weighted=True) # includes the root's length, edit units 
                # get the edit rate
                if 'lambda' in fixed_params:
                    edit_rate = fixed_params['lambda']
                    print(f"Tree height pre-scaling: {tree_height}, input edit rate: {edit_rate}") 
                else:    
                    edit_rate = tree_height/float(args['timescale'])
                    print(f"Tree height pre-scaling: {tree_height}, input timescale: {args['timescale']}") 
                for node in tree.traverse_preorder(): 
                    if node.edge_length:
                        node.edge_length = node.edge_length / edit_rate
                tree_height = tree.height(weighted=True) 
                print(f"Tree height after scaling: {tree_height}")
                if len(tree.root.children) > 1:
                    new_root = Node()
                    new_root.add_child(tree.root)
                    tree.root = new_root
                fout.write(tree.newick() + "\n")
        
        with open(out_params,'w') as fout:
            fout.write("Dropout probability: phi=" + str(opt_params['phi']) + "\n")
            fout.write("Silencing rate: nu=" + str(opt_params['nu']) + "\n") 
            fout.write("Sequencing accuracy: rho=" + str(opt_params['rho']) + "\n") 
            fout.write("Mutation rate: lambda=" +  str(edit_rate) + "\n") 
            fout.write("Negative-llh: " +  str(nllh) + "\n")
                    
        stop_time = timeit.default_timer()
        print("Runtime (s):", stop_time - start_time)

if __name__ == "__main__":
    main()

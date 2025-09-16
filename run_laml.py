#! /usr/bin/env python


from laml_libs.runner import run_from_namespace

import os
import pickle
import laml_libs as laml
#from laml_libs.sequence_lib import read_sequences, read_priors, dedup, add_dup
#from laml_libs.ML_solver import ML_solver
#from laml_libs.EM_solver import EM_solver
#from laml_libs.fastEM_solver import fastEM_solver, parse_data, parse_tree
#from laml_libs.Topology_search_parallel import Topology_search_parallel as Topology_search_parallel
#from laml_libs.Topology_search import Topology_search as Topology_search_sequential
#from laml_libs.starting_tree import build_starting_tree
from math import *
from treeswift import *
import datetime
from datetime import date
import re
import math
import random
import argparse
import mosek
import timeit
from sys import argv,exit,stdout,setrecursionlimit
import sys
from copy import deepcopy
import jax

setrecursionlimit(5000)


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
    requiredNamed.add_argument("-c","--characters",required=True,help="[REQUIRED] The input character matrix. Must have header.")

    inputOptions.add_argument("-t","--topology",required=False,help="The input tree topology in newick format. Branch lengths will be ignored. If you do not provide an input tree topology, we will use NJ with weighted Hamming distances.") 
    inputOptions.add_argument("-p","--priors",required=False, default="uniform", help="The input prior matrix Q. Default: if not specified, use a uniform prior.")
    inputOptions.add_argument("--delimiter",required=False,default="comma",help="The delimiter of the input character matrix. Can be one of {'comma','tab','whitespace'} .Default: 'comma'.")
    inputOptions.add_argument("-m","--missing_data",required=False,default="?",help="Missing data character. Default: if not specified, assumes '?'.")
    
    # output arguments
    outputOptions.add_argument("-o","--output",required=False,help="Output prefix. Default: LAML_output")
    outputOptions.add_argument("-v","--verbose",required=False,action='store_true',help="Show verbose messages.")
   
    # Numerical Optimization Arguments
    numericalOptions.add_argument("--solver",required=False,default="EM",help="Specify a solver. Options are 'Scipy' or 'EM' or 'fastEM-gpu' or 'fastEM-cpu'. Default: EM")
    numericalOptions.add_argument("-L","--compute_llh",required=False,help="Compute likelihood of the input tree using the input (lambda,phi,nu). Will NOT optimize branch lengths, lambda, phi, or nu. The input tree MUST have branch lengths. This option has higher priority than --topology_search and --resolve_search.")
    numericalOptions.add_argument("--timescale",required=False,default=1.0,help="Timeframe of experiment. Scales ultrametric output tree branches to this timescale. To get an accurate estimate of mutation rate, provide timeframe in number of cell generations. Default: 1.0.")
    numericalOptions.add_argument("--noSilence",action='store_true',help="Assume there is no gene silencing, but allow missing data by dropout in single cell sequencing.")
    numericalOptions.add_argument("--noDropout",action='store_true',help="Assume there is no sc-sequencing dropout, but allow missing data by gene silencing.")
    #numericalOptions.add_argument("--noultrametric",action='store_true',help="[Deprecated] But turns OFF the ultrametric constraint.")
    numericalOptions.add_argument("--nInitials",type=int,required=False,default=20,help="The number of initial points. Default: 20.")
    numericalOptions.add_argument("--randseeds",required=False,help="Random seeds for branch length optimization. Can be a single interger number or a list of intergers whose length is equal to the number of initial points (see --nInitials).")
    numericalOptions.add_argument("--gpu",action='store_true',required=False, default=False, help="Runs on discoverable GPUs otherwise throws error.")

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

    args = parser.parse_args()
    #args = vars(parser.parse_args())
    tree_str, params, annotations_path, imputed_matrix_path, log_path = run_from_namespace(args)

if __name__ == "__main__":
    main()

#! /usr/bin/env python
import os
import pickle
import scmail_libs as scmail
from scmail_libs.sequence_lib import read_sequences, read_priors
from scmail_libs.ML_solver import ML_solver
from scmail_libs.EM_solver import EM_solver
from scmail_libs.Topology_search_parallel import Topology_search_parallel as Topology_search_parallel
from scmail_libs.Topology_search import Topology_search as Topology_search_sequential
from math import *
from treeswift import *
import random
import argparse
import timeit
from sys import argv,exit,stdout,setrecursionlimit
from copy import deepcopy

setrecursionlimit(5000)

def main():
    parser = argparse.ArgumentParser()

    # input arguments
    parser.add_argument("-t","--topology",required=True,help="Binary input tree topology in newick format. Branch lengths will be ignored.") 
    parser.add_argument("-c","--characters",required=True,help="The input character matrix. Must have header.")
    parser.add_argument("-p","--priors",required=False, default="uniform", help="The input prior matrix Q. Default: if not specified, use a uniform prior.")
    parser.add_argument("--delimiter",required=False,default="tab",help="The delimiter of the input character matrix. Can be one of {'comma','tab','whitespace'} .Default: 'tab'.")
    parser.add_argument("-m","--maskedchar",required=False,default="-",help="Masked character. Default: if not specified, assumes '-'.")
    parser.add_argument("-o","--output",required=True,help="Output prefix.")
   
    # which problem are you solving? 
    parser.add_argument("--solver",required=False,default="EM",help="Specify a solver. Options are 'Scipy' or 'EM'. Default: EM")
    parser.add_argument("--topology_search",action='store_true', required=False,help="Perform topology search using NNI operations. Always return fully resolved (i.e. binary) tree.")
    parser.add_argument("--resolve_search",action='store_true', required=False,help="Resolve polytomies by performing topology search ONLY on branches with polytomies. This option has higher priority than --topology_search.")
    parser.add_argument("--keep_polytomies",action='store_true', required=False,help="Keep polytomies while performing topology search. This option only works with --topology_search.")
    parser.add_argument("-L","--compute_llh",required=False,help="Compute likelihood of the input tree using the input (phi,nu). Will NOT optimize branch lengths, phi, or nu. The input tree MUST have branch lengths. This option has higher priority than --topology_search and --resolve_search.")

    # problem formulation
    parser.add_argument("--ultrametric",action='store_true',help="Enforce ultrametricity to the output tree.")
    parser.add_argument("--noSilence",action='store_true',help="Assume there is no gene silencing, but allow missing data by dropout in sc-sequencing.")
    parser.add_argument("--noDropout",action='store_true',help="Assume there is no sc-sequencing dropout, but allow missing data by gene silencing.")

    # miscellaneous
    parser.add_argument("-v","--verbose",required=False,action='store_true',help="Show verbose messages.")
    parser.add_argument("--nInitials",type=int,required=False,default=20,help="The number of initial points. Default: 20.")
    parser.add_argument("--randseeds",required=False,help="Random seeds. Can be a single interger number or a list of intergers whose length is equal to the number of initial points (see --nInitials).")
    parser.add_argument("--randomreps", required=False, default=1, type=int, help="Number of replicates to run for the random strategy of topology search.")
    parser.add_argument("--maxIters", required=False, default=500, type=int, help="Maximum number of iterations to run topology search.")
    parser.add_argument("--parallel", required=False,action='store_true', help="Turn on parallel version of topology search.")

    if len(argv) == 1:
        parser.print_help()
        exit(0)

    if 'MOSEKLM_LICENSE_FILE' not in os.environ:
        print("MOSEK license not found in environment variables. Please set the MOSEK license!")
        exit(0)
    
    print("Launching " + scmail.PROGRAM_NAME + " version " + scmail.PROGRAM_VERSION)
    print(scmail.PROGRAM_NAME + " was called as follow: " + " ".join(argv))
    start_time = timeit.default_timer()
    
    args = vars(parser.parse_args())
    
    # preprocessing: read and analyze input
    delim_map = {'tab':'\t','comma':',','whitespace':' '}
    delimiter = delim_map[args["delimiter"]]
    msa, site_names = read_sequences(args["characters"],filetype="charMtrx",delimiter=delimiter,masked_symbol=args["maskedchar"])
    #prefix = '.'.join(args["output"].split('.')[:-1])
    prefix = args["output"]

    with open(args["topology"],'r') as f:
        input_trees = []
        for line in f:
            input_trees.append(line.strip())

    k = len(msa[next(iter(msa.keys()))])
    if args["compute_llh"]:
        fixed_phi,fixed_nu = [float(x) for x in args["compute_llh"].strip().split()]
    else:    
        fixed_phi = 0 if args["noDropout"] else None
        fixed_nu = 0 if args["noSilence"] else None

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
            m_i = len(M_i)
            q = {x:1/m_i for x in M_i}
            q[0] = 0
            Q.append(q)
    else:
        Q = read_priors(args["priors"], site_names)

    selected_solver = EM_solver
    em_selected = True
    if args["solver"].lower() != "em": 
        selected_solver = ML_solver   
        em_selected = False

    # main tasks        
    data = {'charMtrx':msa} 
    prior = {'Q':Q} 
    
    params = {'nu':fixed_nu if fixed_nu is not None else scmail.eps,'phi':fixed_phi if fixed_phi is not None else scmail.eps}  
    Topology_search = Topology_search_sequential if not args["parallel"] else Topology_search_parallel


    myTopoSearch = Topology_search(input_trees, selected_solver, data=data, prior=prior, params=params)


    if args["compute_llh"]:
        print("Compute likelihood of the input tree and specified parameters without any optimization")
        mySolver = myTopoSearch.get_solver()
        nllh = mySolver.negative_llh()
        opt_trees = myTopoSearch.treeTopoList
        opt_params = myTopoSearch.params
        print("Tree negative log-likelihood: " + str(nllh))
        print("Tree log-likelihood: " + str(-nllh))
    else:
        # setup the strategy
        my_strategy = deepcopy(scmail.DEFAULT_STRATEGY)
        # enforce ultrametric or not?
        my_strategy['ultra_constr'] = args["ultrametric"]
        # resolve polytomies or not?
        resolve_polytomies = not args["keep_polytomies"]
        # only resolve polytomies or do full search?
        my_strategy['resolve_search_only'] = args["resolve_search"]
        # full search or local search to only resolve polytomies? 
        if not args["resolve_search"] and not args["topology_search"]:
            print("Optimizing branch lengths, phi, and nu without topology search")
            if em_selected:
                print("Optimization by EM algorithm") 
            else:    
                print("Optimization by generic solver (Scipy-SLSQP)")        
            mySolver = myTopoSearch.get_solver()
            nllh = mySolver.optimize(initials=args["nInitials"],fixed_phi=fixed_phi,fixed_nu=fixed_nu,verbose=args["verbose"],random_seeds=random_seeds,ultra_constr=args["ultrametric"])      
            myTopoSearch.update_from_solver(mySolver)
            opt_trees = myTopoSearch.treeTopoList
            opt_params = myTopoSearch.params
        else:
            if args["resolve_search"]:
                if not resolve_polytomies:
                    print("WARNING: --resolve_search was specified with --keep_polytomies. Program will only optimize numerical parameters WITHOUT any topology search.")
                else:    
                    print("Starting local topology search to resolve polytomies")
            else:
                if not resolve_polytomies:
                    print("Keeping all the polytomies")                 
                else:
                    print("All polytomies will be resolved")    
                print("Starting topology search")
            if args["parallel"]:
                print("Running topology search in parallel...")
            else:
                print("Running topology search sequentially...")
            randval = int(random.random() * 1000)
            checkpoint_file = f"{prefix}._ckpt.{randval}.txt"
            opt_trees,max_score,opt_params = myTopoSearch.search(resolve_polytomies=resolve_polytomies,maxiter=args["maxIters"], verbose=args["verbose"], strategy=my_strategy, nreps=args['randomreps'],checkpoint_file=checkpoint_file) 
            nllh = -max_score        
    
    # post-processing: analyze results and output 
    out_tree = prefix + "_trees.nwk"
    out_annotate = prefix + "_annotations.txt"
    out_params = prefix + "_params.txt"

    with open(out_tree,'w') as fout:
        for tree in opt_trees:
            fout.write(tree + "\n")

    # output annotations
    def format_posterior(p0,p_minus_1,p_alpha,alpha,q):
        if p0 == 1:
            return '0'
        elif p_minus_1 == 1:
            return '-1'
        elif p_alpha == 1:
            return alpha
        else:
            out = ''
            if p0 > 0:
                out += '0:' + str(p0)
            if p_minus_1 > 0:
                if out != '':
                    out += '/'
                out += '-1:' + str(p_minus_1)
            if p_alpha > 0:
                if out != '':
                    out += '/'
                if alpha == '?':
                    out += "/".join([str(y) + ':' + str(round(p_alpha*q[y],3)) for y in q if round(p_alpha*q[y],3)>0])
                else:
                    out += alpha + ":" + str(p_alpha)
            return out

    my_solver = EM_solver(opt_trees,{'charMtrx':msa},{'Q':Q},{'phi':opt_params['phi'],'nu':opt_params['nu']})
    my_solver.az_partition()
    my_solver.Estep()
    idx = 0
    with open(out_annotate,'w') as fout:
        for tree in my_solver.trees:
            # add root
            if len(tree.root.children) > 1:
                root = Node()
                root.label = 'I0'
                k = len(tree.root.alpha)
                root.alpha = ['z']*k
                root.post0 = [0]*k
                root.post1 = [-float("inf")]*k
                idx = 1
                root.add_child(tree.root)
                tree.root = root
            # branch length by expected number of mutations
            all_labels = set()
            for node in tree.traverse_preorder():
                if node.is_root():
                    continue
                # label the node
                if node.label is None or node.label in all_labels:
                    node.label = 'I' + str(idx)
                    idx += 1                    
                all_labels.add(node.label)
                node.edge_length = round(sum(node.S1)+sum(node.S2)+sum(node.S4),3)
            fout.write(tree.newick()+"\n")    
            
            # ancestral labeling and imputation
            for node in tree.traverse_preorder():
                node.posterior = ''
                for j in range(len(node.alpha)):
                    if not node.is_root():
                        if node.alpha[j] == '?' and node.parent.alpha[j] != 'z':
                            node.alpha[j] = node.parent.alpha[j]
                    p0 = round(exp(node.post0[j]),2)
                    p_minus_1 = round(exp(node.post1[j]),2)
                    p_alpha = round(1-p0-p_minus_1,2)
                    if node.posterior != '':
                        node.posterior += ','
                    node.posterior += format_posterior(p0,p_minus_1,p_alpha,str(node.alpha[j]),Q[j])
                fout.write(node.label+"," + str(node.posterior)+"\n")
    
    with open(out_params,'w') as fout:
        fout.write("Dropout rate: " + str(opt_params['phi']) + "\n")
        fout.write("Silencing rate: " + str(opt_params['nu']) + "\n") 
        fout.write("Negative-llh: " +  str(nllh) + "\n")
                
    stop_time = timeit.default_timer()
    print("Runtime (s):", stop_time - start_time)

if __name__ == "__main__":
    main()

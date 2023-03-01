#! /usr/bin/env python
from problin_libs.sequence_lib import read_sequences, read_Q, write_sequences
from problin_libs.ML_solver import ML_solver
from problin_libs.EM_solver import EM_solver
from problin_libs.sim_lib import *
from treeswift import *
import random
import argparse

#random.seed(a=7)

# e.g. python simulate_seq.py -e --height 8 -b 0.5 -k 30 -r 1 --mu 0.1 -d 0.1 -s 0.1 --prefix "tmp" 

parser = argparse.ArgumentParser()

parser.add_argument("-e", "--balanced", action='store_true', help="Generate balanced tree of height h and branch length of length b.")
parser.add_argument("-b", "--branchlen", type=float, required=False, default=0.5, help="Generate balanced tree with branch lengths of length b. Default is 0.5.")
parser.add_argument("--height", type=int, required=False, default=8, help="The height of the balanced tree. Default is 8 levels.")
parser.add_argument("-n", "--numcells", required=False, help="the number of cells (leaves) in the simulated tree.")
parser.add_argument("-k", "--seqlen", type=int, required=True, help="The sequence length in each cell.")
parser.add_argument("--alphabetsize", default=10, type=int, required=False, help="The alphabet size for each site. The default is 10.")

parser.add_argument("-t","--topology",required=False,help="Binary input tree topology in newick format. Required if balanced tree flag is not set.")
parser.add_argument("--priors",required=False,help="Prior mutation probabilities file.")
parser.add_argument("-r","--reps",type=int, required=False,help="The number of replicates to create.")
parser.add_argument("--mu", type=float, required=False, default=0.1,  help="Mutation rate. Defaults to 0.1.")

parser.add_argument("-d", "--dropout", type=float, required=False, default=0.0,  help="The dropout rate, default is 0.")
parser.add_argument("-s", "--silencing", type=float, required=False, default=0.0,  help="The heritable silencing rate, default is 0.")
parser.add_argument("-p","--prefix",required=True,help="The prefix of all output files.")

# TODO: add random seed
#parser.add_argument("--randseeds",required=False,help="Random seeds. Can be a single interger number or a list of intergers whose length is equal to the number of initial points (see --nInitials).")
parser.add_argument("--maskedchar",required=False,default="?",help="Masked character. Default: if not specified, assumes '?'.")


args = vars(parser.parse_args())

tree_outfile = args["prefix"] + "_tree_b" + str(args["branchlen"]) + ".nwk"
if args["balanced"]:
    treeStr = get_balanced_tree(args["height"], args["branchlen"])
    with open(tree_outfile, "w") as fout:
        fout.write(treeStr)
    tree = read_tree_newick(treeStr)
else:
    tree = read_tree_newick(args["topology"])
# add cassiopeia tree option.

if args["priors"]:
    Q = read_Q(args["priors"])
    pass
else:
    prior_outfile = args["prefix"] + "_priors.csv"
    Q = sim_Q(args["seqlen"], args["alphabetsize"], prior_outfile=prior_outfile)

mu = args["mu"]
sim_dropout, sim_silencing = False, False
if args["dropout"] != 0.0:
    sim_dropout = True
    d = args["dropout"]
else:
    d = 0.0
if args["silencing"] != 0.0:
    sim_silencing = True
    s = args["silencing"]
else:
    s = None

mchar = args["maskedchar"]
k = args["seqlen"]
for i in range(args["reps"]):
    char_mtrx = simulate_seqs(tree, Q, mu, with_heritable=sim_silencing, silencing_rate=s, s=None)

    if not sim_dropout and not sim_silencing:
        nomissing_outfile = args["prefix"] + "_r" + str(i) + "_nomissing_charmtrx.csv"
        write_sequences(char_mtrx, k, nomissing_outfile)

    elif not sim_dropout and sim_silencing:
        silencing_outfile = args["prefix"] + "_r" + str(i) + "_s" + str(s) + "_character_matrix.csv"
        write_sequences(char_mtrx, k, silencing_outfile)

    elif not sim_silencing and sim_dropout:
        nmtx = simulate_dropout(char_mtrx, d)
        dropout_outfile = args["prefix"] + "_r" + str(i) + "_d" + str(d) + "_character_matrix.csv"
        write_sequences(char_mtrx, k, dropout_outfile)

    elif sim_silencing and sim_dropout:
        both_outfile = args["prefix"] + "_r" + str(i) + "_s" + str(s) + "_d" + str(d) + "_character_matrix.csv"
        write_sequences(char_mtrx, k, both_outfile)

# TODO: Call functions from unit_test/utils.py to check the proportion of missing data simulated and print statistics?
.

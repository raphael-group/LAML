#! /usr/bin/env python
from problin_libs.sequence_lib import read_sequences, read_Q, write_sequences
from problin_libs.ML_solver import ML_solver
from problin_libs.EM_solver import EM_solver
from problin_libs.sim_lib import *
from treeswift import *
import random
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-k", "--seqlen", type=int, required=True, help="The sequence length in each cell.")
parser.add_argument("--alphabetsize", default=10, type=int, required=False, help="The alphabet size for each site. The default is 10.")

parser.add_argument("-t","--inputTree",required=True,help="Input tree topology with branch length (binary tree required).")
parser.add_argument("--priors",required=False,help="Prior mutation probabilities file.")
parser.add_argument("-r","--reps",type=int, required=False,help="The number of replicates to create.")
parser.add_argument("--mu", type=float, required=False, default=1,  help="Mutation rate. Defaults to 1.")

parser.add_argument("-d", "--dropout", type=float, required=False, default=0.0,  help="The dropout rate, default is 0.")
parser.add_argument("-s", "--silencing", type=float, required=False, default=0.0,  help="The heritable silencing rate, default is 0.")
parser.add_argument("-p","--prefix",required=True,help="The prefix of all output files.")
parser.add_argument("--randseed",required=False,type=int,help="Random seed: an interger number.")
parser.add_argument("--maskedchar",required=False,default="?",help="Masked character. Default: if not specified, assumes '?'.")


args = vars(parser.parse_args())
tree = read_tree_newick(args["inputTree"])

i = 0
for node in tree.traverse_preorder():
    if not node.is_leaf():
        node.label = "I" + str(i)
        i += 1

k = args["seqlen"]

if args["priors"]:
    Q = read_Q(args["priors"])
else:
    prior_outfile = args["prefix"] + "_priors.csv"
    Q = sim_Q(args["seqlen"], args["alphabetsize"])
    with open(prior_outfile, "w") as fout:
        for i in range(k):
            for x in Q[i]:
                fout.write(str(i) + "," + str(x) + "," + str(Q[i][x]) + "\n")

mu = args["mu"]
sim_dropout, sim_silencing = False, False

if args["dropout"] != 0.0:
    sim_dropout = True
    d = args["dropout"]
else:
    d = 0.0

sim_silencing = args["silencing"] != 0    
s = args["silencing"]

if args["randseed"] is not None:
    random.seed(args["randseed"])

nreps = args["reps"]
for i in range(nreps):
    leaf_char_mtrx,all_char_mtrx = simulate_seqs(tree, Q, mu, with_heritable=sim_silencing, silencing_rate=s, dropout_rate=d)
    out_seqs = args["prefix"] + "_r" + str(i+1).rjust(len(str(nreps)),'0') + "_character_matrix.csv"
    out_history = args["prefix"] + "_r" + str(i+1).rjust(len(str(nreps)),'0') + "_all_sequences.csv"
    write_sequences(leaf_char_mtrx, k, out_seqs)
    write_sequences(all_char_mtrx, k, out_history)
    with open(out_history,'a') as fout:
        fout.write("Evolve tree: " + tree.newick())

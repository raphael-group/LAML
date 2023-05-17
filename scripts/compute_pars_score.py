from problin_libs.sequence_lib import read_sequences, read_Q 
from scripts.preprocess import load_pickle
import startle # TODO: REMOVE BEFORE RELEASE!
import argparse
from treeswift import *
from math import log
import pandas as pd
import time

def pars_score_startle(args):
    seed_tree = startle.from_newick_get_nx_tree(args.tree)
    character_matrix = pd.read_csv(args.msa, index_col=[0], sep=args.delimiter, dtype=str)
    character_matrix = character_matrix.replace('-', '-1')
    #print(character_matrix)
    if args.prior != '':
        mutation_prior_dict = load_pickle(args.prior)
        #mutation_prior_dict = read_Q(args.prior)
        weighted_seed_parsimony, _, _ = startle.small_parsimony(mutation_prior_dict, seed_tree, character_matrix, weighted=True)
        return weighted_seed_parsimony
    else:
        mutation_prior_dict = []
        seed_parsimony, _, _ = startle.small_parsimony(mutation_prior_dict, seed_tree, character_matrix, weighted=False)
        return seed_parsimony

def pars_score(T, msa, mchar, norm, priorfile):
    def get_seq(n, nodedict):
        #print(n.label)
        if n.is_leaf():
            return msa[n.label]
        else:
            return nodedict[n]

    def score_fn(c, cidx, use_weighted, q):
        #print("score_fn")
        if use_weighted:
            #print("prior keys", q.keys())
            #print("cidx", cidx, "c", c)
            #print("q[cidx][c]", q[cidx][c])
            if c in q[cidx]:
                if q[cidx][c] > 0:
                    return -log(q[cidx][c])
            else:
                return 0
        else:
            return 1

    m, site_names = read_sequences(msa,delimiter=",",masked_symbol=mchar, suppress_warnings=True)
    msa = m
    site_names = sorted([int(x[1:]) for x in site_names])
    #print("MSA:", msa)
    q = []
    
    use_weighted = False
    if priorfile != "":
        use_weighted = True
        q = load_pickle(priorfile) #read_Q(priorfile)
        #q = read_Q(priorfile)

        # check that site_names == q.keys()

        prior_names = sorted([x for x in q.keys()])
        #if prior_names != site_names: 
        #    print("Provided prior names do not match the character site names.")
        #    print("Prior names:", prior_names)
        #    print("Site names:", site_names)


    nodedict = dict()
    score = 0
    num_internal = 0
    for n in T.traverse_postorder():
        if n.is_leaf():
            nodedict[n] = get_seq(n, nodedict)
            #print(n.label, nodedict[n])
        else:
            num_internal += 1
            a, b = n.child_nodes()
            # print(n.label, "num children", len(n.child_nodes()))
            s1 = get_seq(a, nodedict)
            s2 = get_seq(b, nodedict)
            # print(s1, s2)

            s = []
            # handle missing data
            # TODO: Handle root
            for cidx, xy in enumerate(zip(s1, s2)):
                x, y = xy
                if use_weighted:
                    cidx = prior_names[cidx]
                if x == y: 
                    s.append(x)
                elif x == "?" or y == "?":
                    if x == "?":
                        s.append(y) # the one that's not mchar
                    else:
                        s.append(x)
                else: # x != y
                    # check if one is 0 and alpha
                    if x == 0 or y == 0: #'0' or y == '0':
                        if x == 0:
                            score += score_fn(x, cidx, use_weighted, q)
                        elif y == 0: 
                            score += score_fn(y, cidx, use_weighted, q)
                    else:
                        score += score_fn(x, cidx, use_weighted, q)
                        score += score_fn(y, cidx, use_weighted, q)
                    s.append(0)
            nodedict[n] = s
            # print(n.label, s)

    # ensure the root is 0
    for cidx, c in enumerate(s):
        if c != 0 and c != '?':
            # print("root")
            score += score_fn(c, cidx, use_weighted, q)

    #print("nodedict")
    #for n in nodedict:
    #    print(n.label, nodedict[n])
    #for n in T.traverse_leaves():
    #    if n.label == "TTACCGCAGCAAATCA-1": #AACCATGGTAATGCGG-1":
    #        print("Parent of TTACCGCAGCAAATCA-1:", nodedict[n.get_parent()])
    if norm:
        return score / num_internal
    else:
        return score



def main(args):
    tree = read_tree_newick(args.tree)

    #print(m)
    if args.startle:
        print(pars_score_startle(args))
    else:
        print(pars_score(tree, args.msa, args.mchar, args.norm, args.prior))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--delimiter", type=str,
                        help="Input character matrix delimiter.",
                        default="\t",
                        required=True)
    parser.add_argument("-t1", "--tree", type=str,
                        help="Input file containing tree.",
                        required=True)
    parser.add_argument("-msa", type=str,
                        help="Input file containing character matrix.",
                        required=True)
    parser.add_argument("-mchar", type=str,
                        help="Missing character.",
                        required=True)
    parser.add_argument("--norm", action="store_true",
                        help="Whether to compute the normalized parsimony score.")
    parser.add_argument("--startle", action="store_true",
                        help="Whether to compute the parsimony score using startle methods.")
    parser.add_argument("--prior", type=str,
                        required=False,
                        default="",
                        help="Prior file to help compute weighted parsimony score.")
    main(parser.parse_args())

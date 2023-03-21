#! /usr/bin/env python
from problin_libs.sim_lib import *
from treeswift import *
import random
import argparse

def scale_tree(nwk_str,depth=1.0):
    # scale the input tree to the specified depth
    # if the input tree is not ultrametric, 
    # scale the average root-to-tip to the specified depth
    tree = read_tree_newick(nwk_str)
    tree.root.depth = 0
    for node in tree.traverse_preorder():
        if not node.is_root():
            node.depth = node.parent.depth + node.edge_length
    leaf_depths = [node.depth for node in tree.traverse_leaves()]
    avg_depth = sum(leaf_depths)/len(leaf_depths)
    mu = depth/avg_depth
    for node in tree.traverse_preorder():
        if not node.is_root():
            node.edge_length *= mu
            node.depth = node.parent.depth + node.edge_length
    return tree.newick()           

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-n", "--numcells",type=int,required=False, help="The number of cells (leaves) in the simulated tree.")
    parser.add_argument("-r","--reps",type=int,default=1,required=False,help="The number of replicates to create.")
    parser.add_argument("-d","--depth",type=float,default=1.0,required=False,help="The tree depth. Default: 1.0.")
    parser.add_argument("--randseed",required=False,type=int,help="Random seed: an interger number.")
    parser.add_argument("-o","--outprefix",required=True,help="The prefix of output files.")

    args = vars(parser.parse_args())

    if args["randseed"] is not None:
        random.seed(args["randseed"])

    nreps = args["reps"]
    nleaves = args["numcells"]

    for i in range(nreps):
        nwk_str = simTree_lnorm(nleaves,1,0.1)
        scaled_tree = scale_tree(nwk_str,depth=args["depth"])
        outfile = args["outprefix"] + "_r" + str(i+1).rjust(len(str(nreps)),'0') + ".nwk"
        with open(outfile,'w') as fout:
            fout.write(scaled_tree)

if __name__ == "__main__":
    main()            

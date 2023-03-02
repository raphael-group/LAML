#! /usr/bin/env python

from treeswift import *
from random import sample

tree = read_tree_newick("3724_NT_T1_tree_suppressed_unifurcations.nwk")

prune_set = []
for node in tree.traverse_postorder():
    C1 = [c for c in node.children if c.is_leaf()]
    C2 = [c for c in node.children if not c.is_leaf()]
    if len(C1) > 2 or (len(C1) >= 2 and len(C2) > 0):
        if len(C2) > 0: # only keep 1 child in C1
            prune_set += [c.label for c in sample(C1,len(C1)-1)]
        else: # keep 2 of the children in C1
            prune_set += [c.label for c in sample(C1,len(C1)-2)]

tree_pruned = tree.extract_tree_without(prune_set)            
tree_pruned.write_tree_newick("temp.nwk")
            

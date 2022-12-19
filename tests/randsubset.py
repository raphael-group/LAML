#! /usr/bin/env python

from treeswift import *
from random import sample

n = 100

tree = read_tree_newick("3724_NT_T1_tree_suppressed_unifurcations.nwk")

select_set = sample([node.label for node in tree.traverse_leaves()],n)

tree_pruned = tree.extract_tree_with(select_set)            
tree_pruned.write_tree_newick("rand"+ str(n) + ".nwk")      

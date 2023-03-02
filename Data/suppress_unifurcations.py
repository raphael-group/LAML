#! /usr/bin/env python

from treeswift import *

tree = read_tree_newick("3724_NT_T1_tree.nwk")
tree.suppress_unifurcations()
tree.write_tree_newick("3724_NT_T1_tree_suppressed_unifurcations.nwk")

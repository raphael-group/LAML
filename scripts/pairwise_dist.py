#! /usr/bin/env python

from sys import argv
from treeswift import *

inTree = read_tree_newick("trial_simulated_tree.tree")

mu = float(argv[1])
for node in inTree.traverse_postorder():
	if node.is_leaf():
		node.leaves = set([(node.label,0)])
	else:
		if len(node.children) < 2:
			continue
		c1,c2 = node.children
		for x,dx in c1.leaves:
			for y,dy in c2.leaves:				
				if x < y:
					a,b = x,y
					da,db = dx+c1.edge_length*mu,dy+c2.edge_length*mu
				else:
					a,b = y,x
					da,db = dy+c2.edge_length*mu,dx+c1.edge_length*mu		
				print(a + "_" + b + " " + str(da) + " " + str(db))
		node.leaves = set([])
		for x,dx in c1.leaves:
			node.leaves.add((x,dx+c1.edge_length*mu))
		for y,dy in c2.leaves:
			node.leaves.add((y,dy+c2.edge_length*mu))

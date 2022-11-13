#! /usr/bin/env python

from treeswift import *
from math import *
from random import random
from sys import argv

k = int(argv[1])
mu = 0.025
nstate=9

tree = read_tree_newick("trial_simulated_tree.tree")
tree.root.seq = '0'*k

for node in tree.traverse_preorder():
	if node.is_root():
		continue
	p = 1-exp(-mu*node.edge_length)
	pnode = node.parent
	seq = ''
	for c in pnode.seq:
		if c != '0':
			seq += c
		else:
			r = random()
			if r < p: # change state
				seq += str(int(random()*nstate) + 1) 
			else: # add '0'
				seq += '0'  
	node.seq = seq
	#if node.is_leaf():
	#	print(node.label +  " " + node.seq)	

L = list(tree.traverse_leaves())
n = len(L)
for i in range(n-1):
	for j in range(i+1,n):
		a,sa = L[i].label,L[i].seq	
		b,sb = L[j].label,L[j].seq	
		za = 0
		zb = 0
		z = 0
		for ca,cb in zip(sa,sb):
			za += (ca == '0')	
			zb += (cb == '0')
			z += (ca == '0' or cb == '0' or ca != cb)	
		da = -log(za/z)	
		db = -log(zb/z)
		if a > b:
			a,b = b,a
			da,db = db,da
		print(a + "_" + b + " " + str(da) + " " + str(db))	

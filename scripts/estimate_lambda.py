from treeswift import *
from math import *

h = 15 # generations

tree = read_tree_newick("GSE146712_lg10_tree_hybrid.processed.nwk")
n = 0
z = 0
states = []

for node in tree.traverse_leaves():
	chars = node.parent.label.split("|")
	if len(states) == 0:
		states = [set([x]) for x in chars]
	else:
		for c,site in zip(chars,states):
			site.add(c)	
	n += len(chars)
	z += len([x for x in chars if x == '0'])	

print(-log((n-z)/n)/h)
#for site in states:
#	print(len(site))

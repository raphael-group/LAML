import random 
import treeswift

random.seed(1984)

def pick_child(n, ndict):
	if n.is_leaf():
		# print("n is leaf", n)
		return n
	else:
		# print("n is internal")
		children = n.child_nodes()
		x = random.random()
		# print(x, ndict[n], [y.label for y in children])
		child = children[0] if x < ndict[n][0] else children[1]
		# print("picking child", child)
		return pick_child(child, ndict)

def subsample(T, num_nodes):
	"""
	Takes as input a newick tree file T and a number of nodes n. 
	n must be less than or equal to the number of nodes in tre T.
	Returns a new treeswift object 'newtree' containing exactly n nodes.
	"""
	# make treeswift tree object
	T = treeswift.read_tree_newick(T)
	print("Input tree:", T.newick())	
	ndict = dict()
	for n in T.traverse_postorder():
		if not n.is_leaf():
			# print(n.child_nodes())
			c1, c2 = n.child_nodes()
			c1_num = float(len(list(c1.traverse_leaves())))
			# print(c1_num)
			c2_num = float(len(list(c2.traverse_leaves())))
			# print(c2_num)	
			total = float(c1_num + c2_num)
			ndict[n] = [c1_num/total, c2_num/total]
	
	#print(T.num_nodes(leaves=True, internal=False))
	while T.num_nodes(leaves=True, internal=False) > num_nodes:
		# pick a child to go to and remove_child
		n = pick_child(T.root, ndict)
		# print("child to remove", n)
		parent = n.get_parent()
		# rm nodes
		parent.remove_child(n)
		# print(T.num_nodes(leaves=True, internal=False))
	return T

if __name__ == "__main__":
	filename = ""
	T = "((a:1, b:1)1:1, c:1)2:1;" # filename or string
	n = 2
	print("Output tree:", subsample(T, n))

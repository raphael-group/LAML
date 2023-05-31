from treeswift import *

tree = read_tree_newick("startle_random_resolve.nwk")
tree.collapse_short_branches(0.0055)
tree.write_tree_newick("startle_polytomies.nwk")

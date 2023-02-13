from problin_libs.EM_solver import EM_solver
from problin_libs import conv_eps
import treeswift
import random
#https://github.com/raphael-group/startle/blob/main/src/nni/startle.py

# Input 1
#t = "((a:1,b:1)d:1,c:1)e:1;"
#Q = [{0:0, 1:0.5, 2:0.5}, {0:0, 1:0.5, 2:0.5}] #, 2:0.5}]
#msa = {'a':[0, 0], 'b':[1, 1], 'c':[1, 1]}

# Input 2
t = "((a:1,b:1)e:2,(c:1,d:1)f:2)g:1;"
Q = [{0:0, 1:0.5, 2:0.5}, {0:0, 1:0.5, 2:0.5}] #, 2:0.5}]
msa = {'a':[0, 0], 'b':[1, 1], 'c':[1, 2], 'd':[1, 2]}

def apply_nni(currtree, instr): #parenta, childa, parentb, childb):
#def apply_nni(currtree, parenta, childa, parentb, childb):
    edge1, edge2 = instr['move']
    parenta, childa = edge1
    parentb, childb= edge2
    
    label_to_node = currtree.label_to_node(selection="all")
    # move childa from parenta to parentb
    label_to_node[childa].set_parent(label_to_node[parentb])
    label_to_node[parenta].remove_child(label_to_node[childa])
    label_to_node[parentb].add_child(label_to_node[childa])
    # move childb from parentb to parenta
    label_to_node[childb].set_parent(label_to_node[parenta])
    label_to_node[parentb].remove_child(label_to_node[childb])
    label_to_node[parenta].add_child(label_to_node[childb])
    return currtree
            
def move_set_treecopy(nt, parenta, parentb, childa):
    label_to_node = nt.label_to_node(selection="all")
    label_to_node[childa].set_parent(label_to_node[parentb])
    label_to_node[parenta].remove_child(label_to_node[childa])
    label_to_node[parentb].add_child(label_to_node[childa])

def move_set(parenta, parentb, childa):
    childa.set_parent(parentb)
    parenta.remove_child(childa)
    parentb.add_child(childa)

def undo_move_set(parenta, parentb, childa):
    childa.set_parent(parenta)
    parenta.add_child(childa)
    parentb.remove_child(childa)

def score_tree(msa, Q, tree):
    mySolver = EM_solver(msa, Q, tree.newick())
    optimal_llh = mySolver.optimize(initials=10, fixed_phi=1e-10, fixed_nu=1e-10, verbose=False)
    #print("Optimal tree: " + mySolver.params.tree.newick() + "\n")
    #print("Optimal negative-llh: " + str(optimal_llh) + "\n")
    #print("Optimal dropout rate: " + str(mySolver.params.phi) + "\n")
    #print("Optimal silencing rate: " + str(mySolver.params.nu) + "\n")
    return mySolver.params.tree.newick(), optimal_llh

def get_nnis(u):
    v = u.get_parent()
    print(v)
    u_edges = [(u, w) for w in u.child_nodes()]
    v_edges = [(v, w) for w in v.child_nodes() if w is not u]
    nni_moves = []
    for (u, w) in u_edges:
        for (v, z) in v_edges:
            #nni_moves.append({"move": ((u, w), (v, z))})
            nni_moves.append({"move": ((u.label, w.label), (v.label, z.label))})
    return nni_moves

def stochastic_branch(t):
    num_nodes = t.num_nodes(internal=True, leaves=False)
    # TODO: Random seed
    rn = random.randint(1, num_nodes - 2) # don't pick the root

    nodes = []
    for node in t.traverse_inorder(internal=True, leaves=False):
        if not node.is_root():
            nodes.append(node)
    return nodes[rn] # TODO: return node label

def longest_branch(t): 
    # find longest edge e = (u, v), incident to node u
    edgelens = dict()
    for node in tree.traverse_inorder(internal=True, leaves=False):
        edgelens[node] = node.get_edge_length()
        
    u = max(edgelens, key=edgelens.get)
    maxval = edgelens[u]
    #us = [u for u in edgelens if edgelens[u] == maxval]
    return u

def inplace_nni(tree, u): 
    # this changes the tree object inside! dangerous, but constant time 
    # tree object and branch incident to node object u
    # TODO: Do the NNI, after this, output the newick string
    # TODO: Revert the tree to original state!
    pass

def nni(treestr):
#`def nni(msa, Q, treestr):
    tree = treeswift.read_tree_newick(treestr)
    u = stochastic_branch(tree)
    nni_moves = get_nnis(u)
    v = u.get_parent()
    out_strs = []

    for i, instr in enumerate(nni_moves):
        # make a copy of the tree
        curr_tree = tree.extract_subtree(tree.root)

        curr_tree = apply_nni(curr_tree, instr) #parenta, cladea, parentb, cladeb)
        out_strs.append(curr_tree.newick())
    print(out_strs)

def main():
    nni(t)
# TODO: Set up argument flags parsing

def main2(): 
    tree = treeswift.read_tree_newick(t)
    opt_tree, pre_llh = score_tree(msa, Q, tree)
    em_iter, max_iter = 1, 1
    while 1:
        curr_tree, curr_llh = score_tree(msa, Q, tree) 
        #print([u.label for u in us])
        topo_dict = {}

        # pick random internal branch
        print(tree)
        u = stochastic_branch(tree)
        print("chosen edge is incident to:", u) #, str(u.edge_length()))

        # set up all nni moves
        nni_moves = get_nnis(u)
        v = u.get_parent()

        opt_tree, opt_score = score_tree(msa, Q, tree)
        topo_dict[opt_tree] = opt_score

        print(len(nni_moves), "nni moves")

        # produce new topology
        for i, instr in enumerate(nni_moves):
            print("NNI idx ", str(i))

            # make a copy of the tree
            curr_tree = tree.extract_subtree(tree.root)
            #nt = tree.extract_tree(None, False, False)

            edge1, edge2 = instr['move']
            parenta, cladea = edge1
            parentb, cladeb = edge2

            curr_tree = apply_nni(curr_tree, parenta, cladea, parentb, cladeb)
            #move_set(curr_tree, parenta, parentb, cladea)
            #move_set(curr_tree, parentb, parenta, cladeb)
            print("new tree:", curr_tree)

            # call functions to score each move
            opt_tree, opt_score = score_tree(msa, Q, curr_tree)
            topo_dict[opt_tree] = opt_score

            #print(cladea.get_parent() == parenta)
            #print(cladeb.get_parent() == parentb)

            # TODO: alternatively, use label_to_node to deal with copies of trees
            #undo_move_set(parenta, parentb, cladea)
            #print(cladea.get_parent(), parenta)
            #undo_move_set(parentb, parenta, cladeb)
            #print(cladeb.get_parent(), parentb)
            #print("reset tree:", tree)

        for topo in topo_dict:
            print(topo, topo_dict[topo])

        # TODO: accept the new tree topology, and keep searching
        
        if curr_llh - pre_llh < conv_eps or em_iter < max_iter:
            break
        pre_llh = curr_llh
        em_iter += 1



if __name__ == "__main__":
    main()

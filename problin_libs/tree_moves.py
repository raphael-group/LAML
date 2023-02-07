from problin_libs.EM_solver import EM_solver
import treeswift
#https://github.com/raphael-group/startle/blob/main/src/nni/startle.py

# Input 1
#t = "((a:1,b:1)d:1,c:1)e:1;"
#Q = [{0:0, 1:0.5, 2:0.5}, {0:0, 1:0.5, 2:0.5}] #, 2:0.5}]
#msa = {'a':[0, 0], 'b':[1, 1], 'c':[1, 1]}

# Input 2
t = "((a:1,b:1)e:2,(c:1,d:1)f:2)g:1;"
Q = [{0:0, 1:0.5, 2:0.5}, {0:0, 1:0.5, 2:0.5}] #, 2:0.5}]
msa = {'a':[0, 0], 'b':[1, 1], 'c':[1, 2], 'd':[1, 2]}

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
    u_edges = [(u, w) for w in u.child_nodes()]
    v_edges = [(v, w) for w in v.child_nodes() if w is not u]
    nni_moves = []
    for (u, w) in u_edges:
        for (v, z) in v_edges:
            nni_moves.append({"move": ((u, w), (v, z))})
            # nni_moves.append({"move": ((u.label, w.label), (v.label, z.label))}
    return nni_moves

def main(): 
    tree = treeswift.read_tree_newick(t)

    # find longest edge e = (u, v), incident to node u
    edgelens = dict()
    for node in tree.traverse_inorder(internal=True, leaves=False):
        edgelens[node] = node.get_edge_length()
        
    u = max(edgelens, key=edgelens.get)
    maxval = edgelens[u]
    #us = [u for u in edgelens if edgelens[u] == maxval]

    #print([u.label for u in us])
    topo_dict = {}

    #for u in us:
    print(tree)
    print("max edge is incident to:", u, edgelens[u])

    # set up all nni moves
    nni_moves = get_nnis(u)
    v = u.get_parent()

    opt_tree, opt_score = score_tree(msa, Q, tree)
    topo_dict[opt_tree] = opt_score

    print(len(nni_moves), "nni moves")

    # produce new topology
    for i, instr in enumerate(nni_moves):
        print("NNI idx ", str(i))
        #nt = tree.extract_tree(None, False, False)
        #print("original tree:", tree)
        # print(instr['move'])

        edge1, edge2 = instr['move']
        parenta, cladea = edge1
        #print(parenta.label, cladea.label)
        parentb, cladeb = edge2
        #print(parentb.label, cladeb.label)

        move_set(parenta, parentb, cladea)
        move_set(parentb, parenta, cladeb)
        
        print("new tree:", tree)

        # call functions to score each move
        opt_tree, opt_score = score_tree(msa, Q, tree)
        topo_dict[opt_tree] = opt_score

        #print(cladea.get_parent() == parenta)
        #print(cladeb.get_parent() == parentb)

        # alternatively, use label_to_node to deal with copies of trees
        undo_move_set(parenta, parentb, cladea)
        #print(cladea.get_parent(), parenta)
        undo_move_set(parentb, parenta, cladeb)
        #print(cladeb.get_parent(), parentb)

        print("reset tree:", tree)

    for topo in topo_dict:
        print(topo, topo_dict[topo])

if __name__ == "__main__":
    main()

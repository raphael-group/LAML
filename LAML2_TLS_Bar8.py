#!/bin/python
from laml_libs.IO_handler.DLT_parser import DLT_parser
from laml_libs.Count_model.PMMC_model import PMMC_model
import laml_libs as laml
from laml_libs.TopoSearch.Topology_search import Topology_search
from copy import deepcopy

datafile = "examples/example5/TLS_Bar8.json" 
starting_tree = "examples/example5/Bar8_newick_noMutationlessEdges_Labeled.nwk"
prefix = "output_example5/LAMLPro_TLS_Bar8_opt2"

resolve_polytomies = True

parser = DLT_parser()
parser.get_from_path(datafile=datafile, priorfile=None, max_allele_per_cassette=3)

selected_model = PMMC_model

data = {'DLT_data': parser.DLT_data }
prior = {'Q': parser.priors}

with open(starting_tree, "r") as f:
    input_trees = []
    for line in f:
        input_trees.append(line.strip())

fixed_params = {"nu": 0}
myTopoSearch = Topology_search(input_trees, selected_model, data=data, prior=prior, params=fixed_params)
my_strategy = deepcopy(laml.DEFAULT_STRATEGY)

my_strategy['ultra_constr'] = True
my_strategy['resolve_search_only'] = False
my_strategy['fixed_params'] = deepcopy(fixed_params)


# Full topology search
opt_trees,max_score,opt_params = myTopoSearch.search(resolve_polytomies=resolve_polytomies,maxiter=500,verbose=True, strategy=my_strategy, nreps=1,checkpoint_file=f"{prefix}_ckpt.txt")
nllh = -max_score

# post-processing: analyze results and output
out_tree = prefix + "_trees.nwk"
out_params = prefix + "_params.txt"

with open(out_tree,'w') as fout:
    for tstr in opt_trees:
        tree = read_tree_newick(tstr)
        # get the height of the tree
        tree_height = tree.height(weighted=True) # includes the root's length, edit units
        # get the edit rate
        if 'lambda' in fixed_params:
            edit_rate = fixed_params['lambda']
            print(f"Tree height pre-scaling: {tree_height}, input edit rate: {edit_rate}")
        else:
            edit_rate = tree_height/float(args['timescale'])
            print(f"Tree height pre-scaling: {tree_height}, input timescale: {args['timescale']}")
        for node in tree.traverse_preorder():
            if node.edge_length:
                node.edge_length = node.edge_length / edit_rate
        tree_height = tree.height(weighted=True)
        print(f"Tree height after scaling: {tree_height}")
        if len(tree.root.children) > 1:
            new_root = Node()
            new_root.add_child(tree.root)
            tree.root = new_root
        fout.write(tree.newick() + "\n")

with open(out_params,'w') as fout:
    fout.write("Dropout probability: phi=" + str(opt_params['phi']) + "\n")
    fout.write("Silencing rate: nu=" + str(opt_params['nu']) + "\n")
    fout.write("Sequencing accuracy: rho=" + str(opt_params['rho']) + "\n")
    fout.write("Mutation rate: lambda=" +  str(edit_rate) + "\n")
    fout.write("Negative-llh: " +  str(nllh) + "\n")

stop_time = timeit.default_timer()
print("Runtime (s):", stop_time - start_time)



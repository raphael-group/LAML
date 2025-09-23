import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from laml_libs import run_laml_infer

# remove all generated files
out = run_laml_infer(character_matrix="character_matrix.csv",
                   tree_topology="starting.tree", 
                   output_prefix="example2_LAMLapi",
                   priors="priors.csv",
                   solver="fastEM-cpu",
                   nInitials=1,
                   randomreps=1,
                   topology_search=True,
                   timescale=10,
                   verbose=True
                   )
print(out.params)
print(out.tree_newick[:120], "...")

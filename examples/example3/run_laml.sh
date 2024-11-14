python ../../run_laml2.py -c c30_noise0.02_rep84.msa.txt -t c30_noise0.02_rep84.nj_tree.nwk -p c30_noise0.02_rep84.prior.csv -o example3 --nInitials 1 --randomreps 1 -v --topology_search --readout_model "PMMN"

out1=$(python /Users/gc3045/useful_scripts/compare_two_trees.py -t1 c30_noise0.02_rep84.nj_tree.nwk -t2 c30_noise0.02_rep84.true_tree.nwk)
echo "Compare NJ and true tree: ${out1}"
out2=$(python /Users/gc3045/useful_scripts/compare_two_trees.py -t1 example3_trees.nwk -t2 c30_noise0.02_rep84.true_tree.nwk)
echo "Compare LAML and true tree: ${out2}"

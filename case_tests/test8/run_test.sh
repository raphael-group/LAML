#! /bin/bash

python ../../run_problin.py -c character_matrix.csv -t tree.nwk -p prior_k30_r01.csv  --delimiter comma -o problin_output_ultr.txt --nInitials 1 --randomreps 1 --topology_search -v --ultrametric > problin_output_ultr.log 2>&1

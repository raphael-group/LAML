#! /bin/bash

python ../../run_problin.py -c n64_d0s32_r1_character_matrix.csv -t n64.tre -p n64_d0s32_priors.csv --delimiter comma -o n64_problin.txt --nInitials 1 --randomreps 1 --topology_search -v --ultrametric

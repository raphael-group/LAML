#! /bin/bash

rm -f problin_output.txt
python ../../run_problin.py -c bar10_character_matrix_without_cluster.csv -t t3_r4_stitch_polytomies_optimal_llh_collapsed_0.007.nwk -p bar10_mutation_prior.csv --delimiter comma -o problin_output.txt --nInitials 1 -v

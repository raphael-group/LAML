#! /bin/bash

rm -f n64_problin.txt
python ../../run_problin.py -c n64_d0_r1_character_matrix.csv -t n64.tre -p n64_d0_priors.csv --delimiter comma -o n64_problin --nInitials 1

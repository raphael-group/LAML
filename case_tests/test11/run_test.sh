#! /bin/bash

python ../../run_problin.py -c n64_d0s32_r1_character_matrix.csv -t n64.tre -p n64_d0s32_priors.csv --delimiter comma -o n64_problin.txt -L "0 0.32" --placement --solver Scipy

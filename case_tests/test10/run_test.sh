#! /bin/bash

rm -f problin_output.txt
python ../../run_problin.py -c 3724_NT_All_dedup_mtx.txt -t startle_random_resolve.nwk  -p 3724_NT_All_priors_norm.pkl --delimiter comma -o problin_output.txt --nInitials 1 -v --ultrametric

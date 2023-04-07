#! /bin/bash

rm -f n65_problin.txt
python ../../run_problin.py -c n65_cmtx.csv -t n65.tre -p n65_priors.csv --delimiter comma -o n65_problin.txt --nInitials 1 --randomreps 1 

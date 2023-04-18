#! /bin/bash

rm -f problin_out.txt
python ../../run_problin.py -c cmtx.csv -t test.tre -p priors.csv --delimiter comma -o problin_out.txt --nInitials 1 --randomreps 1 --likelihood 
#python ../../run_problin.py -c cmtx.csv -t test.tre -p priors.csv --delimiter comma -o problin_out.txt --nInitials 1 --randomreps 1 --likelihood --phi 0.0000000001 --nu 0.0000000001 

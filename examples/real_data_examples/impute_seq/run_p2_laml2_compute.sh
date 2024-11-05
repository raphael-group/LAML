#!/bin/bash

topology="/Users/gc3045/laml2_experiments/proc_realdata/baseMemoir/improved_lamlpro_trees/LAML2_p2_rep12_trees.nwk"
input_observations="/Users/gc3045/laml2_experiments/proc_realdata/baseMemoir/baseMemoir.msa.txt"
emission_file="/Users/gc3045/laml2_experiments/proc_realdata/baseMemoir/baseMemoir.emissions.txt"
outputfile="/Users/gc3045/laml2_experiments/proc_realdata/baseMemoir/improved_lamlpro_trees/LAML2_p2_rep12_computellh"

python /Users/gc3045/scmail_v1/LAML/run_laml2.py -c $input_observations -t $topology -p "uniform" -m -1 -y "character_matrix" -M "PMMN" -o ${outputfile} -v --noSilence --e $emission_file --topology_search --compute_llh --fixedParams "lambda=0.13222946088246448 phi=0.29949494949494687 nu=0 rho=0.9878836840443251"

#Dropout probability: phi=0.29949494949494687
#Silencing rate: nu=0
#Sequencing accuracy: rho=0.9878836840443251
#Mutation rate: lambda=0.13222946088246448
#Negative-llh: 14992.479109028593

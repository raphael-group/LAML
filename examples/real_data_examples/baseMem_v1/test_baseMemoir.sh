#!/bin/bash
#

#input_observations="baseMemoir_test_small.json"
#input_observations="baseMemoir_test.json"
topology="baseMem_random_topo.nwk"
# priors are uniform
#outputfile="LAML2_test_baseMem"

# test 1
#echo "python /Users/gc3045/scmail_v1/LAML/run_laml2.py -c $input_observations -t $topology -p "uniform" -m -1 -y "observed_features" -M "PMMN" -o ${outputfile} -v --noSilence" #--topology_search" 
#python /Users/gc3045/scmail_v1/LAML/run_laml2.py -c $input_observations -t $topology -p "uniform" -m -1 -y "observed_features" -M "PMMN" -o ${outputfile} -v --noSilence # --topology_search

input_observations="baseMemoir.msa.small.txt"
emission_file="baseMemoir.emissions.small.txt"

# test 2
#python /Users/gc3045/scmail_v1/LAML/run_laml2.py -c $input_observations -t $topology -p "uniform" -m -1 -y "character_matrix" -M "PMMN" -o ${outputfile} -v --noSilence --e $emission_file # --topology_search

# test 3
outputfile="test_3/LAML2_test_baseMem"
python /Users/gc3045/scmail_v1/LAML/run_laml2.py -c $input_observations -t $topology -p "uniform" -m -1 -y "character_matrix" -M "PMMN" -o ${outputfile} -v --noSilence --e $emission_file --topology_search


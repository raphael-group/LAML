#!/bin/bash
#

input_observations="baseMemoir_test.json"
topology="baseMem_random_topo.nwk"
# priors are uniform
outputfile="LAML2_test_baseMem"

echo "python /Users/gc3045/scmail_v1/LAML/run_laml2.py -c $input_observations -t $topology -p "uniform" -m -1 -y "observed_features" -M "PMMN" -o ${outputfile} -v --noSilence" 
python /Users/gc3045/scmail_v1/LAML/run_laml2.py -c $input_observations -t $topology -p "uniform" -m -1 -y "observed_features" -M "PMMN" -o ${outputfile} -v --noSilence # --topology_search

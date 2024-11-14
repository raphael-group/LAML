#!/bin/bash

python ../../run_laml2.py -c TLS_Bar8.json -t Bar8_newick_noMutationlessEdges_Labeled.nwk -p priors_Bar8.csv -m -1 --readout_model "PMMC" -o LAML2_Bar8_test_wpriors -v -s 3 --topology_search 

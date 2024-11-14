#!/bin/bash

python ../../run_laml2.py -c TLS_Bar8.json -t Bar8_newick_noMutationlessEdges_Labeled.nwk -p uniform -m -1 --readout_model "PMMC" -o output_example5/LAML2_Bar8_test -v -s 3 --topology_search 

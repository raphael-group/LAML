#! /bin/bash

for x in k*_rep*; do python MLP_brlen.py model_tree.nwk $x/characters.txt | sort | sed -e "s/^/$x /g" | sed -e "s/k//g" -e "s/_/ /g"; done > results_MLP_brlen.txt

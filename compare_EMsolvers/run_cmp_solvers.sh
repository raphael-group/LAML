#!/bin/bash

for i in $(seq 1 10); do

    python run_laml.py -c examples/example2/character_matrix.csv -t examples/example2/starting.tree -p examples/example2/priors.csv --nInitials 1 --randomreps 1 --topology_search -v --timescale 10 --solver fastEM -o example2_fastEM_rep${i}
    python run_laml.py -c examples/example2/character_matrix.csv -t examples/example2/starting.tree -p examples/example2/priors.csv --nInitials 1 --randomreps 1 --topology_search -v --timescale 10 --solver EM -o example2_oldEM_rep${i}

done

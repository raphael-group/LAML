#! /bin/bash

#rm n64_problin.txt n64_problin_rf.txt n64_problin.tre
#python ../../run_problin.py -c n64_d0s32_r1_character_matrix.csv -t n64.tre -p n64_d0s32_priors.csv --delimiter comma -o n64_problin.txt --nInitials 1 --randomreps 1 --topology_search -v --ultrametric
grep -A1 "Newick" n64_problin.txt | grep -v "Newick" > n64_problin.tre
java -jar ~/Packages_N_Libraries/TreeCmp_v2.0-b76/bin/treeCmp.jar -r n64_resolved.tre -d rf -i n64_problin.tre -o n64_problin_rf.txt
cat n64_problin_rf.txt

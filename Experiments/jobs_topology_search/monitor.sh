#!/bin/bash

# for k in 50 100 200 
# for k in 20 30 40 50 100 200 
for k in 20 30 40 50 100 200 300 400 500 5000
# for k in 200 300 400 500 
do
	num=$(ls mlpars_results_k$k/ | wc -l)
	echo "k: $k, num done: $num"
	#python find_best_topologies.py $k 1000
	#python find_best_mlpars_topologies.py $k 1000
done

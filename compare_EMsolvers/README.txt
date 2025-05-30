##########################################################################
CPU
##########################################################################
RF dist: from a few samples, fastEM is producing trees further from the old output, but with better likelihood. I don't have the true tree for this one since it was subsampled down to 64 leaves. Repeat this experiment on the server.

####### both_scores ########
grep "llh diff between EMsolvers (old - new):" /Users/gc3045/scmail_v1/LAML/example2_bothscores.log | awk '{sum += $NF; count++} END {print "Average llh diff:", sum/count}'

(base) gc3045@COS-T0264KMR6R compare_EMsolvers % ./compare_both.sh
Average llh diff: -6.47996
Average time diff: 9.13679
Average phi diff: -1.87903e-05
Average nu diff: 6.72325e-05

######## LLH ############
example2_fastEM_params.txt:Negative-llh: 1127.6241169964358
example2_fastEM_rep1_params.txt:Negative-llh: 1125.1478877421807
example2_fastEM_rep2_params.txt:Negative-llh: 1125.1470035015218
example2_oldEM_params.txt:Negative-llh: 1125.1501782313612
example2_oldEM_rep1_params.txt:Negative-llh: 1127.6217602830227

######## RF #############
(base) gc3045@COS-T0264KMR6R useful_scripts % python compare_two_trees.py -t1 /Users/gc3045/scmail_v1/LAML/examples/example2/example2_trees.nwk -t2 /Users/gc3045/scmail_v1/LAML/compare_EMsolvers/example2_fastEM_rep2_trees.nwk
64,62,62,5,5,0.080645
(base) gc3045@COS-T0264KMR6R useful_scripts % python compare_two_trees.py -t1 /Users/gc3045/scmail_v1/LAML/examples/example2/example2_trees.nwk -t2 /Users/gc3045/scmail_v1/LAML/compare_EMsolvers/example2_oldEM_trees.nwk
64,62,62,4,4,0.064516

######## Collect and compare #######
 - The number of NNI iterations
 - The difference in final LLH
 - The runtime difference

(base) gc3045@COS-T0264KMR6R compare_EMsolvers % python collect_stats.py
File: /Users/gc3045/scmail_v1/LAML/compare_EMsolvers/example2_fastEM_rep2.log
Average runtime per tree: 0.344549 seconds
Average number of trees checked: 23.50
Optimal score: -1125.1470035015218

File: /Users/gc3045/scmail_v1/LAML/compare_EMsolvers/example2_oldEM.log
Average runtime per tree: 0.303759 seconds
Average number of trees checked: 32.08
Optimal score: -1125.1501782313612

File: /Users/gc3045/scmail_v1/LAML/compare_EMsolvers/example2_fastEM_rep1.log
Average runtime per tree: 0.373085 seconds
Average number of trees checked: 23.92
Optimal score: -1125.1478877421807

File: /Users/gc3045/scmail_v1/LAML/compare_EMsolvers/example2_oldEM_rep1.log
Average runtime per tree: 0.287914 seconds
Average number of trees checked: 20.61
Optimal score: -1127.6217602830227

File: /Users/gc3045/scmail_v1/LAML/compare_EMsolvers/example2_fastEM.log
Average runtime per tree: 0.350769 seconds
Average number of trees checked: 22.88
Optimal score: -1127.6241169964358

File: /Users/gc3045/scmail_v1/LAML/compare_EMsolvers/example2_oldEM_rep2.log
Average runtime per tree: 0.160119 seconds
Average number of trees checked: 26.96
Optimal score: None

##########################################################################
GPU
##########################################################################

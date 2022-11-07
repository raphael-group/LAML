#!/bin/bash
#SBATCH --job-name=topo_search	 # create a short name for your job
#SBATCH --output=./run_test/slurm-%A.%a.out # stdout file
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=2G         # memory per cpu-core (4G is default)
#SBATCH --time=10:00:00          # total run time limit (HH:MM:SS)



m=10
#for k in 20 30 40 50 100 200 
# not run yet: 100 200 
#do
k=100
outdir="/n/fs/ragr-research/projects/problin/jobs_topology_search/log_ml_results_k${k}"
mkdir -p ${outdir}

for rep in {0..99}
do 
	for idx in {0..14}
	do
		logrun="${outdir}/m${m}_k${k}_rep${rep}_topoidx${idx}.logrun"
		slurmout="${outdir}/m${m}_k${k}_rep${rep}_topoidx${idx}.logout"
		outfile="/n/fs/ragr-research/projects/problin/jobs_topology_search/ml_results_k${k}_rep${rep}/topo${idx}.txt"
		if [ ! -f ${outfile} ]
		then

			sbatch run_test_topologies.sh ${idx} ${k} ${rep} ${logrun}
		else
			echo "${outfile} already exists."
		fi
	done
done
#done

#!/bin/bash
#SBATCH --job-name=test_topo 	 # create a short name for your job
#SBATCH --output=./slurm_mlpars/slurm-%A.%a.out # stdout file
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=2G         # memory per cpu-core (4G is default)
#SBATCH --time=10:00:00          # total run time limit (HH:MM:SS)

# conda activate py3k

k=$1
logrun=$2
s=$3
e=$(($s+99))
rep=$s
while [ "$rep" -le "$e" ]; do
	for idx in {0..14}
	do
		# mkdir -p "/n/fs/ragr-research/projects/problin/jobs_topology_search/mlpars_results_k${k}/rep${rep}"
		outfile="/n/fs/ragr-research/projects/problin/jobs_topology_search/mlpars_results_k${k}/rep${rep}/topo${idx}.txt"
		if [ ! -f ${outfile} ]
		then
			start=`date +%s`
			echo "python test_topologies.py $idx $k $rep"
			python test_mlpars.py $idx $k $rep
			end=`date +%s`
			runtime=$( echo "$end - $start" | bc -l )
			echo "${runtime}" > ${logrun}
		else
			echo "${outfile} already exists."
		fi
			
	done 
	rep=$(($rep+1))
done

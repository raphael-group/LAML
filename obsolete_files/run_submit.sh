#!/bin/bash
#SBATCH --job-name=topo_search	 # create a short name for your job
#SBATCH --output=./run_log_bintree/slurm-%A.%a.out # stdout file
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=2G         # memory per cpu-core (4G is default)
#SBATCH --time=10:00:00          # total run time limit (HH:MM:SS)


mkdir -p run_log_bintree
outdir="/n/fs/ragr-research/projects/problin/problin_libs/log_bintree"
mkdir -p ${outdir} 

m=10
k=50 

# for n in 30 50 100 200 300 500 1000
for n in 50 100 200 300 500 1000
do
#n=30
	logrun="${outdir}/n${n}_k${k}_m${m}_bintre.logrun"
	outfile="/n/fs/ragr-research/projects/problin/problin_libs/n${n}_k${k}_m${m}_bintre.txt"
	if [ ! -f ${outfile} ]
	then

		echo "sbatch run_test_fels.sh ${k} ${m} ${n} ${logrun}"
		sbatch run_test_fels.sh ${k} ${m} ${n} ${logrun} ${outfile}
	else
		echo "${outfile} already exists."
	fi
done

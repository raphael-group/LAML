#!/bin/bash
#SBATCH --job-name=test_fels 	 # create a short name for your job
#SBATCH --output=run_test/slurm-%A.%a.out # stdout file
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=2G         # memory per cpu-core (4G is default)
#SBATCH --time=10:00:00          # total run time limit (HH:MM:SS)

# outdir="/n/fs/ragr-research/projects/problin/jobs/run_true_branches"
outdir="/n/fs/ragr-research/projects/problin/jobs/run_no_branches"

mkdir -p ${outdir}

k=5000
m=10
initials=15
#t="((a:0.0360971597765934,b:3.339535381892265)e:0.0360971597765934,(c:0.0360971597765934,d:3.339535381892265)f:0.0360971597765934)r:0.0;"
t="((a,b)e,(c,d)f)r;"

for i in {0..200}
do
	logfile="${outdir}/m${m}_k${k}_idx${i}.logfile"
	logrun="${outdir}/m${m}_k${k}_idx${i}.logrun"
	slurmout="${outdir}/m${m}_k${k}_idx${i}.logout"
	sbatch run_test_fels.sh ${k} ${m} ${i} ${initials} ${t} ${logfile} ${logrun} 

done


#!/bin/bash
#SBATCH --job-name=test_topo 	 # create a short name for your job
#SBATCH --output=./run_log_bintree/slurm-%A.%a.out # stdout file
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=2G         # memory per cpu-core (4G is default)
#SBATCH --time=10:00:00          # total run time limit (HH:MM:SS)

# conda activate py3k

k=$1
m=$2
n=$3
logrun=$4
logfile=$5

start=`date +%s`
echo "python run_test_fels.py $k $m $n"
python run_test_fels.py $k $m $n > ${logfile}
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "${runtime}" > ${logrun}

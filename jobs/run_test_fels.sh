#!/bin/bash
#SBATCH --job-name=test_fels 	 # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=2G         # memory per cpu-core (4G is default)
#SBATCH --time=10:00:00          # total run time limit (HH:MM:SS)

# conda activate py3k

k=$1 
m=$2
i=$3
initials=$4
t=$5
logfile=$6
logrun=$7

#./run_test_fels.sh ${k} ${m} ${i} ${initials} ${t} ${logfile} ${logrun}

start=`date +%s`
echo "python run_fels.py $k $m $i $initials $t > ${logfile}" > ${logfile} 
python run_fels.py $k $m $i $initials $t > ${logfile}
end=`date +%s`
runtime=$( echo "$end - $start" | bc -l )
echo "${runtime}" > ${logrun}

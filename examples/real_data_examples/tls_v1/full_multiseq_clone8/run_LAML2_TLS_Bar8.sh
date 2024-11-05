#!/bin/bash
#SBATCH --job-name=LAML2_Bar8_no_priors_no_top_search
#SBATCH --output=/Genomics/chanlab/richard/running_LAML2/slurm_out2/LAML2_Bar8_%A_%a.out
#SBATCH --error=/Genomics/chanlab/richard/running_LAML2/slurm_out2/LAML2_Bar8_%A_%a.err
#SBATCH --array=1-10
#SBATCH --mail-type=all
#SBATCH --mail-user=rz3427@princeton.edu
#SBATCH --mem=4G
#SBATCH --time=168:00:00 

module purge
module load conda
conda activate cell_fate_mapping

batch=$((SLURM_ARRAY_TASK_ID - 1))
seed=$(( 3291*batch + 381 ))
echo $seed

cd /Genomics/chanlab/richard/running_LAML2

python LAML2/run_laml2.py -c TLS_Bar8.json -t /Genomics/chanlab/richard/trees/Bar8_newick_noMutationlessEdges_Labeled.nwk -p uniform -m -1 -y allele_counts -M PMMC -o no_priors_no_top_search/LAML2_Bar8_$batch -v --randseeds $seed
# python LAML2/run_laml2.py -c TLS_Bar8.json -t /Genomics/chanlab/richard/trees/Bar8_newick_noMutationlessEdges_Labeled.nwk -p uniform -m -1 -y allele_counts -M PMMC -o no_priors_top_search/LAML2_Bar8_$batch -v --topology_search --randseeds $seed
# python LAML2/run_laml2.py -c TLS_Bar8.json -t /Genomics/chanlab/richard/trees/Bar8_newick_noMutationlessEdges_Labeled.nwk -p priors_Bar8.csv -m -1 -y allele_counts -M PMMC -o priors_no_top_search/LAML2_Bar8_$batch -v --randseeds $seed
# python LAML2/run_laml2.py -c TLS_Bar8.json -t /Genomics/chanlab/richard/trees/Bar8_newick_noMutationlessEdges_Labeled.nwk -p priors_Bar8.csv -m -1 -y allele_counts -M PMMC -o priors_top_search/LAML2_Bar8_$batch -v --randseeds $seed --topology_search
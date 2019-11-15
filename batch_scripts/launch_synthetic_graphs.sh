#!/bin/bash
#
#SBATCH --job-name=submit_brain
#SBATCH --output=/scratch/users/cdonnat/Bayesian_ICA/synthetic/logs/synth_experiment_%A.out
#SBATCH --error=/scratch/users/cdonnat/Bayesian_ICA/synthetic/logs/synth_experiment_%A.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --qos=normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cdonnat@stanford.edu
#SBATCH --output=R-brain_launch.log
#SBATCH --partition=hns,normal
# load modules
ml R
ml r-rstan



job_directory=/scratch/users/cdonnat/Bayesian_ICA/synthetic/logs
RDM_NB=$$
sig=$1 
prior=$2
for graph in {1..20}
do
           seed=${graph}
   	   name=DasMultiplication_synthetic_graph_exp${graph}_${RDM_NB}_${sig}_${prior}
	   echo "name: ${name}"
	   job_file="${job_directory}/${name}.job"
	   echo ${job_file}
	   if [ -f ${job_file} ]
	   then 
		  rm ${job_file}
	   fi
	   echo "#!/bin/bash
#
#SBATCH --job-name=submit_R_synthetic_experiment_multiple_graphs
#SBATCH --output=/scratch/users/cdonnat/Bayesian_ICA/synthetic/logs/R_synthetic_experiment_multiple_graphs${name}.out
#SBATCH --error=/scratch/users/cdonnat/Bayesian_ICA/synthetic/logs/R_synthetic_experiment_multiple_graphs${name}.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --qos=normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cdonnat@stanford.edu
#SBATCH --partition=hns,stat,normal,owners

# load modules
ml R
ml r-rstan

 
# execute script
cd $SCRATCH/Bayesian_ICA
FILENAME=${name}
echo ${name}
Rscript synthetic_experiments.R ${seed} ${name} ${sig} ${prior}
">>${job_file}
			sbatch ${job_file}

done

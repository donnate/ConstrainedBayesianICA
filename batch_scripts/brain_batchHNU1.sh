#!/bin/bash
#
#SBATCH --job-name=submit_brain
#SBATCH --output=/scratch/users/cdonnat/Bayesian_ICA/scripts_real_experiments/HNU1/logs/overall_launch_%A.out
#SBATCH --error=/scratch/users/cdonnat/Bayesian_ICA/scripts_real_experiments/HNU1/logs/overall_launch_%A.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --qos=normal
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cdonnat@stanford.edu
#SBATCH --partition=hns,stat,normal,owners
# load modules
ml R
ml r-rstan

job_directory=/scratch/users/cdonnat/Bayesian_ICA/scripts_real_experiments/HNU1/logs
results_dir=/scratch/users/cdonnat/Bayesian_ICA/scripts_real_experiments/HNU1/results
SCALE=$1
TRANSFORM_PRIOR=$2
TRANSFORM_TS=$3
ALPHA=$4
prior_lowerdiag=1
readheader=1


for subj in {25427..25456}
do
  echo $subj
  for session in {1..10}
  do    
    data_ts="/scratch/users/cdonnat/data/HNU1/ts_sub_${subj}_ses_${session}.csv"
    data_prior="/scratch/users/cdonnat/data/HNU1/dwi_graph_sub_${subj}_ses_${session}.csv"
    if [ ${SCALE} == 1 ]; then
        name="final_scaled_elementwise_mult_${subj}_run-${session}_tp${TRANSFORM_PRIOR}_tts${TRANSFORM_TS}_alpha_${ALPHA}"
    else
        name="final_unscaled_elementwise_mult_${subj}_run-${session}_tp${TRANSFORM_PRIOR}_tts${TRANSFORM_TS}_alpha_${ALPHA}"
    fi
    echo ${name}
  
    job_file="${job_directory}/${name}.job"
		       echo ${job_file}
                       if [ -f ${job_file} ]
                       then 
                          rm ${job_file}
			  rm "${job_directory}/${name}.err"
			  rm "${job_directory}/${name}.out"
                       fi
			   echo "#!/bin/bash
#
#SBATCH --job-name=/scratch/users/cdonnat/Bayesian_ICA/scripts_real_experiments/HNU1/logs/${name}.job
#SBATCH --output=/scratch/users/cdonnat/Bayesian_ICA/scripts_real_experiments/HNU1/logs/${name}.out
#SBATCH --error=/scratch/users/cdonnat/Bayesian_ICA/scripts_real_experiments/HNU1/logs/${name}.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=15GB
#SBATCH --mail-type=ALL
#SBATCH --mail-user=$USER@stanford.edu
#SBATCH --partition=hns,stat,owners

cd $SCRATCH/Bayesian_ICA
ml R
ml r-rstan
Rscript brain_ICA_analysis.R  ${subj} ${session} ${data_ts} ${data_prior} ${TRANSFORM_TS} ${SCALE} ${results_dir} ${prior_lowerdiag}  ${TRANSFORM_PRIOR} ${readheader} ${ALPHA}
">>${job_file}
			sbatch ${job_file}
   
  done
done

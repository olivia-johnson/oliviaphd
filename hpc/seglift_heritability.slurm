#!/bin/bash -l
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem=1G
#SBATCH --gres=tmpfs:1G
#SBATCH --mail-type=END                                         # Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL                                        # Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=olivia.johnson@adelaide.edu.au          # Email to which notifications will be sent


params=${1} #parameters, pass from text file outside slurm file (could use for loop to cycle through multiple sets outside of this script)
module load Anaconda3/2020.07
conda activate sim_env

# directories
param_dir="/hpcfs/users/a1704225/parameters"
p_name=$(basename -s ".txt" ${params});

results_dir="/hpcfs/users/a1704225/results/${p_name}"
mkdir -p ${results_dir}

#python script ${paramter file} ${array id} ${job id for tmpdir path in SLiM}
python heritability_run.py ${params} ${SLURM_ARRAY_TASK_ID} ${SLURM_JOB_ID}

#copy files from tmpdir to results folder
cp ${TMPDIR}/pedigree_${params}_${SLURM_ARRAY_TASK_ID}.txt ${results_dir}

conda deactivate

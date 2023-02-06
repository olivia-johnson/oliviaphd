#!/bin/bash -l
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=2G
#SBATCH --mail-type=END                                         # Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL                                        # Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=olivia.johnson@adelaide.edu.au          # Email to which notifications will be sent


sim_type=${1}
params=${2} #parameters, pass from text file outside slurm file (could use for loop to cycle through multiple sets outside of this script)
module load Anaconda3/2020.07
conda activate sim_env

# directories
param_dir="/hpcfs/users/a1704225/parameters/${sim_type}"
p_name=$(basename -s ".txt" ${params});

results_dir="/hpcfs/users/a1704225/results/${sim_type}/${p_name}"

#python script ${paramter file} ${array id} ${job id for tmpdir path in SLiM}
python analysis_run.py ${params} ${SLURM_ARRAY_TASK_ID} ${results_dir} ${sim_type}

conda deactivate

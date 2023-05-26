#!/bin/bash -l
#SBATCH --array=1
#SBATCH -p skylake,icelake,a100cpu,v100cpu
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=15:00:00
#SBATCH --mem=45G
#SBATCH --mail-type=END                                         # Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL                                        # Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=olivia.johnson@adelaide.edu.au          # Email to which notifications will be sent

export HOME=/hpcfs/users/a1704225

params=${1} #parameters, pass from text file outside slurm file (could use for loop to cycle through multiple sets outside of this script)
sim_type=${2}
module load Anaconda3/2020.07
conda activate slim_4

# directories
param_dir="/hpcfs/users/a1704225/parameters/single_locus/${sim_type}/"
p_name=$(basename -s ".txt" ${params});

results_dir="/hpcfs/users/a1704225/results/single_locus/${sim_type}/${p_name}"

#python script ${paramter file} ${array id} ${job id for tmpdir path in SLiM}
python3 /hpcfs/users/a1704225/oliviaphd/hpc/single_locus_analysis.py ${params} ${SLURM_ARRAY_TASK_ID} ${results_dir} ${sim_type}

conda deactivate

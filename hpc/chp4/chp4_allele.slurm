#!/bin/bash -l
#SBATCH -p skylake,icelake,a100cpu,v100cpu
#SBATCH --array=2-9
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=60:00:00
#SBATCH --mem=3G
#SBATCH --gres=tmpfs:2G
#SBATCH --mail-type=END                                         # Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL                                        # Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=olivia.johnson@adelaide.edu.au          # Email to which notifications will be sent

export HOME=/hpcfs/users/a1704225

params=${1} #parameters, pass from text file outside slurm file (could use for loop to cycle through multiple sets outside of this script)
module purge
module use /apps/skl/modules/all
module load Anaconda3/2020.07
conda activate slim_4

# directories
param_dir="/hpcfs/users/a1704225/parameters/chp4_alleles/"
p_name=$(basename -s ".txt" ${params});

results_dir="/hpcfs/users/a1704225/results/chp4/${p_name}"
mkdir -p ${results_dir}

#python script ${paramter file} ${array id} ${job id for tmpdir path in SLiM}
python /hpcfs/users/a1704225/oliviaphd/hpc/chp4/chp4_allele.py ${params} ${SLURM_ARRAY_TASK_ID} ${TMPDIR}

#copy files from tmpdir to results folder
cp ${TMPDIR}/al_freq_${params}_${SLURM_ARRAY_TASK_ID}.txt ${results_dir}
cp ${TMPDIR}/fitness_${params}_${SLURM_ARRAY_TASK_ID}.txt ${results_dir}

#!/bin/bash -l
#SBATCH --array=1-20
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=16:00:00
#SBATCH --mem=52G
#SBATCH --gres=tmpfs:6G
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
mkdir -p ${results_dir}

#python script ${paramter file} ${array id} ${job id for tmpdir path in SLiM}
python /hpcfs/users/a1704225/oliviaphd/hpc/single_locus_run.py ${params} ${SLURM_ARRAY_TASK_ID} ${TMPDIR} ${results_dir} ${sim_type}

#copy files from tmpdir to results folder
cp ${TMPDIR}/rec_map_${params}.txt ${results_dir}
cp ${TMPDIR}/burnin_seglift_${params}_${SLURM_ARRAY_TASK_ID}.trees ${results_dir}
cp ${TMPDIR}/treeseq_${params}_${SLURM_ARRAY_TASK_ID}.trees ${results_dir}
cp ${TMPDIR}/al_freq_${params}_${SLURM_ARRAY_TASK_ID}.txt ${results_dir}
cp ${TMPDIR}/pedigree_${params}_${SLURM_ARRAY_TASK_ID}.txt ${results_dir}
cp ${TMPDIR}/slimlog_${params}_${SLURM_ARRAY_TASK_ID}.txt ${results_dir}

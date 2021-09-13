
#!/bin/bash
#SBATCH --array=1-1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time= 1:00:00
#SBATCH --mem=6G
#SBATCH --gres=tmpfs:2G
#SBATCH -o /parameters/slurm_out-%j
#SBATCH --mail-type=END                                         # Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL                                        # Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=olivia.johnson@adelaide.edu.au          # Email to which notifications will be sent


params=${1}
module load Anaconda3/2020.07
conda activate sim_env

# directories
param_dir="/hpcfs/users/a1704225/parameters"
p_name=$(basename -s ".txt" ${params});

results_dir="/hpcfs/users/a1704225/results/${p_name}"
mkdir -p ${results_dir}


python allele_run.py ${param_dir}/${params} ${SLURM_ARRAY_TASK_ID} ${TMPDIR}

cp ${TMPDIR}/reproductiveOutput_${params}_${SLURM_ARRAY_TASK_ID}.txt ${results_dir}
cp ${TMPDIR}/al_freq_${params}_${SLURM_ARRAY_TASK_ID}.txt ${results_dir}
cp ${TMPDIR}/fitness_${params}_${SLURM_ARRAY_TASK_ID}.txt ${results_dir}

#!/bin/bash -l
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=3:00:00
#SBATCH --mem=2G
#SBATCH --mail-type=END                                         # Send a notification email when the job is done (=END)
#SBATCH --mail-type=FAIL                                        # Send a notification email when the job fails (=FAIL)
#SBATCH --mail-user=olivia.johnson@adelaide.edu.au          # Email to which notifications will be sent

export HOME=/hpcfs/users/a1704225

module load Anaconda3/2020.07
conda activate slim_4

#python script to run
slim /hpcfs/users/a1704225/oliviaphd/benchmarking/burnin_mut.slim

conda deactivate

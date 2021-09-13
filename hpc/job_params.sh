param_dir="/hpcfs/users/a1704225/parameters"

for params in $(ls ${param_dir}/*.txt) {
do sbatch seglift_slurm.slurm ${params}
done
    }

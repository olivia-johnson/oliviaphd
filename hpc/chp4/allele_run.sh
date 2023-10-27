
params=${1} #parameters, pass from text file outside slurm file (could use for loop to cycle through multiple sets outside of this script)
conda activate olivia_chp4

# directories
param_dir="/localscratch/olivia/parameters/chp4_alleles/"
p_name=$(basename -s ".txt" ${params});

TMPDIR="/localscratch/olivia/chp4/${p_name}"
mkdir -p ${results_dir}
for i in seq '1 20'; do
#python script ${paramter file} ${array id} ${job id for tmpdir path in SLiM}
python /localscratch/oliviaphd/hpc/chp4/chp4_allele.py ${params} ${i} ${TMPDIR}; done;

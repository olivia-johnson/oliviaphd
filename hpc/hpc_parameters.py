import os
import yaml

os.chdir("/Users/olivia/phd_data/hpc_parameters/")

param_label = "CP_EG_L10_d65_y05_allele"

parameters= {"params":param_label,
"slim_sim": "hpc_seglift_allele",
"nChrom":10,
"chromSize":1,
"recRate":0,
"mutRate":0,
"s_pop":1e4,
"w_pop":1e4,
"l":10,
"y":0.5,
"d":0.65,
"rGen":750,
"fitness_on":1,
"sum_gen":15,
"win_gen":15,
"winpChrom":0,
"burnin_Ne":0}



with open('{0}.txt'.format(param_label), 'w') as outfile:
    yaml.dump(parameters, outfile, default_flow_style=False)

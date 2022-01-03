import os
import yaml

os.chdir("/Users/olivia/phd_data/hpc_parameters/")

##param_label = "CP_UG_L1000_d65_y05_allele"
group =84

parameters= {"group":group,
"slim_sim": "hpc_seglift_l10",
"nChrom":21,
"chromSize":5e5,
"recRate":1e-8,
"mutRate":3e-9,
"s_pop":1e4,
"w_pop":1e4,
"l":1000,
"y":4,
"d":0.65,
"rGen":3000,
"fitness_on":1,
"sum_gen":15,
"win_gen":15,
"winpChrom":0,
"burnin_Ne":1e4}



with open('group_{0}.txt'.format(group), 'w') as outfile:
    yaml.dump(parameters, outfile, default_flow_style=False)

import os
import yaml

os.chdir("/Users/olivia/phd_data/hpc_parameters/")

##param_label = "CP_UG_L1000_d65_y05_allele"
group =29
sim="background"

parameters= {"group":group,
"slim_sim": sim,
"nChrom":21,
"chromSize":5e5,
"recRate":1e-6,
"mutRate":1e-7,
"s_pop":1e4,
"w_pop":1e4,
"l":10,
"y":4,
"d":0.65,
"rGen":7500,
"fitness_on":1,
"sum_gen":13,
"win_gen":2,
"winpChrom":51,
"burnin_Ne":1e4}



with open('{1}/group_{0}.txt'.format(group,sim), 'w') as outfile:
    yaml.dump(parameters, outfile, default_flow_style=False)


group =4
sim="capped"

parameters= {"group":group,
"slim_sim": sim,
"nChrom":21,
"chromSize":5e5,
"recRate":1e-6,
"mutRate":1e-7,
"s_pop":1e4,
"w_pop":1e4,
"l":10,
"y":4,
"d":0.65,
"rGen":7500,
"fitness_on":1,
"sum_gen":15,
"win_gen":15,
"winpChrom":51,
"burnin_Ne":1e4,
"offCap":350}



with open('{1}/group_{0}.txt'.format(group,sim), 'w') as outfile:
    yaml.dump(parameters, outfile, default_flow_style=False)
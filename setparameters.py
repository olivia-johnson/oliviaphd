import yaml
import os

group = 1
sim_type = "seglift_linked_allele"
##"seglift_long"
os.chdir("/Users/olivia/phd_data/{0}".format(sim_type))
path  = "./group_" + str(group)
os.mkdir(path)

os.chdir("/Users/olivia/phd_data/{0}/group_{1}".format(sim_type, group))

parameters= {"sim_type": sim_type,
"slim_sim": "seglift_linked_allele",
    ##"seglift_long",
"group":group, 
"runs" : 8,
"nChrom":50,
"chromSize":2,
"recRate":0,
"mutRate":3e-9,
"s_pop":1e4,
"w_pop":1e4,
"l":100,
"y":0.5,
"d":0.65,
"rGen":3000,
"fitness_on":1,
"sum_gen":10,
"win_gen":10,
"winpChrom":51,
"burnin_Ne":0}

with open('parameters.yml', 'w') as outfile:
    yaml.dump(parameters, outfile, default_flow_style=False)

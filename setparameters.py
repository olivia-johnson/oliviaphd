import yaml
import os

group = 13
sim_type = "seglift_long"
os.chdir("/Users/olivia/oliviaphd/{0}".format(sim_type))
path  = "./group_" + str(group)
os.mkdir(path)

os.chdir("/Users/olivia/oliviaphd/{0}/group_{1}".format(sim_type, group))

parameters= {"sim_type": sim_type,
"slim_sim": "seglift_long_unselchrom",
"group":group, 
"runs" : 8,
"nChrom":21,
"chromSize":5e5,
"recRate":1e-8,
"mutRate":3e-9,
"s_pop":1e4,
"w_pop":1e3,
"l":10,
"y":0.5,
"d":0.65,
"rGen":10000,
"fitness_on":1,
"sum_gen":13,
"win_gen":2,
"winpChrom":51}

with open('parameters.yml', 'w') as outfile:
    yaml.dump(parameters, outfile, default_flow_style=False)
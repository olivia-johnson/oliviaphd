import os
os.chdir("/Users/olivia/oliviaphd/")
import multiprocess as mp
import msprime
import pyslim
import numpy as np
import tskit
import pandas as pd
import time
import itertools
import allel
import yaml
import recombination
import ts_analysis

sim_code= "seglift_long"  ## DONT CHANGE
os.chdir("/Users/olivia/oliviaphd/{0}".format(sim_code))
import seglift_long
##"seglift_win"

sim_type = "seglift_long"
#"seglift_substitute"
#
group = 9

os.chdir("/Users/olivia/oliviaphd/{0}/group_{1}".format(sim_type, group))


####  READ IN PARAMETERS
    # load in parameter fiile
with open('./parameters.yml', 'r') as f:
    parameters = yaml.load(f, Loader=yaml.FullLoader)
    
    
    #set parameters from file
group = int(parameters["group"])
slim_sim = parameters["slim_sim"]
runs = int(parameters["runs"])
nChrom = int(parameters["nChrom"])
chromSize = int(parameters["chromSize"])
recRate = parameters["recRate"]
mutRate=parameters["mutRate"]
s_pop = int(parameters["s_pop"])
w_pop = int(parameters["w_pop"])
l = int(parameters["l"])
y = parameters["y"]
d = parameters["d"]
rGen=int(parameters["rGen"])
fitness_on = parameters["fitness_on"]
sum_gen = int(parameters["sum_gen"])
win_gen = int(parameters["win_gen"])
winpChrom = parameters["winpChrom"]

####  GENERATE RECOMBINATION MAP

rec_map = recombination.recombination_map(group, nChrom, chromSize, recRate)

####  SIMULATE BURNIN
processes=[]
for sim_run in range(runs):
    p=mp.Process(target=seglift_long.simulate_burnin, args=[group, sim_run, rec_map, s_pop])
    p.start()
    processes.append(p)
    
for process in processes:
    process.join()

####  SIMULATE SEGLIFT
processes=[]

for sim_run in range(runs):
    p=mp.Process(target=seglift_long.simulate_seglift, args=[sim_type, slim_sim, group, sim_run, recRate, nChrom, chromSize, s_pop, w_pop, l, y, d, rGen, fitness_on, sum_gen, win_gen])
    p.start()
    processes.append(p)
    
for process in processes:
    process.join()



#### ANALYSE TREE SEQUENCE WITH TS_ANALYSIS
nWin=winpChrom*nChrom

processes=[]
start_time = time.time()
for run in range(runs):
    sim_run = run
    print(sim_run)
    p=mp.Process(target=ts_analysis.analyse, args=[group, sim_type , sim_run, mutRate, l, nChrom, nWin, sum_gen, win_gen])
    p.start()
    processes.append(p)
    
for process in processes:
    process.join()

print("Time = ", (time.time() - start_time))
import os
import sys
import msprime
import pyslim
import numpy as np
import yaml
import tskit
import pandas as pd
import time
import allel
sys.path.insert(1, '/home/a1704225/oliviaphd/')
import seglift_hpc

params=sys.argv[1]
sim_run = sys.argv[2]
tmpdir =str(sys.argv[3])
results_dir=str(sys.argv[4])
####  READ IN PARAMETERS
    # load in parameter file
with open('/hpcfs/users/a1704225/parameters/background/{0}.txt'.format(params), 'r') as f:
    parameters = yaml.load(f, Loader=yaml.FullLoader)

    #set parameters from file
slim_sim = parameters["slim_sim"]
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
group=parameters["group"]
burnin_Ne = int(parameters["burnin_Ne"])


start_time = time.time()
## GENERATE RECOMBINATION MAP
rec_map = seglift_hpc.recombination_map(tmpdir, group, l, nChrom, chromSize, recRate)

####  SIMULATE BURNIN
if os.path.exists('{0}_{1}.trees'.format(results_dir, sim_run))==False:
    seglift_hpc.simulate_burnin(tmpdir, results_dir, group, l, sim_run, rec_map, s_pop, burnin_Ne, chromSize, nChrom)

####  SIMULATE SEGLIFT
if slim_sim =="background":
    seglift_hpc.simulate_seglift(tmpdir, group, sim_run, recRate, nChrom, chromSize, s_pop, w_pop, l, y, d, rGen, fitness_on, sum_gen, win_gen)
if slim_sim == "background_mfit":
    seglift_hpc.simulate_seglift_mfit(tmpdir, results_dir, group, sim_run, recRate, nChrom, chromSize, s_pop, w_pop, l, y, d, rGen, fitness_on, sum_gen, win_gen)

#### ANALYSE TREE SEQUENCE WITH TS_ANALYSIS
#nWin=winpChrom*nChrom

#seglift_hpc.analyse(tmpdir, group, sim_run, mutRate, l, nChrom, nWin, sum_gen, win_gen)

print("Time = ", (time.time() - start_time))

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
sys.path.insert(1, '/hpcfs/users/a1704225/scripts/oliviaphd/')
import single_locus_hpc

params=sys.argv[1]
sim_run = sys.argv[2]
tmpdir =str(sys.argv[3])
results_dir=str(sys.argv[4])
sim_type=str(sys.argv[5])
####  READ IN PARAMETERS
    # load in parameter file
with open('/hpcfs/users/a1704225/parameters/single_locus/{1}/{0}.txt'.format(params, sim_type), 'r') as f:
    parameters = yaml.load(f, Loader=yaml.FullLoader)

    #set parameters from file
genomeSize = int(parameters["genomeSize"])
recRate = parameters["recRate"]
s_pop = int(parameters["s_pop"])
w_pop = int(parameters["w_pop"])
h_s = parameters["h_s"]
h_w = parameters["h_w"]
s_s = parameters["s_s"]
s_w = parameters["s_w"]
rGen=int(parameters["rGen"])
fitness_on = parameters["fitness_on"]
sum_gen = int(parameters["sum_gen"])
win_gen = int(parameters["win_gen"])
group=parameters["group"]
winSize = parameters["winSize"]
freq=int(parameters["f"])


start_time = time.time()

## Calculate Burnin Ne
burnin_Ne = round((sum_gen+win_gen)/(((1/s_pop)*sum_gen)+((1/w_pop)*win_gen)))

####  SIMULATE BURNIN
rate_map, sequenceSize = single_locus_hpc.recombination_map(tmpdir, group, genomeSize, recRate, winSize)
if os.path.exists('{0}/burnin_seglift_group_{1}_{2}.trees'.format(results_dir, group,sim_run))==False:
    single_locus_hpc.single_locus_burnin(tmpdir, group, sim_run, sequenceSize, s_pop, burnin_Ne, recRate)

####  SIMULATE SEGLIFT

single_locus_hpc.simulate_single_locus(tmpdir, results_dir, group, sim_run, recRate, genomeSize, s_pop, w_pop, h_s, h_w, s_s, s_w, rGen, fitness_on, sum_gen, win_gen, freq, sim_type, winSize)



print("Time = ", (time.time() - start_time))


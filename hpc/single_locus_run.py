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
import single_locus_hpc

params=sys.argv[1]
sim_run = sys.argv[2]
tmpdir =str(sys.argv[3])
results_dir=str(sys.argv[4])
####  READ IN PARAMETERS
    # load in parameter file
with open('/hpcfs/users/a1704225/parameters/single_locus/{0}.txt'.format(params), 'r') as f:
    parameters = yaml.load(f, Loader=yaml.FullLoader)

    #set parameters from file
genomeSize = int(parameters["chromSize"])
recRate = parameters["recRate"]
mutRate=parameters["mutRate"]
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
winpChrom = parameters["winpChrom"]
group=parameters["group"]
burnin_Ne = int(parameters["burnin_Ne"])


start_time = time.time()

####  SIMULATE BURNIN
if os.path.exists('{0}/burnin_seglift_group_{1}_{2}.trees'.format(results_dir, group,sim_run))==False:
    single_locus_hpc.single_locus_burnin(tmpdir, group, sim_run, genomeSize, s_pop, burnin_Ne, recRate)

####  SIMULATE SEGLIFT

single_locus_hpc.simulate_single_locus(tmpdir, results_dir, group, sim_run, recRate, genomeSize, s_pop, w_pop, h_s, h_w, s_s, s_w, rGen, fitness_on, sum_gen, win_gen)

#### ANALYSE TREE SEQUENCE WITH TS_ANALYSIS
#nWin=winpChrom

#seglift_hpc.analyse(tmpdir, group, sim_run, mutRate, l, nChrom, nWin, sum_gen, win_gen)

print("Time = ", (time.time() - start_time))


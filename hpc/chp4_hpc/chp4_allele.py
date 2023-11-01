import os
import sys
import msprime
import pyslim
import numpy as np
import yaml
import tskit
import pandas as pd
import time
sys.path.insert(1, '/hpcfs/users/a1704225/oliviaphd/hpc/')
import chp4_functions

params=sys.argv[1]
sim_run = sys.argv[2]
tmpdir =str(sys.argv[3])
####  READ IN PARAMETERS
    # load in parameter file
with open('/hpcfs/users/a1704225/parameters/{0}.txt'.format(params), 'r') as f:
    parameters = yaml.load(f, Loader=yaml.FullLoader)

    #set parameters from file
s_pop = int(parameters["s_pop"])
w_pop = int(parameters["w_pop"])
l = int(parameters["l"])
y = parameters["y"]
rGen=int(parameters["rGen"])
fitness_on = parameters["fitness_on"]
sum_gen = int(parameters["sum_gen"])
win_gen = int(parameters["win_gen"])
group=parameters["group"]

start_time = time.time()
####  SIMULATE SEGLIFT
chp4_functions.simulate_alleles(tmpdir, group, sim_run, s_pop, w_pop, l, y, rGen, fitness_on, sum_gen, win_gen)

print("Time = ", (time.time() - start_time))

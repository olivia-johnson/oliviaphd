import os
import sys
import msprime
import pyslim
import numpy as np
import tskit
import pandas as pd
import time
import allel
sys.path.insert(1, '/home/a1704225/oliviaphd/hpc/')
import seglift_hpc

params=sys.getenv('params')
sim_run = sys.getenv('SLURM_ARRAY_TASK_ID')
tmpdir = sys.getenv('TMPDIR')
####  READ IN PARAMETERS
    # load in parameter file
parameters = open('/hpcfs/users/a1704225/parameters/{0}.txt'.format(params), 'r')
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
burnin_Ne = parameters["burnin_Ne"]


####  SIMULATE SEGLIFT
seglift_hpc.simulate_seglift(tmpdir, slim_sim, params, sim_run, recRate, nChrom, chromSize, s_pop, w_pop, l, y, d, rGen, fitness_on, sum_gen, win_gen)



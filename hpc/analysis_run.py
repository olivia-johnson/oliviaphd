#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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
sim_type =str(sys.argv[4])


with open('/hpcfs/users/a1704225/parameters/{1}/{0}.txt'.format(params, sim_type), 'r') as f:
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
#### ANALYSE TREE SEQUENCE WITH TS_ANALYSIS

print("Group: " + str(group))
nWin=winpChrom*nChrom

seglift_hpc.analyse(tmpdir, group, sim_run, mutRate, l, nChrom, nWin, sum_gen, win_gen)

print("Time = ", (time.time() - start_time))
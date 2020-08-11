#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 14:17:48 2020

@author: olivia
"""

import os
os.chdir("/Users/olivia/oliviaphd/")

import msprime
import pyslim
import tskit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import itertools


genomeSize = int(1e6)
popnSize = int(1e4)
mutRate = 1e-5
recRate = 1e-8
y = 2.0
d = 0.6
nWin = 100
rGen = 100
sum_gen = 8#no. summer generations
win_gen = 3#no. winter generations

## COALESCENT BURN IN

#burnin = msprime.simulate(sample_size=2*popnSize,Ne=popnSize, length=genomeSize, mutation_rate=mutRate, recombination_rate=recRate)
# issue is having mutations in this tree. Turn mut rate to 0. ALso, makes sense, as we only need the geneaology from the burn in, and we add all muations at the end of combined coalescent + forward time.
burnin = msprime.simulate(sample_size=popnSize,Ne=int(1e6), length=genomeSize, mutation_rate=mutRate, recombination_rate=recRate)
burnin_ts = pyslim.annotate_defaults(burnin, model_type="WF", slim_generation=1)
burnin_ts.dump("burnin_SL.trees")

## FORWARD SIMULATION
# for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
cmd = "slim -d GenomeSize=" + str(int(genomeSize)) + " -d N=" + str(int(popnSize)) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=" + str(mutRate) +" -d rr=" + str(recRate) +" ~/oliviaphd/SL.slim"
print(cmd)
os.system(cmd)

slim_ts = pyslim.load("./treeseq_SL.trees").simplify()
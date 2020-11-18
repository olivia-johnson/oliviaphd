import os
os.chdir("/Users/olivia/oliviaphd/seligt_no_treeseq/")


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import itertools
import random


group=1
runs = 1
genomeSize = int(1e6)
popnSize = int(5e3)
mutRate = 1e-6
recRate = 1e-8
l = 20
y = 2.0
d = 0.6
nWin = 20
sum_gen = 8#no. summer generations
win_gen = 3#no. winter generations


## BURN IN ##

for i in range(5):
    sim_runt = time.time()
    cmd = "slim -d n_burnin=" + str(i)+" -d GenomeSize=" + str(int(genomeSize)) + " -d N=" + str(int(popnSize)) + " -d mut=" +str(mutRate)+" -d rr=" + str(recRate) +" ~/oliviaphd/seglift_no_treeseq/sim_burnin.slim"
    print(cmd)
    os.system(cmd)
    print("Time for burnin = ", (time.time()- sim_runt), "\n")


for x in range(runs):
    #print(x)
    sim_run = str(x+1)
    sim_runt = time.time()

    n_burnin = random.randint(0,4)
    
    ## FORWARD SIMULATION
    
    cmd = "slim -d n_burnin="+str(n_burnin)+" -d group="+str(group)+" -d sim_run=" + str(sim_run)+" -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d N=" + str(int(popnSize)) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=" + str(mutRate) + " -d rr=" + str(recRate) +" ~/oliviaphd/seglift_no_treeseq/seglift_nts.slim"
    print(cmd)
    os.system(cmd)



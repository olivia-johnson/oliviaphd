import os
os.chdir("/Users/olivia/oliviaphd/seglift_even_dist/wittmann_rep")

import msprime
import pyslim
import tskit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import itertools
import random
import allel

group =0 ## identifier for a set of parameters
runs = 10
genomeSize = int(100)
popnSize = int(1e3)
mutRate = 1e-6
recRate = 0.5
l = 100
y = 4
d = 0.6
nWin = 20
rGen = 15
fitness_on = 1
sum_gen = 15#no. summer generations
win_gen = 15#no. winter generations


for d in (0.15, 0.5, 0.65):
    print (d)
    # for y in (0.5, 1, 2 , 4):
    #     print (y)
        
        
    group = group +1
    path  = "./ed_" + str(group)
    os.mkdir(path)
    for x in range(runs): 
        #print(x)
        sim_run = str(x+1)
        sim_runt = time.time()
        #print("./data/seglift_ts/burnin_seglift_ts{0}.trees".format(sim_run))
        
        
    ## COALESCENT BURN IN    
        start_time = time.time()
        burnin = msprime.simulate(sample_size=2*popnSize,Ne=popnSize, length=genomeSize, mutation_rate=0, recombination_rate=recRate)
        burnin_ts = pyslim.annotate_defaults(burnin, model_type="WF", slim_generation=1)
        burnin_ts.dump("./ed_{0}/burnin_seglift_ed_{1}_{2}.trees".format(group,group,sim_run))
        print("Time for burnin = ", (time.time()- start_time))
        
        ## FORWARD SIMULATION
        # for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
        start_time = time.time()
        cmd = "slim -d fit="+ str(fitness_on)+" -d group=" + str(group) + " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)+" -d sim_run=" + str(sim_run) + " -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d N=" + str(int(popnSize)) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=0.0 -d rr=" + str(recRate) + " -d rGen="+ str(rGen) +" ~/oliviaphd/seglift_even_dist/wittmann_rep/seglift_witt_rep.slim"
        print(cmd)
        os.system(cmd)
        print("Time for SLiM sim = ", (time.time()- start_time))
        
        #slim_ts = pyslim.load("./ed_{0}/treeseq_seglift_ed_{1}_{2}.trees".format(group,group,sim_run)).simplify()

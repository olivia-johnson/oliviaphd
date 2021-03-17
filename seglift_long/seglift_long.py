import os
os.chdir("/Users/olivia/oliviaphd/seglift_long")

import multiprocess as mp
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

group =9 ## identifier for a set of parameters
runs = 8
chrom = 21
chrom_size =int(5e5)
genomeSize = chrom_size*chrom
popnSize = int(1e4)
mutRate = 1e-6
recRate = 1e-8
l = 10
y = 4
d = 0.65
rGen = 10000
fitness_on = 1
sum_gen = 15#no. summer generations
win_gen = 15#no. winter generations

        
group = group +1
path  = "./group_" + str(group)
os.mkdir(path)

## GENERATE CHROMOSMES ##
## create recombination map for msprime and slim to simulate unlinked chromsomes
rec_rows=[]
for c in range(chrom):
    print(c)
    rec_dict = {}
    rec_dict.update({"positions": c*chrom_size})
    rec_dict.update({"rates":recRate})
    rec_rows.append(rec_dict)
    rec_dict = {}
    rec_dict.update({"positions": (c+1)*chrom_size-1})
    rec_dict.update({"rates":0.5})
    rec_rows.append(rec_dict)
rec_data = pd.DataFrame(rec_rows)


## generate slim recombination map
slim_rec = rec_data.positions[1:]
slim_rec = slim_rec.reset_index(drop=True)
slim_rec= pd.concat([slim_rec,rec_data.rates[0:-1]],axis=1)


## formulate msprime recombination map
rec_data.positions.iloc[-1]=rec_data.positions.iloc[-1]+1
rec_map = msprime.RecombinationMap(positions = list(rec_data.positions), rates= list(rec_data.rates), num_loci = chrom*chrom_size-1)



slim_rec.to_csv("./group_{0}/rec_map.txt".format(group), index=False, header = False, sep = "\t")

# for x in range(runs): 
#     #print(x)
#     sim_run = str(x+1)
#     sim_runt = time.time()
#     #print("./data/seglift_ts/burnin_seglift_ts{0}.trees".format(sim_run))
def simulate_seglift(sim_run): 
    ## COALESCENT BURN IN    
    start_time = time.time()
    burnin = msprime.simulate(sample_size=2*popnSize,Ne=popnSize, mutation_rate=0, recombination_map=rec_map)
    burnin_ts = pyslim.annotate_defaults(burnin, model_type="WF", slim_generation=1)
    burnin_ts.dump("./group_{0}/burnin_seglift_long_{1}_{2}.trees".format(group,group,sim_run))
    print("Time for burnin = ", (time.time()- start_time))
    
    ## FORWARD SIMULATION
    # for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
    #start_time = time.time()
    cmd = "slim -d fit="+ str(fitness_on)+" -d group=" + str(group) + " -d nChrom=" + str(chrom)+" -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)+" -d sim_run=" + str(sim_run) + " -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d N=" + str(int(popnSize)) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=0.0 -d rr=" + str(recRate) + " -d rGen="+ str(rGen) +" ~/oliviaphd/seglift_long/seglift_long.slim"
    print(cmd)
    os.system(cmd)
    #print("Time for SLiM sim = ", (time.time()- start_time))
    print("Simulations took ", (time.time()-start_time) ,  " seconds")
    

processes=[]

for sim_run in range(runs):
    p=mp.Process(target=simulate_seglift, args=[sim_run])
    p.start()
    processes.append(p)
    
for process in processes:
    process.join()
  
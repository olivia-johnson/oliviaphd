import os
os.chdir("/Users/olivia/oliviaphd/seglift_treeseq/")

import msprime
import pyslim
import tskit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import itertools
import random


genomeSize = int(1e6)
popnSize = int(1e4)
runs = 20
mutRate = 1e-6
recRate = 1e-8
l = 20
y = 2.0
d = 0.6
nWin = 20
rGen = 100
sum_gen = 8#no. summer generations
win_gen = 3#no. winter generations


## FUNCTIONS  ##
def allele_counts(ts, sample_sets=None):
    if sample_sets is None: 
        sample_sets = [ts.samples()]
    def f(x):
        return x
    return ts.sample_count_stat(sample_sets, f, len(sample_sets),
                       span_normalise=False, windows='sites',
                       polarised=True, mode='site', strict=False)


def sum_stats(ts, wins, inds, mut_tab, sample_sets=None):
        if sample_sets is None:
           sample_sets = [ts.samples()]
        rows_list3 = []
       
        #create windows
        win = np.linspace(0, ts.sequence_length, num=wins+1)
        
        ## seperated into ts for each timepoint
        for t in np.unique(inds.time):
            sample = inds.nodes[inds.time == t]
            samples= list(itertools.chain(*sample))
            # samp_ts = ts.subset(samples) ## may not be keeping history of nodes
            samp_ts = ts.simplify(samples = samples) ## should keep history of nodes
            
        # calculate Tajima's D for windows
            tajds =  samp_ts.Tajimas_D(sample_sets=None, windows=win, mode="site")
            tajdb =  samp_ts.Tajimas_D(sample_sets=None, windows=win, mode="branch") ## doesn change between windows by site more informative
    
    
            div = samp_ts.diversity(sample_sets = None, windows = win)
            for w in range(wins):
                # print(w)
                s_mutcount = np.sum((mut_tab.mut_pos[mut_tab.mut_type=="fluctuating"] >= win[w]) & (mut_tab.mut_pos[mut_tab.mut_type=="fluctuating"] < (win[w+1]-1)))
                n_mutcount = np.sum((mut_tab.mut_pos[mut_tab.mut_type=="neutral"] >= win[w]) & (mut_tab.mut_pos[mut_tab.mut_type=="neutral"] < (win[w+1]-1)))
               
                dict3={}
                dict3.update({"time":t})
                dict3.update({"n_win":w})
                dict3.update({"win_start" : win[w]})
                dict3.update({"win_end" : win[w+1]-1})
                dict3.update({"n_s_mut" : s_mutcount})
                dict3.update({"n_n_mut" : n_mutcount})
                dict3.update({"tajimas_d_site":tajds[w]})
                dict3.update({"tajimas_d_branch":tajdb[w]})
                dict3.update({"diversity": div[w]})
                rows_list3.append(dict3)
    ### REMOVE BRACKETS
        s_stats = pd.DataFrame(rows_list3) 

        return(s_stats)
    

##  RUN ##
for x in range(runs):
   # print(x)
    sim_run = str(x)
    sim_runt = time.time()
    #print("./data/seglift_ts/burnin_seglift_ts{0}.trees".format(sim_run))
    
    
## COALESCENT BURN IN    
    #burnin = msprime.simulate(sample_size=2*popnSize,Ne=popnSize, length=genomeSize, mutation_rate=mutRate, recombination_rate=recRate)
    # issue is having mutations in this tree. Turn mut rate to 0. ALso, makes sense, as we only need the geneaology from the burn in, and we add all muations at the end of combined coalescent + forward time.
    burnin = msprime.simulate(sample_size=popnSize,Ne=int(1e6), length=genomeSize, mutation_rate=0, recombination_rate=recRate)
    burnin_ts = pyslim.annotate_defaults(burnin, model_type="WF", slim_generation=1)
    burnin_ts.dump("./burnin/burnin_seglift_ts_{0}.trees".format(sim_run))
    
    ## FORWARD SIMULATION
    # for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
    cmd = "slim -d sim_run=" + str(sim_run) + " -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d N=" + str(int(popnSize)) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=0.0 -d rr=" + str(recRate) +" ~/oliviaphd/scripts/seglift_ts.slim"
    print(cmd)
    os.system(cmd)
    
    slim_ts = pyslim.load("./slim_out/treeseq_seglift_ts_{0}.trees".format(sim_run)).simplify()
    
    ## mutation check
    if (slim_ts.num_sites != l):
        print("Less than " + str(l) + "introduced mutations")
    else:
        print (str(l) + " introduced mutations")
        
    ## SUMMARISE MUTATIONS
    
    start_time = time.time()
    mut_met = pd.DataFrame({"mut_num": [], "mut_pos": []})
    for mut in slim_ts.mutations():
        #print(mut)
        mut_met = mut_met.append({"mut_num" : mut.site, "mut_pos" : mut.position}, ignore_index=True)
    
    mut_met = mut_met.loc[mut_met.astype(str).drop_duplicates(subset=("mut_num")).index]
    
    mut_freq= np.loadtxt("./slim_out/sim_data_{0}.txt".format(sim_run), delimiter = ",", skiprows=(7 + l)) ## import freqs from SLiM
    #for pos in mut_met.mut_pos:  ## label non-segregating seasonal loci
        #print(pos)
       #mts = mut_freq[:,2 == pos]
        #if mut_freq[i,2] == 0 | mut_freq[i,2] == 1:
            
    print("Time for initial mut info = ", (time.time()- start_time))
    
    
    ## SUMMARISE INDIVIDUALS
    start_time = time.time()
    rows_list = []
    
    for ind in slim_ts.individuals():
     
      #  print(ind)
        if ind.time % 10 == 0:
            ind_season = "S"
        else:
            ind_season = "W"
    
        dict1 = {}
            # get input row in dictionary format
            # key = col_name
        dict1.update({"id" : ind.id})
        dict1.update({"time" : ind.time})
        dict1.update({"pop" : ind.population})
        dict1.update({"season" : ind_season})
        dict1.update({"nodes": ind.nodes}) 
    
        rows_list.append(dict1)
    
    ind_met = pd.DataFrame(rows_list) 
    print("Time for ind table = ", (time.time()- start_time))
    
    
        
    ## ADD NEUTRAL MUTATIONS
    mut_ts = pyslim.SlimTreeSequence(msprime.mutate(slim_ts, rate=mutRate, keep=True))
    
    #mut_gm = mut_ts.genotype_matrix()  ##  IF WANT TO USE GENOTYPE MATRIX
    
    start_time = time.time()
    
    
    rows_list2 = []
    a_list = list(mut_met.mut_pos)
    for mut in mut_ts.mutations():
        #print(mut)
    
            #print(mut)
            
        # given_value = mut.position
            
        # absolute_difference_function = lambda list_value : abs(list_value - given_value)
    
        # closest_value = min(a_list, key=absolute_difference_function)
            
        # dist = given_value - closest_value
        
        #num = mut_met.mut_num[mut_met.mut_pos == closest_value]
        if mut.metadata == []:
            m_type = "neutral"
        else:
            m_type = "fluctuating"
            
        dict2 = {}
            # get input row in dictionary format
            # key = col_name
        dict2.update({"mut_id":mut.id})
        dict2.update({"mut_num" : mut.site})
        dict2.update({"mut_pos" : mut.position})
        dict2.update({"mut_type":m_type})
        # dict2.update({"nearest_dist": abs(dist)}) 
        # dict2.update({"nearest_mut_pos": closest_value})
       # dict2.update({"nearest_mut_num": int(num)})
    
        rows_list2.append(dict2)
    
    n_met = pd.DataFrame(rows_list2) 
    
    
    print("Time for mut table = ", (time.time()- start_time))
    
    
    o_len = len(n_met)
    s_len = len(n_met[n_met.mut_type == 'fluctuating'])
    
    mut_ud = n_met.loc[n_met.astype(str).drop_duplicates(subset=("mut_num")).index]
    n_len = len(mut_ud)
    
    if n_len == (o_len - s_len + l):
        print("Only duplicate selected mutations removed")
    else:
        print("More than duplicate selected mutations removed")
    
    ## REMOVE SELECTED SITES FROM TS
    
    list_m = mut_ud.mut_pos[mut_ud.mut_type=="fluctuating"].astype(int)
    list_mu = list_m .tolist()
    stat_ts = pyslim.SlimTreeSequence(mut_ts.delete_sites(list_mu, record_provenance=True))
      
    ## CALCULATE SUMMARY STATISTICS
    start_time = time.time()
    s_stat = sum_stats(stat_ts, nWin, ind_met, mut_ud, sample_sets=None)
    print("Time for sum stats = ", (time.time()- start_time))
    
    s_stat.to_string(buf = "~/oliviaphd/seglift_treeseq/py_out/sim_s_stat_{0}.txt".format(sim_run))
    
    print("Time for sim run = ", (time.time()- sim_runt))

#### PLOTS


# plt.scatter(s_stat.n_win[s_stat.Gen==9800.0], s_stat.tajimas_d_branch[s_stat.Gen==9800.0], s=2)


# plt.scatter(n_met.nearest_dist[n_met.mut_type == "neutral"], n_met.allele_freq[n_met.mut_type == "neutral"], s=2)
# plt.show()

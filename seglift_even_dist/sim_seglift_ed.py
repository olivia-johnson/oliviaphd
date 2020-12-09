import os
os.chdir("/Users/olivia/oliviaphd/seglift_even_dist/")

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

    ##  PARAMETERS  ##
group =0                ## identifier for a set of parameters
runs = 10               ## number of replicates
genomeSize = int(1e6)   ## size of simulated genome segment
popnSize = int(1e4)     ## number of individuals in population
mutRate = 1e-6          ## mutation rate (used to add mutations onto tree sequence)
recRate = 1e-8          ## recombination rate (used in sims)
l = 10                  ## no. selected loci
y = 2                   ## value for pataemeter that confers epistatic effects in fitness function
d = 0.6                 ## dominace of selected loci
nWin = 100               ## number of windows statistics are calculated over
rGen = 500              ## number of generations between which individuals are remebered in the tree sequence
fitness_on = 1          ## turn fitness fucntion on/off (=0 fitness off)
sum_gen = 10            ## no. summer generations
win_gen = 10            ## no. winter generations


    ##  RUN ##

    
for x in range(runs):       
        ## run replicates
    #print(x)   
        ## set run replicate identifier
    sim_run = str(x+1)   
    #sim_runt = time.time()
    
## COALESCENT BURN IN    
        ## simulate in msprime    
    burnin = msprime.simulate(sample_size=2*popnSize,Ne=popnSize, length=genomeSize, mutation_rate=0, recombination_rate=recRate)
    burnin_ts = pyslim.annotate_defaults(burnin, model_type="WF", slim_generation=1)
        ## write as .trees file
    burnin_ts.dump("./burnin/burnin_seglift_ed_{0}_{1}.trees".format(group,sim_run))
    
    ## FORWARD SIMULATION
    # for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
        ## command for SLiM simulation, parameters set and defined above
    cmd = "slim -d group=" + str(group) + " -d sim_run=" + str(sim_run) + " -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d N=" + str(int(popnSize)) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=0.0 -d rr=" + str(recRate) +" ~/oliviaphd/seglift_even_dist/seglift_ed.slim"
    print(cmd)
        ## run command
    os.system(cmd)
        ## load tree sequence from SLiM simulation
    slim_ts = pyslim.load("./slim_out/treeseq_seglift_ed_{0}_{1}.trees".format(group,sim_run)).simplify()
    
    ## Check number of mutations that were introduced
    if (slim_ts.num_sites != l):
        print("Less than " + str(l) + "introduced mutations")
    else:
        print (str(l) + " introduced mutations")
        
    ## SUMMARISE MUTATIONS
    
    #start_time = time.time()
        ## obtain data for selected sites (used later to remove from tree sequence)
    mut_met = pd.DataFrame({"mut_num": [], "mut_pos": []})
    for mut in slim_ts.mutations():## run through muts in ts
        #print(mut)
        mut_met = mut_met.append({"mut_num" : mut.site, "mut_pos" : mut.position}, ignore_index=True)
        ## remove replicated mutations (as mutations added in SLiM by randomly sampling individuals, goes against tree structure of burnin hence sometimes have duplicates)
    mut_met = mut_met.loc[mut_met.astype(str).drop_duplicates(subset=("mut_num")).index]
   # print("Time for initial mut info = ", (time.time()- start_time))
    
    
    ## SUMMARISE INDIVIDUALS
    #start_time = time.time()
    rows_list = []
        ## obtain data for individuals saved in tree sequence 
    for ind in slim_ts.individuals(): ## run through inds in ts
        #print(ind)
            ## determine what season inds were saved from
        if ind.time % sum_gen == 0:
            ind_season = "W"
        else:
            ind_season = "S"
    
        dict1 = {}
            # get input row in dictionary format
            # key = col_name
        dict1.update({"id" : ind.id})           ## individuals id in ts
        dict1.update({"time" : 10000-ind.time}) ## generation they are from in SLiM time (SLiM and ts count time differently)
        dict1.update({"pop" : ind.population})  ## population individual is from
        dict1.update({"season" : ind_season})   ## season individuals was in
        dict1.update({"nodes": ind.nodes})      ## genotypes individuals contained
    
        rows_list.append(dict1)
    
    ind_met = pd.DataFrame(rows_list) ## convert from dictionary to data frame (fast)
    #print("Time for ind table = ", (time.time()- start_time))
    
    
        
    ## ADD NEUTRAL MUTATIONS
    #start_time = time.time()
    mut_ts = pyslim.SlimTreeSequence(msprime.mutate(slim_ts, rate=mutRate, keep=True))
    #print("Time to ad neutral muts = ", (time.time()- start_time))
    
   # start_time = time.time()
    
    ## SUMMARISE NEUTRAL MUTATIONS
    
    rows_list2 = []
    a_list = list(mut_met.mut_pos)
        ## obtain data for neutral mutations added to ts
    for mut in mut_ts.mutations(): ## run through muts in ts
        #print(mut)
    
            #print(mut)
            ## calculate distance from nearest seleted mutation (commented to increase speed)
        # given_value = mut.position
        # absolute_difference_function = lambda list_value : abs(list_value - given_value)
        # closest_value = min(a_list, key=absolute_difference_function)
        # dist = given_value - closest_value
        #num = mut_met.mut_num[mut_met.mut_pos == closest_value]
        
            ## differentiate mutation type
        if mut.metadata == []:
            m_type = "neutral"
        else:
            m_type = "fluctuating"
            
        dict2 = {}
            # get input row in dictionary format
            # key = col_name
        dict2.update({"mut_id":mut.id})                     ## id of mutation in ts
        dict2.update({"mut_num" : mut.site})                ## number of mutation
        dict2.update({"mut_pos" : mut.position})            ## position of mutation
        dict2.update({"mut_type":m_type})                   ## fluctating or neutral
        # dict2.update({"nearest_dist": abs(dist)})         ## distance between neutral mut and nearest selected loci
        # dict2.update({"nearest_mut_pos": closest_value})  ## position of nearest selected loci
       # dict2.update({"nearest_mut_num": int(num)})        ## number of nearest seelcted loci
    
        rows_list2.append(dict2)
    n_met = pd.DataFrame(rows_list2) ## convert dict to df
    
    #print("Time for mut table = ", (time.time()- start_time))
    

    
    ## REMOVE SELECTED SITES FROM TS
        ##format fluctating loci into list
    list_m = n_met.mut_pos[n_met.mut_type=="fluctuating"].astype(int)
    list_mu = list_m .tolist()
        ## remove selected loci from ts for analysis of neutral loci
    stat_ts = pyslim.SlimTreeSequence(mut_ts.delete_sites(list_mu, record_provenance=True))
      
    ## CALCULATE SUMMARY STATISTICS
    
    rows_list3 = []
    
    #start_time = time.time()
        ## create windows for stats to be calculated in
    win = np.linspace(0, stat_ts.sequence_length, num=nWin+1)
    #print("create windows = ", (time.time()- start_time))
    start_time = time.time()
        ## seperate into ts for each timepoint
    for t in np.unique(ind_met.time):   ## run through time points
        #s_t = time.time()
            ## nodes of indviduals at time t
        sample = ind_met.nodes[ind_met.time == t]
            ## convert to list
        samples= list(itertools.chain(*sample))
            ## create ts of just samples at time t using list
            ## used by tskit to calculate states based on tree structure
        #samp_ts = stat_ts.subset(samples)              ## may not be keeping history of nodes
        samp_ts = stat_ts.simplify(samples = samples)   ## should keep history of nodes
        #print("get samples = ", (time.time()- s_t))
            ## create genotype array to put into scikit-allel
        # genotypes=samp_ts.genotype_matrix()  ## cant read in ms output for ts sims as no neutral mutations
        # geno_count = genotypes.count_alleles()
    ## calculate Tajima's D 
            ## using tskit branches
        tajdb =  samp_ts.Tajimas_D(sample_sets=None, windows=win, mode="branch") ## using ts branches
            ## using allel
        # taja = allel.windowed_tajima_d(genotypes, windows = win)
        

    ## calculate diversity
            ## using tskit
        div = samp_ts.diversity(sample_sets = None, windows = win)
            ## Collate data
        for w in range(nWin):           ##run through widnows
            #w_t = time.time()
            # print(w)
                ##count of selected and neutral mutations per window (commented to speed up)
            # s_mutcount = np.sum((mut_ud.mut_pos[mut_ud.mut_type=="fluctuating"] >= win[w]) & (mut_ud.mut_pos[mut_ud.mut_type=="fluctuating"] < (win[w+1]-1)))
            # n_mutcount = np.sum((mut_ud.mut_pos[mut_ud.mut_type=="neutral"] >= win[w]) & (mut_ud.mut_pos[mut_ud.mut_type=="neutral"] < (win[w+1]-1)))
           
            dict3={}
            dict3.update({"time":t})                        ## generation
            dict3.update({"n_win":w})                       ## identifier for window
            dict3.update({"win_start" : win[w]})            ## window start position
            dict3.update({"win_end" : win[w+1]-1})          ## window end position
            # dict3.update({"n_s_mut" : s_mutcount})        ## no. selected mutations in window
            # dict3.update({"n_n_mut" : n_mutcount})        ## no. neutral mutatiions per window
            dict3.update({"tajimas_d_branch":tajdb[w]})     ## tajima's D (calculated with tskit)
            dict3.update({"diversity": div[w]})             ## diversity (tajimas pi; calculated with tskit)
            rows_list3.append(dict3)
            #print("per win = ", (time.time()- w_t))
       # print("per time point = ", (time.time()- s_t))
### REMOVE BRACKETS
    #start_time = time.time()   
    s_stats = pd.DataFrame(rows_list3) ## convert dict to df
   # print("convert to pd = ", (time.time()- start_time))
    
    
    #start_time = time.time()   
        ## write statistic df to text file
    s_stats.to_string(buf = "~/oliviaphd/seglift_even_dist/py_out/sim_s_stat_{0}_{1}.txt".format(group, sim_run), index=False)
   # print("pd to txt = ", (time.time()- start_time))
    
    #print("Time for sim run = ", (time.time()- sim_runt), "\n")

#### PLOTS
# plt.scatter(s_stat.n_win[s_stat.Gen==9800.0], s_stat.tajimas_d_branch[s_stat.Gen==9800.0], s=2)
# plt.scatter(n_met.nearest_dist[n_met.mut_type == "neutral"], n_met.allele_freq[n_met.mut_type == "neutral"], s=2)
# plt.show()

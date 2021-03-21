
sim_type = "seglift_long"
##"seglift_win"

import os
os.chdir("/Users/olivia/oliviaphd/{0}".format(sim_type))

import multiprocess as mp
import msprime
import pyslim
import numpy as np
import tskit
import pandas as pd
import time
import itertools
import allel

#groups = (4, 5)
group =4 ## identifier for a set of parameters
runs = 2
#genomeSize = int(1e6)
nChrom = 5
popnSize = int(1e4)
mutRate = 1e-6
recRate = 1e-8
l = 5
nWin = 101 * nChrom
sum_gen = 15#no. summer generations
win_gen = 15#no. winter generations


sim_run = 1

## STAT FUNCTIONS ##


# def FW_H(self):
#     ts = self
#     W = np.zeros((ts.get_num_samples(), 4))

#     def f(x):
#         H = np.sum(x * x), np.prod(x + np.arange(x.shape[0]))
#         return y

#     sigma = ts.general_stat(W, f, 2, ts.get_breakpoints(), mode="branch")
#     assert sigma.shape == (ts.get_num_trees(), 2)




####   LOAD IN  TS
# for group in groups:
#     print("Group:" + str(group))
#     for x in range(runs):
#         #print(x)
        
#         start_time = time.time() 
#         sim_run = str(x+1)
#         print(sim_run)
        
def ts_analysis(sim_run):
        slim_ts = pyslim.load("./group_{0}/treeseq_{1}_{2}_{3}.trees".format(group,sim_type,group,sim_run)).simplify()
        genomeSize = slim_ts.sequence_length
        ## Check number of mutations that were introduced
        if (slim_ts.num_sites != l):
            print("Less than " + str(l) + "introduced mutations")
        else:
            print (str(l) + " introduced mutations")
            
        ## SUMMARISE MUTATIONS
        
        #start_time = time.time()
            ## obtain data for selected sites (used later to remove from tree sequence)
        mut_met = pd.DataFrame({"mut_num": [], "mut_pos": [], "site_id":[]})
        for mut in slim_ts.mutations():## run through muts in ts
            #print(mut)
            mut_met = mut_met.append({"mut_num" : mut.site, "mut_pos" : mut.position, "site_id" : mut.id}, ignore_index=True)
            ## remove replicated mutations (as mutations added in SLiM by randomly sampling individuals, goes against tree structure of burnin hence sometimes have duplicates)
       # t = slim_ts.dump_tables()
       # t.deduplicate_sites()   ### cannot get deduplicate to work
       
        #ts = slim_ts.delete_sites(list(mut_met.site_id.astype(int)))
       
           
        
       #mut_met = mut_met.loc[mut_met.astype(str).drop_duplicates(subset=("mut_num")).index] ## doesnt remove duplicates from tree only table
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
            dict1.update({"time" : slim_ts.slim_generation-ind.time}) ## generation they are from in SLiM time (SLiM and ts count time differently)
            dict1.update({"pop" : ind.population})  ## population individual is from
            dict1.update({"season" : ind_season})   ## season individuals was in
            dict1.update({"nodes": ind.nodes})      ## genotypes individuals contained
        
            rows_list.append(dict1)
        
        ind_met = pd.DataFrame(rows_list) ## convert from dictionary to data frame (fast)
        #print("Time for ind table = ", (time.time()- start_time))
        
        
            
        ## ADD NEUTRAL MUTATIONS
        #start_time = time.time()
        mut_ts = pyslim.SlimTreeSequence(msprime.mutate(slim_ts, rate=mutRate, keep=True))
        ##discrete = True
        #print("Time to ad neutral muts = ", (time.time()- start_time))
        
           # start_time = time.time()
        
        ## SUMMARISE NEUTRAL MUTATIONS
        
        rows_list2 = []
       # a_list = list(mut_met.mut_pos)
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
            if mut.metadata["mutation_list"] == []:
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
        list_m = n_met.mut_num[n_met.mut_type=="fluctuating"].astype(int)
        list_mu = list_m .tolist()
        #mut_met = mut_met.loc[mut_met.astype(str).drop_duplicates(subset=("mut_num")).index]
            ## remove selected loci from ts for analysis of neutral loci
        #stat_ts = pyslim.SlimTreeSequence(mut_ts.delete_sites(list_mu, record_provenance=True))
        stat_ts = mut_ts.delete_sites(list_mu, record_provenance=True)
          
        ## CALCULATE SUMMARY STATISTICS
        
        rows_list3 = []
        
        #start_time = time.time()
            ## create windows for stats to be calculated in
        win = np.linspace(0, stat_ts.sequence_length, num=nWin+1).astype(int)
        #print("create windows = ", (time.time()- start_time))
        start_time = time.time()
            ## seperate into ts for each timepoint
        for t in np.unique(ind_met.time):   ## run through time points
           # s_t = time.time()
                ## nodes of indviduals at time t
            sample = ind_met.nodes[ind_met.time == t]
                ## convert to list
            samples= list(itertools.chain(*sample))
                ## create ts of just samples at time t using list
                ## used by tskit to calculate states based on tree structure
            #samp_ts = stat_ts.subset(samples)              ## may not be keeping history of nodes
            samp_ts = stat_ts.simplify(samples = samples)   ## should keep history of nodes
            #print("get samples = ", (time.time()- s_t))
                ## create genotype array to put into scikit-allel  shape ==(var, ind)
            samp_gm=samp_ts.genotype_matrix()  
                ## convert samp_gm to shape(ind, var) ## NOT REQUIRED
            # samp_g2h = []
            # for var in range(len(samp_gm[1,])):
            #    #print(var)
            #    samp_g2h.append(samp_gm[:,var])
            #     #print(samp_gm[:,var])
            #samp_hm = np.asarray(samp_gm)
            
            
            ##Check all sites of gm are polymorphic
            for i in range(len(samp_gm)):
                if sum(samp_gm[i]) == 0:
                    print (str(i) + "is not polymoprhic")
                    
                    
            ## convert to haplotyoe array
            h= allel.HaplotypeArray(samp_gm)
            samp_ac = h.count_alleles()
         
            #h_stat = allel.moving_garud_h(h, size = int(genomeSize/nWin), step=None)
            
            size_factor = 1e15 # HARD coded
            positions = []
            for mut in samp_ts.mutations():
                positions.append(int(mut.position*size_factor))
                
            # list comprehension
            positions = [mut.position * size_factor for mut in samp_ts.mutations()]
                
            mut_positions = np.asarray(positions)  ##rounding mutation positions
            
            #if len(np.unique(mut_positions.astype(int))) != len(mut_positions):
               # print("Multiple mutattions at single pos")
            
            hap_stats = allel.windowed_statistic(mut_positions,h,allel.garud_h,size=int((genomeSize/nWin)*size_factor), start = 0, stop = int((genomeSize-1)*size_factor))
            
        
        ## calculate Tajima's D 
                ## using tskit branches
            tajdb =  samp_ts.Tajimas_D(sample_sets=None, windows=win, mode="branch") ## using ts branches
                ## using allel
            # taja = allel.windowed_tajima_d(genotypes, windows = win)
            
            theta_w= allel.windowed_statistic(mut_positions, (mut_positions, samp_ac), allel.watterson_theta, size=int((genomeSize/nWin)*size_factor), start = 0, stop = int((genomeSize-1)*size_factor))
            
            #theta_w = allel.windowed_watterson_theta(mut_positions, samp_ac,size = int(genomeSize/nWin)*1e6, start = 0) 
        
        ## calculate diversity
                ## using tskit
            div = samp_ts.diversity(sample_sets = None, windows = win)
                ## Collate data
            for w in range(nWin):           ##run through widnows
               # w_t = time.time()
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
                dict3.update({"theta_w": theta_w[0][w]})        ## watterson's theta
                dict3.update({"H1": hap_stats[0][w,0]})              ## H1
                dict3.update({"H12":hap_stats[0][w,1]})             ## H12
                dict3.update({"H123": hap_stats[0][w,2]})            ## H123
                dict3.update({"H2/H1": hap_stats[0][w,3]})           ## H2/H1
              
                rows_list3.append(dict3)
               # print("per win = ", (time.time()- w_t))
            #print("per time point = ", (time.time()- s_t))
        ### REMOVE BRACKETS
        #start_time = time.time()   
        s_stats = pd.DataFrame(rows_list3) ## convert dict to df
           # print("convert to pd = ", (time.time()- start_time))
        
        
        #start_time = time.time()   
            ## write statistic df to text file
        s_stats.to_string(buf = "./group_{0}/sim_s_stat_2sps_{1}_{2}.txt".format(group,group,sim_run), index=False)
        print("run= ", (time.time()- start_time))
        
        #print("Time for sim run = ", (time.time()- sim_runt), "\n")
    

processes=[]

for run in range(runs):
    sim_run = run+1
    print(sim_run)
    p=mp.Process(target=ts_analysis, args=[sim_run])
    p.start()
    processes.append(p)
    
for process in processes:
    process.join()




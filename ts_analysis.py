
# sim_type = "seglift_long"
# ##"seglift_win"

# import os
# os.chdir("/Users/olivia/oliviaphd/{0}".format(sim_type))

import multiprocess as mp
import msprime
import pyslim
import numpy as np
import tskit
import pandas as pd
import time
import itertools
import allel

# #groups = (4, 5)
# group =9 ## identifier for a set of parameters
# runs = 8
# #genomeSize = int(1e6)
# nChrom = 21
# popnSize = int(1e4)
# mutRate = 3e-9
# recRate = 1e-8 ##change to 2e-8 in future sims
# l = 10
# nWin = 51 * nChrom
# sum_gen = 15#no. summer generations
# win_gen = 15#no. winter generations


# # sim_run = 1

      
def analyse(tmpdir, group, sim_type, sim_run, mutRate, l, nChrom, nWin, sum_gen, win_gen):
    ## input  group(parameter identifier), simulation type, simulation run), mutation rate, 
    ##  number of selected loci, number of chromosomes and number of windows
    
 ## INPUT DATA
         # read in treesequence (ts) generated in SLiM
    slim_ts = pyslim.load("{2}/treeseq_{0}_{1}.trees".format(group,sim_run, tmpdir)).simplify()
        # extract the length of the simulate seqeunce from slim_ts
    genomeSize = slim_ts.sequence_length
        # check number of mutations that were introduced in slim simulation
    if (slim_ts.num_sites != l):
        print("Less than " + str(l) + "introduced mutations")
    else:
        print (str(l) + " introduced mutations")
        
## SUMMARISE MUTATIONS - obtain data for selected sites (used later to remove from tree sequence)

        #set up pd dataframe to store metadata for mutations generated in slim (selected only)
    mut_met = pd.DataFrame({"mut_site": [], "mut_pos": [], "mut_id":[]})
    ## run through mutations (muts) in ts
    for mut in slim_ts.mutations():
        #print(mut)
        mut_met = mut_met.append({"mut_site" : mut.site, "mut_pos" : slim_ts.site(mut.site).position, "mut_id" : mut.id}, ignore_index=True)

## SUMMARISE INDIVIDUALS - obtain metadata for individuals 'remembered' in ts
    rows_list = []
        # run through individuals (inds) in ts
    for ind in slim_ts.individuals():
        #print(ind)
            # determine what season ind was saved in
        if ind.time % (sum_gen+win_gen) == 0 or ind.time % (sum_gen+win_gen) > sum_gen:  
            ind_season = "W"
        else:
            ind_season = "S"
    
        dict1 = {}
            # individual's id in ts
        dict1.update({"id" : ind.id})
            # generation ind is from in SLiM time (SLiM and ts count time differently)
        dict1.update({"time" : slim_ts.slim_time(ind.time)}) 
            # population individual is from
        dict1.update({"pop" : ind.population})  
            # season individuals was remembered in
        dict1.update({"season" : ind_season})
            # genotypes individuals contained (2 nodes, each directing to a haploid genotype)
        dict1.update({"nodes": ind.nodes})      
    
        rows_list.append(dict1)
        # convert from dictionary to data frame (fast)
    ind_met = pd.DataFrame(rows_list) 
    
## REMOVE SEELCTED MUTATIONS - remove from ts so will not interfere with statistics
    ##no_mut_ts = slim_ts.delete_sites(list(mut_met.mut_site.astype(int))) 
    
## ADD NEUTRAL MUTATIONS - simulations run without neutral mutations, need to put on tree to generate summary statistics removes current muts on tree
    mut_ts = pyslim.SlimTreeSequence(msprime.sim_mutations(slim_ts, rate=mutRate, keep=False))

# ## SUMMARISE NEUTRAL MUTATIONS - obtain data for neutral mutations added to ts
    
#     rows_list2 = []
#         # cycle through neutral mutations on ts
#     for mut in mut_ts.mutations():
#         #print(mut)
  
#         dict2 = {}
#              # id of mutation in ts
#         dict2.update({"mut_id":mut.id})              
#             # site of mutation (1 to l)
#         dict2.update({"mut_site" : mut.site})  
#             # position of mutation              
#         dict2.update({"mut_pos" : mut.position})            
    
#         rows_list2.append(dict2)
#     n_met = pd.DataFrame(rows_list2) ## convert dict to df
    
      
## CALCULATE SUMMARY STATISTICS
    
    rows_list3 = []
    # rows_list4 = []

        # create windows for stats to be calculated in
        # windows for tskit statistics
    win3 = np.linspace(0, mut_ts.sequence_length, num=nWin+1).astype(int)
    # alwin = int(nWin/4)
    # win4 = np.linspace(0, mut_ts.sequence_length, num=alwin+1).astype(int)
    #     # windoes for scikit.allel statistics
        
    # al_win4 = []
    # for w in range(len(win4)-1):
    #      if w == nWin:
    #          window = [win4[w], int(genomeSize)]
    #      else:
    #          window = [win4[w], win4[w+1]-1]
    #      al_win4.append(window)


        # cycle throught timepoints for which data has been collected
    for t in np.unique(ind_met.time): 
        
            # collate nodes (geotype identifiers) of indviduals at time t
        sample = ind_met.nodes[ind_met.time == t]
        
            ## convert to list
        samples= list(itertools.chain(*sample))
        
            ## create ts of just samples at time t using list
        samp_ts = mut_ts.simplify(samples = samples)
        
            ## create genotype array to put into scikit-allel  shape ==(var, ind)
        samp_gm=samp_ts.genotype_matrix()  
                
            # convert genotype matrix to haplotyoe array for haplotype statistics
        h= allel.HaplotypeArray(samp_gm)
            # allele count for scikit.allel stats
        samp_ac = h.count_alleles()
            # positions of mutations in samp_ts for scikit.allel windowed_statistic function
        mut_positions = [mut.site for mut in samp_ts.mutations()]
        
            # generate haplotype statistics (H1, H12, H123, H2/H1)
        hap_stats = allel.windowed_statistic(mut_positions,h,allel.garud_h, windows = win3)
        
            # tajimas D using tskit and branches of ts
        tajdb =  samp_ts.Tajimas_D(sample_sets=None, windows=win3, mode="branch")
            
            # wattersons theta using scikit.allel
        theta_w= allel.windowed_statistic(mut_positions, (mut_positions, samp_ac), allel.watterson_theta, windows = win3)
        
        # tajimas D using scikit.allel
        tajda= allel.windowed_statistic(mut_positions, (mut_positions, samp_ac), allel.tajima_d, windows = win3)

            # calculate diversity (tajima's pi) using tskit
        div = samp_ts.diversity(sample_sets = None, windows = win3)  ##fix windows
        
        
            # check that all stats have ben calculated over the correct number of windows
        ts_tests = [div, tajdb]
        # al_tests = [theta_w[0],hap_stats[0]]
        for test in ts_tests:
            if len(test)!= nWin:
                print("error in ts test ", test, ", number of values does not match number of windows")
                
        # for test in al_tests:
        #     if len(test)!= alwin:
        #         print("error in al test ", test, ", number of values does not match number of windows")
                
        ## Collate summary statics into dataframe
            # loop over windows
        for w in range(nWin): 
            
            try:
                h1 = hap_stats[0][w][0]
            except TypeError:
                h1="NaN"
            try:
                h12 = hap_stats[0][w][1]
            except TypeError:
                h12 = "NaN"
            try:
                h123 = hap_stats[0][w][2]
            except TypeError:
                h123 = "NaN"
            try:
                h2h1 = hap_stats[0][w][3]
            except TypeError:
                  h2h1="NaN"
           
            dict3={}
            dict3.update({"time":t})                        ## generation
            dict3.update({"n_win":w})                       ## identifier for window
            dict3.update({"win_start" : win3[w]})            ## window start position
            dict3.update({"win_end" : win3[w+1]-1})          ## window end position
            dict3.update({"chrom" : int(w/(nWin/nChrom))+1})  ## chromosome
            dict3.update({"tajimas_d_branch":tajdb[w]})     ## tajima's D (calculated with tskit)
            dict3.update({"diversity": div[w]})             ## diversity (tajimas pi; calculated with tsk
            dict3.update({"tajimas_d_allel": tajda[0][w]})        ## watterson's theta
            dict3.update({"theta_w": theta_w[0][w]})        ## watterson's theta
            dict3.update({"H1": h1})         ## H1
            dict3.update({"H12":h12})         ## H12
            dict3.update({"H123": h123})       ## H123
            dict3.update({"H2H1": h2h1})       ## H2/H1
          
            rows_list3.append(dict3)
            
            # convert dictionary to datafram
    ts_stats = pd.DataFrame(rows_list3) 
    
            # write statistic df to text file
    ts_stats.to_string(buf = "{2}/sim_ts_stat_{0}_{1}.txt".format(group,sim_run,tmpdir), index=False)
    
    
        
    # al_stats = pd.DataFrame(rows_list4) 
    
            # write statistic df to text file
    # al_stats.to_string(buf = "./sim_al_stat_{0}_{1}.txt".format(group,sim_run), index=False)
    
        

# processes=[]
# start_time = time.time()
# for run in range(runs):
#     sim_run = run
#     print(sim_run)
#     p=mp.Process(target=ts_analysis, args=[group, sim_type, sim_run, mutRate, l, nChrom, nWin])
#     p.start()
#     processes.append(p)
    
# for process in processes:
#     process.join()

# print("Time = ", (time.time() - start_time))


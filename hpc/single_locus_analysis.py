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
import itertools
import statistics
sys.path.insert(1, '/hpcfs/users/a1704225/scripts/oliviaphd/')
import single_locus_hpc

params=sys.argv[1]
sim_run = sys.argv[2]
results_dir =str(sys.argv[3])
####  READ IN PARAMETERS
    # load in parameter file
with open('/hpcfs/users/a1704225/parameters/single_locus/{0}.txt'.format(params), 'r') as f:
    parameters = yaml.load(f, Loader=yaml.FullLoader)

    #set parameters from file
genomeSize = int(parameters["genomeSize"])
recRate = parameters["recRate"]
mutRate=parameters["mutRate"]
s_pop = int(parameters["s_pop"])
w_pop = int(parameters["w_pop"])
h_s = parameters["h_s"]
h_w = parameters["h_w"]
s_s = parameters["s_s"]
s_w = parameters["s_w"]
rGen=int(parameters["rGen"])
fitness_on = parameters["fitness_on"]
sum_gen = int(parameters["sum_gen"])
win_gen = int(parameters["win_gen"])
nWin = parameters["winpChrom"]
group=parameters["group"]
freq=int(parameters["f"])
burnin_Ne = int(parameters["burnin_Ne"])



start_time = time.time()
 ## INPUT DATA
         # read in treesequence (ts) generated in SLiM
slim_ts = tskit.load("{0}/treeseq_group_{1}_{2}.trees".format(results_dir,group,sim_run)).simplify()
    # extract the length of the simulate seqeunce from slim_ts
    # check number of mutations that were introduced in slim simulation
if (slim_ts.num_sites != 1):
    print("ERROR? " + str(slim_ts.num_sites) + "introduced mutations")
else:
    print (str(1) + " introduced mutations")

## SUMMARISE INDIVIDUALS - obtain metadata for individuals 'remembered' in ts
rows_list = []
    # run through individuals (inds) in ts
for ind in slim_ts.individuals():
    #print(ind)
    dict1 = {}
        # individual's id in ts
    dict1.update({"id" : ind.id})

    #     # genotypes individuals contained (2 nodes, each directing to a haploid genotype)
    dict1.update({"nodes": ind.nodes})

    rows_list.append(dict1)
    # convert from dictionary to data frame (fast)
ind_met = pd.DataFrame(rows_list)

ind_times = np.unique(slim_ts.individual_times).astype(int)


## ADD NEUTRAL MUTATIONS - simulations run without neutral mutations, need to put on tree to generate summary statistics removes current muts on tree
mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)


## CALCULATE SUMMARY STATISTICS

rows_list3 = []

    # create windows for stats to be calculated in
    # windows for tskit statistics

win3 = np.linspace(0, mut_ts.sequence_length, num=nWin+1).astype(int)
   # windows for scikit.allel statistics

al_win3 = []
for w in range(len(win3)-1):
      if w == nWin:
          window = [win3[w], int(mut_ts.sequence_length)]
      else:
          window = [win3[w], win3[w+1]-1]
      al_win3.append(window)


    # cycle throught timepoints for which data has been collected
for t in ind_times:

        # collate nodes (geotype identifiers) of indviduals at time t
    sample_ind = pyslim.individuals_alive_at(slim_ts, t)
    sample = ind_met["nodes"].iloc[sample_ind]

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
    mut_positions = [int(var.position+1) for var in samp_ts.variants()]
    
        ## crete genotype array for LD
    odds = h[:,::2]
    evens = h[:,1::2]
    
    gn=odds+evens
    gn=gn.view('int8')
            
    ehh=allel.ehh_decay(h)
    
    ehh_win=allel.windowed_statistic(mut_positions, ehh, statistics.mean, windows=al_win3)
    
    r2=allel.windowed_r_squared(mut_positions, gn, windows=al_win3)
    
        # generate haplotype statistics (H1, H12, H123, H2/H1)
    hap_stats = allel.windowed_statistic(mut_positions,h,allel.garud_h, windows = al_win3)

    hap_div = allel.windowed_statistic(mut_positions,h,allel.haplotype_diversity, windows = al_win3)
    
        # tajimas D using tskit and branches of ts
    tajdb =  samp_ts.Tajimas_D(sample_sets=None, windows=win3, mode="site")

    ## no. segregatiig sites
    n_seg =samp_ts.segregating_sites(sample_sets=None, windows=win3, mode="site", span_normalise=False)

        # wattersons theta using scikit.allel
    tw_a= allel.windowed_watterson_theta(mut_positions, samp_ac, windows=al_win3)

    # tajimas D using scikit.allel
    tajda= allel.windowed_tajima_d(mut_positions, samp_ac, windows=al_win3)

        # calculate diversity (tajima's pi) using tskit
    div = samp_ts.diversity(sample_sets = None, windows = win3)  ##fix windows
    ts_tests = [div, tajdb, r2[0]]
    for test in ts_tests:
        if len(test)!= nWin:
            print("error in test ", test, ", number of values does not match number of windows")

    # for test in al_tests:
    #     if len(test)!= alwin:
    #         print("error in al test ", test, ", number of values does not match number of windows")

    ## Collate summary statics into dataframe
        # loop over windows
    for w in range(nWin):


        dict3={}
        dict3.update({"time":slim_ts.slim_time(t)})                        ## generation
        dict3.update({"n_win":w})                       ## identifier for window
        dict3.update({"win_start" : win3[w]})            ## window start position
        dict3.update({"win_end" : win3[w+1]-1})          ## window end position
        dict3.update({"seg_sites" : n_seg[w]})   ## number fo segregating sites
        dict3.update({"diversity": div[w]})             ## diversity (tajimas pi; calculated with tsk
        dict3.update({"tajimas_d_branch":tajdb[w]})     ## tajima's D (calculated with tskit)
        dict3.update({"tajimas_d_allel": tajda[0][w]})        
        dict3.update({"theta_w_allele": tw_a[0][w]})        ## watterson's theta (allele with scikit allel)
        dict3.update({"ehh": ehh_win[0][w]})        ## mean ehh per window
        dict3.update({"r2": r2[0][w]})        ## r2 per win
        dict3.update({"haplotype_diversity": hap_div[0][w]})        ## haplotype diversity
        dict3.update({"H1": hap_stats[0][w][0]})         ## H1
        dict3.update({"H12":hap_stats[0][w][1]})         ## H12
        dict3.update({"H123": hap_stats[0][w][2]})       ## H123
        dict3.update({"H2H1": hap_stats[0][w][3]})       ## H2/H1

        rows_list3.append(dict3)

        # convert dictionary to datafram
ts_stats = pd.DataFrame(rows_list3)

        # write statistic df to text file
ts_stats.to_string(buf = "{2}/sim_stat_{0}_{1}.txt".format(group,sim_run,results_dir), index=False)    
   
print("Time = ", (time.time() - start_time))

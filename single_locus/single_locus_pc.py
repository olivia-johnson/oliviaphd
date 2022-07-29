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


## run on computer
group=12
sys.path.insert(1,"/Users/olivia/oliviaphd/hpc/") 
import single_locus_hpc
import multiprocess as mp
for group in range(29, 37):
    with open('/Users/olivia/phd_data/hpc_parameters/single_locus/group_{0}.txt'.format(group), 'r') as f:
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
    winpChrom = parameters["winpChrom"]
    group=parameters["group"]
    burnin_Ne = int(parameters["burnin_Ne"])
    
   
    tmpdir="/Users/olivia/phd_data/Results/single_locus/group_{0}".format(group)
    os.mkdir(tmpdir)
    runs=10
     
    # BURNIN
    processes=[]
    for sim_run in range(runs):
        p=mp.Process(target=single_locus_hpc.single_locus_burnin, args=[tmpdir, group, sim_run, genomeSize, s_pop, burnin_Ne, recRate])
        p.start()
        processes.append(p)
        
    for process in processes:
        process.join()
    
    tmpdir="~/phd_data/Results/single_locus/group_{0}".format(group)
    results_dir="~/phd_data/Results/single_locus/group_{0}".format(group)
    processes=[]
    
    for sim_run in range(runs):
        #sim_run=sim_run+(runs)
        p=mp.Process(target=single_locus_hpc.simulate_single_locus, args=[tmpdir, results_dir, group, sim_run, recRate, genomeSize, s_pop, w_pop, h_s, h_w, s_s, s_w, rGen, fitness_on, sum_gen, win_gen])
        p.start()
        processes.append(p)
        
    for process in processes:
        process.join()
    
    
        
    
    #### ANALYSE TREE SEQUENCE WITH TS_ANALYSIS
    
    tmpdir="/Users/olivia/phd_data/Results/single_locus/group_{0}".format(group)
    nWin=winpChrom
    
    for sim_run in range(0, runs): 
        start_time = time.time()
        # sim_run= sim_run +1
     ## INPUT DATA
             # read in treesequence (ts) generated in SLiM
        slim_ts = pyslim.SlimTreeSequence.load("{0}/treeseq_group_{1}_{2}.trees".format(tmpdir,group,sim_run)).simplify()
            # extract the length of the simulate seqeunce from slim_ts
            # check number of mutations that were introduced in slim simulation
        if (slim_ts.num_sites != 1):
            print("ERROR? " + str(slim_ts.num_sites) + "introduced mutations")
        else:
            print (str(1) + " introduced mutations")

    ## SUMMARISE MUTATIONS - obtain data for selected sites (used later to remove from tree sequence)

            #set up pd dataframe to store metadata for mutations generated in slim (selected only)
        # mut_met = pd.DataFrame({"mut_site": [], "mut_pos": [], "mut_id":[]})
        # ## run through mutations (muts) in ts
        # for mut in slim_ts.mutations():
        #     #print(mut)
        #     mut_met = mut_met.append({"mut_site" : mut.site, "mut_pos" : slim_ts.site(mut.site).position, "mut_id" : mut.id}, ignore_index=True)

    ## SUMMARISE INDIVIDUALS - obtain metadata for individuals 'remembered' in ts
        rows_list = []
            # run through individuals (inds) in ts
        for ind in slim_ts.individuals():
            #print(ind)
            #     # determine what season ind was saved in
            # if ind.time % (sum_gen+win_gen) == 0 or ind.time % (sum_gen+win_gen) > sum_gen:
            #     ind_season = "W"
            # else:
            #     ind_season = "S"

            dict1 = {}
                # individual's id in ts
            dict1.update({"id" : ind.id})
                # generation ind is from in SLiM time (SLiM and ts count time differently)
           # dict1.update({"time" : slim_ts.slim_time(ind.time)})
                # population individual is from
            # dict1.update({"pop" : ind.population})
            #     # season individuals was remembered in
            # dict1.update({"season" : ind_season})
            #     # genotypes individuals contained (2 nodes, each directing to a haploid genotype)
            dict1.update({"nodes": ind.nodes})

            rows_list.append(dict1)
            # convert from dictionary to data frame (fast)
        ind_met = pd.DataFrame(rows_list)

        ind_times = np.unique(slim_ts.individual_times).astype(int)

    ## REMOVE SEELCTED MUTATIONS - remove from ts so will not interfere with statistics
        ##no_mut_ts = slim_ts.delete_sites(list(mut_met.mut_site.astype(int)))

    ## ADD NEUTRAL MUTATIONS - simulations run without neutral mutations, need to put on tree to generate summary statistics removes current muts on tree
        mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)
        mut_ts = msprime.sim_mutations(slim_ts,discrete_genome=False, rate=mutRate, keep=False )

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
            sample_ind = slim_ts.individuals_alive_at(t)
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
            mut_positions = [int(var.position) for var in samp_ts.variants()]
            
                ## crete genotype array for LD
            odds = h[:,::2]
            evens = h[:,1::2]
            
            gn=odds+evens
            ld_calc = tskit.LdCalculator(samp_ts)
            r2 = ld_calc.r2_array(50000, direction=1)
            
            allel.windowed_r_squared(mut_positions, gn, windows=al_win3)
            
            r2_win_val, r2_win, r2_win_n = allel.windowed_statistic(mut_positions,gn,allel.rogers_huff_r, windows = al_win3)
            
            
                # generate haplotype statistics (H1, H12, H123, H2/H1)
            hap_stats = allel.windowed_statistic(mut_positions,h,allel.garud_h, windows = al_win3)

            hap_div = allel.windowed_statistic(mut_positions,h,allel.haplotype_diversity, windows = al_win3)
                # tajimas D using tskit and branches of ts
            tajdb =  samp_ts.Tajimas_D(sample_sets=None, windows=win3, mode="site")

                # wattersons theta using ts
            # tw_b= theta_w(samp_ts, win3)
        

                # wattersons theta using scikit.allel
            tw_a= allel.windowed_watterson_theta(mut_positions, samp_ac, windows=al_win3)

            # tajimas D using scikit.allel
            tajda= allel.windowed_tajima_d(mut_positions, samp_ac, windows=al_win3)

                # calculate diversity (tajima's pi) using tskit
            div = samp_ts.diversity(sample_sets = None, windows = win3)  ##fix windows

                # check that all stats have ben calculated over the correct number of windows
            ts_tests = [div, tajdb]
            # al_tests = [theta_w[0],hap_stats[0]]
            for test in ts_tests:
                if len(test)!= nWin:
                    print("error in test ", test, ", number of values does not match number of windows")

            # for test in al_tests:
            #     if len(test)!= alwin:
            #         print("error in al test ", test, ", number of values does not match number of windows")

            ## Collate summary statics into dataframe
                # loop over windows
            for w in range(nWin):

                # h1 = hap_stats[0][w][0]

                # h12 = hap_stats[0][w][1]

                # h123 = hap_stats[0][w][2]

                # h2h1 = hap_stats[0][w][3]


                # try:
                #     h1 = hap_stats[0][w][0]
                # except TypeError:
                #     h1="NaN"
                # try:
                #     h12 = hap_stats[0][w][1]
                # except TypeError:
                #     h12 = "NaN"
                # try:
                #     h123 = hap_stats[0][w][2]
                # except TypeError:
                #     h123 = "NaN"
                # try:
                #     h2h1 = hap_stats[0][w][3]
                # except TypeError:
                #       h2h1="NaN"

                dict3={}
                dict3.update({"time":slim_ts.slim_time(t)})                        ## generation
                dict3.update({"n_win":w})                       ## identifier for window
                dict3.update({"win_start" : win3[w]})            ## window start position
                dict3.update({"win_end" : win3[w+1]-1})          ## window end position
                # dict3.update({"chrom" : int(w/(nWin/nChrom))+1})  ## chromosome
                dict3.update({"tajimas_d_branch":tajdb[w]})     ## tajima's D (calculated with tskit)
                dict3.update({"diversity": div[w]})             ## diversity (tajimas pi; calculated with tsk
                dict3.update({"tajimas_d_allel": tajda[0][w]})        
                # dict3.update({"theta_w_branch": tw_b[w]})        ## watterson's theta (branch with tskit)
                dict3.update({"theta_w_allele": tw_a[0][w]})        ## watterson's theta (allele with scikit allel)
                dict3.update({"H1": hap_stats[0][w][0]})         ## H1
                dict3.update({"H12":hap_stats[0][w][1]})         ## H12
                dict3.update({"H123": hap_stats[0][w][2]})       ## H123
                dict3.update({"H2H1": hap_stats[0][w][3]})       ## H2/H1

                rows_list3.append(dict3)

                # convert dictionary to datafram
        ts_stats = pd.DataFrame(rows_list3)

                # write statistic df to text file
        ts_stats.to_string(buf = "{2}/sim_stat_{0}_{1}.txt".format(group,sim_run,tmpdir), index=False)

        print("Time = ", (time.time() - start_time))
    
    
    processes=[]
    start_time = time.time()
    for run in range(runs):
        sim_run = run
        print(sim_run)
        p=mp.Process(target=single_locus_hpc.analyse, args=[tmpdir, group, sim_run, genomeSize, s_pop, burnin_Ne, recRate, sum_gen, win_gen])
        p.start()
        processes.append(p)
        
    for process in processes:
        process.join()
    
    print("Time = ", (time.time() - start_time))

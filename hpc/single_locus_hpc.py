import os
import msprime
import pyslim
import numpy as np
import tskit
import pandas as pd
import time
import itertools
import allel
import statistics


def recombination_map(tmpdir, group, genomeSize, recRate, winSize):
    ## input group (parameter set identifier, number of chromosomes, chromosome size, recombination rate)

    ## GENERATE CHROMOSMES ##
## create recombination map for msprime and slim to simulate unlinked section
    unlinked_region = 10*winSize    
    rec_pos=[0, genomeSize-1,genomeSize, (genomeSize+unlinked_region)]
    rates=[recRate, 0.5, recRate]
    rate_map = msprime.RateMap(position=rec_pos,rate=rates)
    rate_map
    sequenceSize=rate_map.sequence_length
    print("Sequence length is ", str(sequenceSize))
## generate slim recombination map
        # reformat to comply with recombinate map required in slim
    slim_rates=[recRate,recRate, 0.5, recRate]
    slim_pos=[0, genomeSize-1,genomeSize, (genomeSize+unlinked_region)-1]
    rec_data = {'ends': slim_pos, 'rates': slim_rates}
    slim_rec = pd.DataFrame(rec_data)
## output slim recombiation map to text file to be read into forward slim simulation
    slim_rec.to_csv("{0}/rec_map_group_{1}.txt".format(tmpdir, group), index=False, header = True, sep = "\t")

    return rate_map, sequenceSize


def single_locus_burnin(tmpdir, group, sim_run, sequenceSize, s_pop, burnin_Ne, rate_map):

    ## COALESCENT BURN IN
    start_time = time.time()
    print("Burnin Ne is ", str(burnin_Ne))

    ##daiquiri.setup(level="DEBUG") ##debug
    burnin = msprime.sim_ancestry(samples=s_pop, population_size=burnin_Ne, recombination_rate=rate_map, sequence_length=sequenceSize)
    ##check burnin size = genome size
    burn_length =burnin.get_sequence_length()
    if burn_length!=sequenceSize:
        print("Burnin sequence length not equal to sequence size!")
    burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
    burnin_ts.dump("{0}/burnin_seglift_group_{1}_{2}.trees".format(tmpdir,group,sim_run))
    print("Time for burnin = ", (time.time()- start_time))
    
def simulate_single_locus(tmpdir, results_dir, group, sim_run, recRate, genomeSize, s_pop, w_pop, h_s, h_w, s_s, s_w, rGen, fitness_on, sum_gen, win_gen, f, sim_type, winSize):

    ## FORWARD SIMULATION
    # for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
    start_time = time.time()
    tmpdir_call = "tmpdir='" + str(tmpdir)+ "'"
    results = "results_dir='" + str(results_dir)+ "'"
    sequnceSize=genomeSize+(10*winSize)
    cmd = 'slim -d "' +str(tmpdir_call)+ '" -d "' +str(results)+ '" -d fit='+ str(fitness_on)+" -d group=" + str(group) + " -d winSize=" + str(winSize) +" -d f=" + str(f) + " -d h_s=" + str(h_s)+ " -d h_w=" + str(h_w)+ " -d s_s=" + str(s_s)+ " -d s_w=" + str(s_w)+" -d g_s=" + str(sum_gen)+" -d g_w=" + str(win_gen)+" -d sim_run=" + str(sim_run) + " -d GenomeSize=" + str(int(sequnceSize)) +  " -d n_s=" + str(int(s_pop)) + " -d n_w=" + str(int(w_pop)) + " -d mut=0.0 -d rr=" + str(recRate) +   " -d rGen="+ str(rGen) +" ~/oliviaphd/hpc/"+str(sim_type)+"_single_locus.slim"

    print(cmd)
    os.system(cmd)
    # #print("Time for SLiM sim = ", (time.time()- start_time))
    print("Simulations took ", (time.time()-start_time) ,  " seconds")
    
def analyse(tmpdir, group, sim_run, mutRate, l, genomeSize, nWin, sum_gen, win_gen):
    ## input  group(parameter identifier), simulation type, simulation run), mutation rate,
    ##  number of selected loci, number of chromosomes and number of windows


    # def Theta_W(ts, sample_sets=None, windows=None, mode="site"):
    #         """
    #         Computes Theta W of sets of nodes from ``sample_sets`` in windows.
    #         Please see the :ref:`one-way statistics <sec_stats_sample_sets_one_way>`
    #         section for details on how the ``sample_sets`` argument is interpreted
    #         and how it interacts with the dimensions of the output array.
    #         See the :ref:`statistics interface <sec_stats_interface>` section for details on
    #         :ref:`windows <sec_stats_windows>`, :ref:`mode <sec_stats_mode>`,
    #         and :ref:`return value <sec_stats_output_format>`.
    #         Operates on ``k = 1`` sample sets at a
    #         time. For a sample set ``X`` of ``n`` nodes, if ``S`` is the number of
    #         sites segregating in ``X`` (computed with :meth:`segregating sites
    #         <.TreeSequence.segregating_sites>`, respectively, both not span
    #         normalised), then Theta W is

    #         .. code-block:: python

    #             D = (T - S / h) / sqrt(a * S + (b / c) * S * (S - 1))
    #             a = 1 + 1 / 2 + ... + 1 / (n - 1)

    #         :param list sample_sets: A list of lists of Node IDs, specifying the
    #             groups of nodes to compute the statistic with.
    #         :param list indexes: A list of 2-tuples, or None.
    #         :param list windows: An increasing list of breakpoints between the windows
    #             to compute the statistic in.
    #         :param str mode: A string giving the "type" of the statistic to be computed
    #             (defaults to "site").
    #         :return: A ndarray with shape equal to (num windows, num statistics).
    #         """
    #         def tw_func(sample_set_sizes, flattened, **kwargs):
    #             n = sample_set_sizes
    #             S = ts.ll_tree_sequence.segregating_sites(n, flattened, **kwargs)
    #             a = np.array([np.sum(1 / np.arange(1, nn)) for nn in n])
    #             with np.errstate(invalid="ignore", divide="ignore"):
                    
    #                 w = S / a ## / #np.sqrt(a * S + (b / (h ** 2 + g)) * S * (S - 1))
    #             return w

    #         return ts.sample_count_stat(
    #             [samples],
    #             tw_func,
    #             1,
    #             sample_sets, windows=windows, mode=mode, span_normalise=False
    #         )
   
        
        
   # def theta_w(ts, wins):
   #     samples = ts.samples()
   #     def f(x):
   #         n = samples
   #         S = ts.segregating_sites()
   #         a = np.array([np.sum(1 / np.arange(1, nn)) for nn in n])
   #         w = S / a 
   #         return w
   #     return ts.sample_count_stat(
   #         [samples],
   #         f,
   #         1,
   #         windows=wins,
   #         polarised=False,
   #         mode="site",
   #         span_normalise=False,
       # )
   
    # def beta(ts, wins):
    #     samples = ts.samples()
      
 ## INPUT DATA
         # read in treesequence (ts) generated in SLiM
    ts=tskit.load("{0}/treeseq_group_{1}_{2}.trees".format(tmpdir,group,sim_run)).simplify()
    slim_ts = pyslim.update(ts) ##update ts from slim 3.7
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
        ehh, ehh_windows, ehh_vars=allel.windowed_statistic(mut_positions, h, allel.ehh_decay, windows=al_win3)
        
        ehh_win=allel.windowed_statistic(mut_positions, ehh, statistics.mean, windows=al_win3)
        
        r2=allel.windowed_r_squared(mut_positions, gn, windows=al_win3)
        
        # r2_win_val, r2_win, r2_win_n = allel.windowed_statistic(mut_positions,gn,allel.rogers_huff_r, windows = al_win3)
        
            # generate haplotype statistics (H1, H12, H123, H2/H1)
        hap_stats = allel.windowed_statistic(mut_positions,h,allel.garud_h, windows = al_win3)

        hap_div = allel.windowed_statistic(mut_positions,h,allel.haplotype_diversity, windows = al_win3)
        
        # ihs = allel.windowed_statistic(mut_positions,(h, mut_positions),allel.ihs, windows = al_win3)
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
              dict3.update({"time":pyslim.slim_time(slim_ts,t)})                        ## generation
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
    ts_stats.to_string(buf = "{2}/sim_stat_{0}_{1}.txt".format(group,sim_run,tmpdir), index=False)




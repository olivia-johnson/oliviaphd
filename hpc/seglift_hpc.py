import os
import msprime
import pyslim
import numpy as np
import tskit
import pandas as pd
import time
import itertools
import allel

##to debug burnin
##import daiquiri

def recombination_map(tmpdir, group, l, nChrom, chromSize, recRate):
    ## input group (parameter set identifier, number of chromosomes, chromosome size, recombination rate)

    ## GENERATE CHROMOSMES ##
## create recombination map for msprime and slim to simulate unlinked chromsomes
    rec_rows=[]
        # cycle through each chromosome
    for c in range(nChrom):
        print(c)
        rec_dict = {}
            # generate start position of chromsome
        rec_dict.update({"positions": c*chromSize})
            # assign recombination rate of chromosome
        rec_dict.update({"rates":recRate})
        rec_rows.append(rec_dict)
        rec_dict = {}
            # assign end position of chromosome
        rec_dict.update({"positions": int((c+1)*chromSize-1)})
            # assign recombination rate of 0.5 to create break between chromosomes
        rec_dict.update({"rates":0.5})
        rec_rows.append(rec_dict)

    if l > 10: ## add additional seasonal loci to the end so loci contribute to fitness
            rec_dict = {}
            # assign end position of site
            rec_dict.update({"positions": (nChrom*chromSize+l-10-1)})
            # assign recombination rate of 0.5 to create break between chromosomes
            rec_dict.update({"rates":0.5})
            rec_rows.append(rec_dict)

    rec_data = pd.DataFrame(rec_rows)


## generate slim recombination map
        # reformat to comply with recombinate map required in slim
    slim_rec = rec_data.positions[1:]
    slim_rec = slim_rec.reset_index(drop=True)
    slim_rec= pd.concat([slim_rec,rec_data.rates[0:-1]],axis=1)


## formulate msprime recombination map
    rec_data.positions.iloc[-1]=rec_data.positions.iloc[-1]+1
    rec_map = msprime.RateMap(position = list(rec_data.positions.astype(int)), rate= list(rec_data.rates[:-1]))

## output slim recombiation map to text file to be read into forward slim simulation
    slim_rec.to_csv("{0}{1}/rec_map_group_{1}.txt".format(tmpdir, group), index=False, header = False, sep = "\t")

    return rec_map


def simulate_burnin(tmpdir, group, l, sim_run, rec_map, s_pop, burnin_Ne, chromSize, nChrom):

    ## COALESCENT BURN IN
    start_time = time.time()
    if l > 10:
        genomeSize = chromSize*nChrom+(l-10)
    else: 
        genomeSize = chromSize*nChrom
        
    ##daiquiri.setup(level="DEBUG") ##debug
    burnin = msprime.sim_ancestry(samples=s_pop,population_size=burnin_Ne, recombination_rate=rec_map)
    ##check burnin size = genome size
    burn_length =burnin.get_sequence_length()
    if burn_length!=genomeSize:
        print("Burnin sequence length not equal to genome size!")
    burnin_ts = pyslim.annotate_defaults(burnin, model_type="WF", slim_generation=1)
    burnin_ts.dump("{0}{1}/burnin_seglift_group_{1}_{2}.trees".format(tmpdir,group,sim_run))
    print("Time for burnin = ", (time.time()- start_time))

def simulate_seglift_cap(tmpdir, slim_sim, group, sim_run, recRate, nChrom, chromSize, s_pop, w_pop, l, y, d, rGen, fitness_on, sum_gen, win_gen, offCap):
    if l > 10:
        genomeSize = chromSize*nChrom+(l-10)
    else: 
        genomeSize = chromSize*nChrom
        
    ## FORWARD SIMULATION
    # for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
    start_time = time.time()
    tmpdir_call = "tmpdir='" + str(tmpdir)+ "'"
    
    cmd = 'slim -d "' +str(tmpdir_call)+ '" -d fit='+ str(fitness_on)+" -d group=" + str(group) + " -d cap=" +str(offCap)+" -d nChrom=" + str(nChrom)+" -d g_s=" + str(sum_gen)+" -d g_w=" + str(win_gen)+" -d sim_run=" + str(sim_run) + " -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d n_s=" + str(int(s_pop)) + " -d n_w=" + str(int(w_pop)) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=0.0 -d rr=" + str(recRate) +   " -d rGen="+ str(rGen) +" ~/oliviaphd/hpc/hpc_seglift_capped.slim"

    print(cmd)
    os.system(cmd)
    # #print("Time for SLiM sim = ", (time.time()- start_time))
    print("Simulations took ", (time.time()-start_time) ,  " seconds")
    
def simulate_seglift_complex(tmpdir, slim_sim, group, sim_run, recRate, nChrom, chromSize, s_pop, w_pop, l, y, rGen, fitness_on, sum_gen, win_gen):
        genomeSize = chromSize*nChrom
            
        ## FORWARD SIMULATION
        # for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
        start_time = time.time()
        tmpdir_call = "tmpdir='" + str(tmpdir)+ "'"
        
        cmd = 'slim -d "' +str(tmpdir_call)+ '" -d fit='+ str(fitness_on)+" -d group=" + str(group) + " -d nChrom=" + str(nChrom)+" -d g_s=" + str(sum_gen)+" -d g_w=" + str(win_gen)+" -d sim_run=" + str(sim_run) + " -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d n_s=" + str(int(s_pop)) + " -d n_w=" + str(int(w_pop)) + " -d y=" + str(y) +  " -d mut=0.0 -d rr=" + str(recRate) +   " -d rGen="+ str(rGen) +" ~/oliviaphd/hpc/hpc_seglift_complex.slim"

        print(cmd)
        os.system(cmd)
        # #print("Time for SLiM sim = ", (time.time()- start_time))
        print("Simulations took ", (time.time()-start_time) ,  " seconds")

def simulate_seglift(tmpdir, slim_sim, group, sim_run, recRate, nChrom, chromSize, s_pop, w_pop, l, y, d, rGen, fitness_on, sum_gen, win_gen):

    if l > 10:
        genomeSize = chromSize*nChrom+(l-10)
    else: 
        genomeSize = chromSize*nChrom
        
    ## FORWARD SIMULATION
    # for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
    start_time = time.time()
    tmpdir_call = "tmpdir='" + str(tmpdir)+ "'"

    cmd = 'slim -d "' +str(tmpdir_call)+ '" -d fit='+ str(fitness_on)+" -d group=" + str(group) + " -d nChrom=" + str(nChrom)+" -d g_s=" + str(sum_gen)+" -d g_w=" + str(win_gen)+" -d sim_run=" + str(sim_run) + " -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d n_s=" + str(int(s_pop)) + " -d n_w=" + str(int(w_pop)) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=0.0 -d rr=" + str(recRate) +   " -d rGen="+ str(rGen) +" ~/oliviaphd/hpc/seglift_background.slim"

    print(cmd)
    os.system(cmd)
    # #print("Time for SLiM sim = ", (time.time()- start_time))
    print("Simulations took ", (time.time()-start_time) ,  " seconds")


def analyse(tmpdir, group, sim_run, mutRate, l, nChrom, nWin, sum_gen, win_gen):
    ## input  group(parameter identifier), simulation type, simulation run), mutation rate, 
    ##  number of selected loci, number of chromosomes and number of windows
    
 ## INPUT DATA
         # read in treesequence (ts) generated in SLiM
    slim_ts = pyslim.load("/{0}{1}/treeseq_{1}_{2}.trees".format(tmpdir,group,sim_run)).simplify()
        # extract the length of the simulate seqeunce from slim_ts
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
    win3 = np.linspace(0, mut_ts.sequence_length - (l-10), num=nWin+1).astype(int)
    win3=np.append(win3, mut_ts.sequence_length)
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
        # samp_gm=samp_ts.genotype_matrix()  
                
            # convert genotype matrix to haplotyoe array for haplotype statistics
        # h= allel.HaplotypeArray(samp_gm)
            # allele count for scikit.allel stats
        # samp_ac = h.count_alleles()
            # positions of mutations in samp_ts for scikit.allel windowed_statistic function
        # mut_positions = [mut.site for mut in samp_ts.mutations()]
        
            # generate haplotype statistics (H1, H12, H123, H2/H1)
        # hap_stats = allel.windowed_statistic(mut_positions,h,allel.garud_h, windows = al_win4)
        
            # tajimas D using tskit and branches of ts
        tajdb =  samp_ts.Tajimas_D(sample_sets=None, windows=win3, mode="branch")
            
            # wattersons theta using scikit.allel
        # theta_w= allel.windowed_statistic(mut_positions, (mut_positions, samp_ac), allel.watterson_theta, windows = al_win4)

            # calculate diversity (tajima's pi) using tskit
        div = samp_ts.diversity(sample_sets = None, windows = win3)  ##fix windows
        
        
            # check that all stats have ben calculated over the correct number of windows
        ts_tests = [div, tajdb]
        # al_tests = [theta_w[0],hap_stats[0]]
        for test in ts_tests:
            if len(test)!= nWin+1:
                print("error in ts test ", test, ", number of values does not match number of windows")
                
        # for test in al_tests:
        #     if len(test)!= alwin:
        #         print("error in al test ", test, ", number of values does not match number of windows")
                
        ## Collate summary statics into dataframe
            # loop over windows
        for w in range(nWin): 
           
            dict3={}
            dict3.update({"time":t})                        ## generation
            dict3.update({"n_win":w})                       ## identifier for window
            dict3.update({"win_start" : win3[w]})            ## window start position
            dict3.update({"win_end" : win3[w+1]-1})          ## window end position
            dict3.update({"chrom" : int(w/(nWin/nChrom))+1})  ## chromosome
            dict3.update({"tajimas_d_branch":tajdb[w]})     ## tajima's D (calculated with tskit)
            dict3.update({"diversity": div[w]})             ## diversity (tajimas pi; calculated with tskit)
            
            rows_list3.append(dict3)
            
        # for w in range(alwin): 
        #     try:
        #         h1 = hap_stats[0][w][0]
        #     except TypeError:
        #         h1="NaN"
        #     try:
        #         h12 = hap_stats[0][w][1]
        #     except TypeError:
        #         h12 = "NaN"
        #     try:
        #         h123 = hap_stats[0][w][2]
        #     except TypeError:
        #         h123 = "NaN"
        #     try:
        #         h2h1 = hap_stats[0][w][3]
        #     except TypeError:
        #          h2h1="NaN"
        #     dict4={}
        #     dict4.update({"time":t})                        ## generation
        #     dict4.update({"n_win":w})                       ## identifier for window
        #     dict4.update({"win_start" : win4[w]})            ## window start position
        #     dict4.update({"win_end" : win4[w+1]-1})          ## window end position
        #     dict4.update({"chrom" : int(w/(alwin/nChrom))+1})  ## chromosome
        #     dict4.update({"theta_w": theta_w[0][w]})        ## watterson's theta
        #     dict4.update({"H1": h1})         ## H1
        #     dict4.update({"H12":h12})         ## H12
        #     dict4.update({"H123": h123})       ## H123
        #     dict4.update({"H2H1": h2h1})       ## H2/H1
          
        #     rows_list4.append(dict4)
            
            # convert dictionary to datafram
    ts_stats = pd.DataFrame(rows_list3) 
    
            # write statistic df to text file
    ts_stats.to_string(buf = "/{0}{1}/sim_ts_stat_{1}_{2}.txt".format(tmpdir,group,sim_run), index=False)
    
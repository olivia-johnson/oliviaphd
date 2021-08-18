import os
import msprime
import pyslim
import numpy as np
import tskit
import pandas as pd
import time
import allel

def allele_recombination(tmpdir, params, l, nChrom, chromSize, recRate): ## only seasonal alleles simulated
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

    rec_data = pd.DataFrame(rec_rows)

## generate slim recombination map
        # reformat to comply with recombinate map required in slim
    slim_rec = rec_data.positions[1:]
    slim_rec = slim_rec.reset_index(drop=True)
    slim_rec= pd.concat([slim_rec,rec_data.rates[0:-1]],axis=1)

## output slim recombiation map to text file to be read into forward slim simulation
    slim_rec.to_csv("/{0}/rec_map_{1}.txt".format(tmpdir, params), index=False, header = False, sep = "\t")



def recombination_map(tmpdir, params, sim_type, group, l, nChrom, chromSize, recRate):
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

    if sim_type == "seglift_l10": ## add additional seasonal loci to the end so loci contribute to fitness
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
    rec_map = msprime.RecombinationMap(positions = list(rec_data.positions.astype(int)), rates= list(rec_data.rates), num_loci = int(nChrom*chromSize-1))

## output slim recombiation map to text file to be read into forward slim simulation
    slim_rec.to_csv("/{0}/rec_map_{1}.txt".format(tmpdir, params), index=False, header = False, sep = "\t")

    return rec_map


def simulate_burnin(tmpdir, params, sim_run, rec_map, s_pop, burnin_Ne):

    ## COALESCENT BURN IN
    start_time = time.time()
    burnin = msprime.simulate(sample_size=2*s_pop,Ne=burnin_Ne, mutation_rate=0, recombination_map=rec_map)
    burnin_ts = pyslim.annotate_defaults(burnin, model_type="WF", slim_generation=1)
    burnin_ts.dump("/{0}/burnin_seglift_{1}_{2}.trees".format(tmpdir,params,sim_run))
    print("Time for burnin = ", (time.time()- start_time))


def simulate_seglift(tmpdir,slim_sim, params, sim_run, recRate, nChrom, chromSize, s_pop, w_pop, l, y, d, rGen, fitness_on, sum_gen, win_gen):

    genomeSize = chromSize*nChrom
    ## FORWARD SIMULATION
    # for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
    start_time = time.time()
    cmd = "slim -d tmpdir=" +str(tmpdir)+"-d fit="+ str(fitness_on)+" -d group=" + str(params) + " -d nChrom=" + str(nChrom)+" -d g_s=" + str(sum_gen)+" -d g_w=" + str(win_gen)+" -d sim_run=" + str(sim_run) + " -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d n_s=" + str(int(s_pop)) + " -d n_w=" + str(int(w_pop)) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=0.0 -d rr=" + str(recRate) +   " -d rGen="+ str(rGen) +" ~/home/a1704225/oliviaphd/hpc/" + slim_sim + ".slim"
    print(cmd)
    os.system(cmd)
    # #print("Time for SLiM sim = ", (time.time()- start_time))
    print("Simulations took ", (time.time()-start_time) ,  " seconds")

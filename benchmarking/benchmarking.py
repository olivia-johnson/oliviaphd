#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 14:10:28 2022

@author: olivia
"""
import sys
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
import tracemalloc

tmpdir="/Users/olivia/phd_data/benchmarking/"
s_pop=10000
sequenceSize=1e6
recRate=1e-6
mutRate=1e-7
group=0

### coalescent burn in
burnin_time=[]
memory=[]
for i in range(10):
# starting the monitoring
    tracemalloc.start()
    
    # function call
    start_time = time.time()
    # msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)

    burnin = msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)
    burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
    burnin_ts.dump("{0}/burnin_seglift_group_{1}_{2}.trees".format(tmpdir,group,i))
  
    end_time=time.time()
    # displaying the memory  IN BYTES
    memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()

    burnin_time.append(end_time-start_time)


data = {'time': burnin_time, 'memory': memory}
bench_data = pd.DataFrame(data)
bench_data.to_csv("{0}/burnin_benchmarking_msprime_coalescence.txt".format(tmpdir), index=False, header = True, sep = "\t")


### coalescent burn in  WITH MUTATIONS
m_burnin_time=[]
m_memory=[]
for i in range(10):
# starting the monitoring
    tracemalloc.start()
    
    # function call
    start_time = time.time()
    # msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)

    burnin = msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)
    burnin_m=msprime.sim_mutations(burnin, rate=mutRate)
    burnin_ts = pyslim.annotate(burnin_m, model_type="WF", tick=1,    stage="late")
    burnin_ts.dump("{0}/mutburnin_seglift_group_{1}_{2}.trees".format(tmpdir,group,i))
  
    end_time=time.time()
    # displaying the memory  IN BYTES
    m_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()

    m_burnin_time.append(end_time-start_time)

mdata = {'time': m_burnin_time, 'memory': m_memory}
mbench_data = pd.DataFrame(mdata)
mbench_data.to_csv("{0}/burnin_benchmarking_msprime_mutations.txt".format(tmpdir), index=False, header = True, sep = "\t")

sm_burnin_time=[]
sm_memory=[]
for i in range(2):
# starting the monitoring
    # tracemalloc.start()
    
    # function call
    start_time = time.time()
    # msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)

    os.system("slim ~/oliviaphd/benchmarking/burnin_mut.slim")
  
    end_time=time.time()
    # displaying the memory  IN BYTES
    # sm_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    # tracemalloc.stop()

    sm_burnin_time.append(end_time-start_time)
    



sts_burnin_time=[]
sts_memory=[]
for i in range(10):
# starting the monitoring
    # tracemalloc.start()
    
    # function call
    start_time = time.time()
    # msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)

    os.system("slim ~/oliviaphd/benchmarking/burnin_ts.slim")
  
    end_time=time.time()
    # displaying the memory  IN BYTES
    # sm_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    # tracemalloc.stop()

    sts_burnin_time.append(end_time-start_time)


# smdata = {'time': sm_burnin_time, 'memory': m_memory}
# smbench_data = pd.DataFrame(smdata)
# smbench_data.to_csv("{0}/burnin_benchmarking_slim_mutations.txt".format(tmpdir), index=False, header = True, sep = "\t")



## FORWARD SIMS  #####

fts_burnin_time=[]
fts_sim_type=[]
fts_model=[]
for i in range(10):
# starting the monitoring
    # tracemalloc.start()
    
    # function call
    tracemalloc.start()
    
    # function call
    start_time = time.time()
    # msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)

    burnin = msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)
    burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
    burnin_ts.dump("{0}/burnints_seglift_group_{1}_{2}.trees".format(tmpdir,group,i))
    
    os.system("slim -d sim_run={0} ~/oliviaphd/benchmarking/forward_single_locus.slim".format(str(i)) )
  
    slim_ts=tskit.load("{0}/treeseq_single_group_{1}_{2}.trees".format(tmpdir,group,i)).simplify()
    
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)


    end_time=time.time()
    # displaying the memory  IN BYTES
    # fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()

    fts_burnin_time.append(end_time-start_time)
    fts_sim_type.append("ts")
    fts_model.append("single_locus")
    
for i in range(10):
# starting the monitoring
    # tracemalloc.start()
    
    # function call
    tracemalloc.start()
    
    # function call
    start_time = time.time()
    # msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)

    burnin = msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)
    burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
    burnin_ts.dump("{0}/burnints_seglift_group_{1}_{2}.trees".format(tmpdir,group,i))
    
    os.system("slim -d sim_run={0} ~/oliviaphd/benchmarking/forward_multi_locus.slim".format(str(i)))
  
    slim_ts=tskit.load("{0}/treeseq_multi_group_{1}_{2}.trees".format(tmpdir,group,i)).simplify()
    
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)


    end_time=time.time()
    # displaying the memory  IN BYTES
    # fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()

    fts_burnin_time.append(end_time-start_time)
    fts_sim_type.append("ts")
    fts_model.append("multilocus")

smdata = {'sim_type': fts_sim_type,  'model': fts_model,'time': fts_burnin_time}
smbench_data = pd.DataFrame(smdata)
smbench_data.to_csv("{0}/forward_benchmarking_time_multi.txt".format(tmpdir), index=False, header = True, sep = "\t")


## Analysis  #####

fts_burnin_time=[]
fts_stat_type=[]
fts_model=[]
fts_memory=[]
fts_stat=[]
for i in range(10):
# starting the monitoring
    # tracemalloc.start()
    slim_ts=tskit.load("{0}/treeseq_single_group_{1}_{2}.trees".format(tmpdir,group,i)).simplify()
 
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)

    # function call
    tracemalloc.start()
    
    # function call
    start_time = time.time()
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)
    mut_ts.diversity()
   

    end_time=time.time()
    # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()

    fts_burnin_time.append(end_time-start_time)
    fts_stat_type.append("ts")
    fts_stat.append("div")
    fts_model.append("single_locus")
    
    tracemalloc.start()
   
   # function call
    start_time = time.time()
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)
    mut_ts.Tajimas_D()
  

    end_time=time.time()
   # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
   # stopping the library
    tracemalloc.stop()

    fts_burnin_time.append(end_time-start_time)
    fts_stat_type.append("ts")
    fts_stat.append("tajimas d")
    fts_model.append("single_locus")
   
    tracemalloc.start()
    
    # function call
    start_time = time.time()
    gm=mut_ts.genotype_matrix()

        # convert genotype matrix to haplotyoe array for haplotype statistics
    h= allel.HaplotypeArray(gm)
        # allele count for scikit.allel stats
    ac = h.count_alleles()
     
    
    mut_positions = [int(var.position+1) for var in mut_ts.variants()]

    allel.sequence_diversity(mut_positions, ac)
   

    end_time=time.time()
    # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()

    fts_burnin_time.append(end_time-start_time)
    fts_stat_type.append("allele")
    fts_stat.append("div")
    fts_model.append("single_locus")
    
    tracemalloc.start()
   
   # function call
    start_time = time.time()

    gm=mut_ts.genotype_matrix()

       # convert genotype matrix to haplotyoe array for haplotype statistics
    h= allel.HaplotypeArray(gm)
       # allele count for scikit.allel stats
    ac = h.count_alleles()
    
   
    mut_positions = [int(var.position+1) for var in mut_ts.variants()]

    allel.tajima_d( ac, mut_positions)
  

  

    end_time=time.time()
   # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
   # stopping the library
    tracemalloc.stop()

    fts_burnin_time.append(end_time-start_time)
    fts_stat_type.append("allele")
    fts_stat.append("tajimas d")
    fts_model.append("single_locus")
   
    
for i in range(10):

 
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)

    # function call
    tracemalloc.start()
    
    # function call
    start_time = time.time()
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)
    mut_ts.diversity()
   

    end_time=time.time()
    # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()

    fts_burnin_time.append(end_time-start_time)
    fts_stat_type.append("ts")
    fts_stat.append("div")
    fts_model.append("multilocus")
    
    tracemalloc.start()
   
   # function call
    start_time = time.time()
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)
    mut_ts.Tajimas_D()
  

    end_time=time.time()
   # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
   # stopping the library
    tracemalloc.stop()

    fts_burnin_time.append(end_time-start_time)
    fts_stat_type.append("ts")
    fts_stat.append("tajimas d")
    fts_model.append("multilocus")
   
    tracemalloc.start()
    
    # function call
    start_time = time.time()
    gm=mut_ts.genotype_matrix()

        # convert genotype matrix to haplotyoe array for haplotype statistics
    h= allel.HaplotypeArray(gm)
        # allele count for scikit.allel stats
    ac = h.count_alleles()
     
    
    mut_positions = [int(var.position+1) for var in mut_ts.variants()]

    allel.sequence_diversity(mut_positions, ac)
   

    end_time=time.time()
    # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()

    fts_burnin_time.append(end_time-start_time)
    fts_stat_type.append("allele")
    fts_stat.append("div")
    fts_model.append("multilocus")
    
    tracemalloc.start()
   
   # function call
    start_time = time.time()

    gm=mut_ts.genotype_matrix()

       # convert genotype matrix to haplotyoe array for haplotype statistics
    h= allel.HaplotypeArray(gm)
       # allele count for scikit.allel stats
    ac = h.count_alleles()
    
   
    mut_positions = [int(var.position+1) for var in mut_ts.variants()]

    allel.tajima_d( ac, mut_positions)
  

  

    end_time=time.time()
   # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
   # stopping the library
    tracemalloc.stop()

    fts_burnin_time.append(end_time-start_time)
    fts_stat_type.append("allele")
    fts_stat.append("tajimas d")
    fts_model.append("multilocus")

smdata = {'stat_type': fts_stat_type,  'model': fts_model,'time': fts_burnin_time, 'stat': fts_stat, 'memory': fts_memory}
smbench_data = pd.DataFrame(smdata)
smbench_data.to_csv("{0}/analysis_benchmarking.txt".format(tmpdir), index=False, header = True, sep = "\t")



    
    
    
    
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
import scipy


mutRate=5e-7
recRate=1e-6
sequenceSize=1e6
winSize=10000
popSize=1000


def dupes(L): ##https://stackoverflow.com/questions/9835762/how-do-i-find-the-duplicates-in-a-list-and-create-another-list-with-them
    seen = set()
    seen2 = set()
    seen_add = seen.add
    seen2_add = seen2.add
    for item in L:
        if item in seen:
            seen2_add(item)
        else:
            seen_add(item)
    return list(seen2)

# burn-in
burnin = msprime.sim_ancestry(samples=popSize, population_size=popSize, recombination_rate=recRate, sequence_length=sequenceSize)
burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
burnin_ts.dump("/Users/olivia/oliviaphd/r2_test/burnin.trees")
## run slim
os.system('slim ~/oliviaphd/r2_test/sim.slim')
#read in ts generated in slim
slim_ts=tskit.load("/Users/olivia/oliviaphd/r2_test/ts.trees").simplify()

## overlay mutations (binary model)
mut_ts =msprime.sim_mutations(slim_ts, rate=mutRate, model = 'binary', keep=False)

print("Number of mutations:", mut_ts.num_mutations) ## number of mutations 
print("Number of sites with mutations:",mut_ts.num_sites)     ## number of sites that have mutations
all_mutations = mut_ts.mutations_site  ## sites of all mutations (includes duplicate positions)
duplicate_mutations = dupes(all_mutations)  ##identify duplicate positions
dupe_no=mut_ts.num_mutations-mut_ts.num_sites  #numebr of duplicates
print("No. duplicate mutations:", str(dupe_no))

# if dupe_no==len(duplicate_mutations):
#     print("No. duplicate mutations:", str(dupe_no))



## create windows for allele frequency based stats
win3 = np.arange(winSize/2, mut_ts.sequence_length, winSize).astype(int)
win3 = np.append(win3, mut_ts.sequence_length)
win3=np.insert(win3, 0, 0)

al_win3 = []
for w in range(len(win3)-1):
      if w == len(win3)-1:
          window = [win3[w], int(mut_ts.sequence_length)]
      else:
          window = [win3[w], win3[w+1]-1]
      al_win3.append(window)
      
## create input for r2 function
     ## create genotype array to put into scikit-allel  shape ==(var, ind)
samp_gm=mut_ts.genotype_matrix()
     # convert genotype matrix to haplotyoe array for haplotype statistics
h= allel.HaplotypeArray(samp_gm)
     # allele count for scikit.allel stats
samp_ac = h.count_alleles()
     # positions of mutations in samp_ts for scikit.allel windowed_statistic function
mut_positions = [int(var.position+1) for var in mut_ts.variants()]
 
     ## crete genotype array for LD
odds = h[:,::2]
evens = h[:,1::2]
 
gn=odds+evens
gn=gn.view('int8')  ##input for r2 a gentoype array made up of 0,1,2

np.argwhere(gn>2) ## values greater than 2

         
print("Length of genotype data:",len(gn))
print("Length of position vector:",len(mut_positions))
print("Length of unique in position vector:",len(np.unique(mut_positions)))


r2=allel.windowed_r_squared(mut_positions, gn, windows=al_win3, fill=None)
print(r2[0])

nan_start=np.take(r2[1][:, 0],np.where(np.isnan(r2[0]))) ## start of nan windows
nan_end=np.take(r2[1][:, 1],np.where(np.isnan(r2[0]))) ##end of nan windows
nans ={'start': nan_start, 'end': nan_end} ## start and end of vector

print("No. sites with nan val:", sum(np.isnan(r2[0])))
print("Max seg sites in segment:", max(r2[2]))
print("Mainseg sites in segment:", min(r2[2]))
print("Seg sites in window with nan val: ", np.take((r2[2]),np.where(np.isnan(r2[0]))) ) ##the seg sites of windows with nan values
print("Windows of nan val:", nans)  ##the start position of windows with nan values





## Compared to previous method using infinite_alleles  (This will overright variables above)
mut_ts =msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)

print("Number of mutations (ia):", mut_ts.num_mutations) ## number of mutations 
print("Number of sites with mutations (ia):",mut_ts.num_sites)     ## number of sites that have mutations
all_mutations = mut_ts.mutations_site  ## sites of all mutations (includes duplicate positions)
duplicate_mutations = dupes(all_mutations)  ##identify duplicate positions
dupe_no=mut_ts.num_mutations-mut_ts.num_sites  #numebr of duplicates
print("No. duplicate mutations (ia):", str(dupe_no))
    
samp_gm=mut_ts.genotype_matrix()
     # convert genotype matrix to haplotyoe array for haplotype statistics
h= allel.HaplotypeArray(samp_gm)
     # allele count for scikit.allel stats
samp_ac = h.count_alleles()
     # positions of mutations in samp_ts for scikit.allel windowed_statistic function
mut_positions = [int(var.position+1) for var in mut_ts.variants()]
 
     ## crete genotype array for LD
odds = h[:,::2]
evens = h[:,1::2]
 
gn=odds+evens
gn=gn.view('int8')  ##input for r2 a gentoype array made up of 0,1,2
         
print("Length of genotype data (ia):",len(gn))
print("Length of position vector (ia):",len(mut_positions))
print("Length of unique in position vector (ia):",len(np.unique(mut_positions)))


r2=allel.windowed_r_squared(mut_positions, gn, windows=al_win3)
print(r2[0])

nan_start=np.take(r2[1][:, 0],np.where(np.isnan(r2[0]))) ## start of nan windows
nan_end=np.take(r2[1][:, 1],np.where(np.isnan(r2[0]))) ##end of nan windows
nans =np.reshape( np.concatenate((nan_start, nan_end), axis=0), (len(nan_start[0]), 2))  ## start and end of vector

print("No. sites with nan val (ia):", sum(np.isnan(r2[0])))
print("Seg sites in window with nan val (ia): ", np.take((r2[2]),np.where(np.isnan(r2[0]))) ) ##the seg sites of windows with nan values
print("Windows of nan val (ia):", nans)  ##the start position of windows with nan values




## testing infinitie sites  -NOT WORKING

##*mut_pos by 1e10 makes unique muattion positions, messes up statistics
mut_ts =msprime.sim_mutations(slim_ts, rate=mutRate, discrete_genome=False, keep=False)

print("Number of mutations (ia):", mut_ts.num_mutations) ## number of mutations 
print("Number of sites with mutations (ia):",mut_ts.num_sites)     ## number of sites that have mutations
all_mutations = mut_ts.mutations_site  ## sites of all mutations (includes duplicate positions)
duplicate_mutations = dupes(all_mutations)  ##identify duplicate positions
dupe_no=mut_ts.num_mutations-mut_ts.num_sites  #numebr of duplicates
print("No. duplicate mutations (ia):", str(dupe_no))
    
al_win3 = []
for w in range(len(win3)-1):
      if w == len(win3)-1:
          window = [win3[w]*1e10, int(mut_ts.sequence_length)*1e10]
      else:
          window = [win3[w]*1e10, win3[w+1]*1e10-1]
      al_win3.append(window)

samp_gm=mut_ts.genotype_matrix()
     # convert genotype matrix to haplotyoe array for haplotype statistics
h= allel.HaplotypeArray(samp_gm)
     # allele count for scikit.allel stats
samp_ac = h.count_alleles()
     # positions of mutations in samp_ts for scikit.allel windowed_statistic function
mut_positions = [int(var.position*1e10+1) for var in mut_ts.variants()]  ## 
dupes(mut_positions) ## positions that are dupicated
 
     ## crete genotype array for LD
odds = h[:,::2]
evens = h[:,1::2]
 
gn=odds+evens
gn=gn.view('int8')  ##input for r2 a gentoype array made up of 0,1,2
         
print("Length of genotype data (ia):",len(gn))
print("Length of position vector (ia):",len(mut_positions))
print("Length of unique in position vector (ia):",len(np.unique(mut_positions)))


r2=allel.windowed_r_squared(mut_positions, gn, windows=al_win3)
print(r2[0])

tajdb =  mut_ts.Tajimas_D(sample_sets=None, windows=win3, mode="branch")
tajds =  mut_ts.Tajimas_D(sample_sets=None, windows=win3, mode="site")
tajda= allel.windowed_tajima_d(mut_positions, samp_ac, windows=al_win3)
tajds-tajda[0] ##differnce between tajimas d site from tskit and scikit allele tajimas d


nan_start=np.take(r2[1][:, 0],np.where(np.isnan(r2[0]))) ## start of nan windows
nan_end=np.take(r2[1][:, 1],np.where(np.isnan(r2[0]))) ##end of nan windows
nans =np.reshape( np.concatenate((nan_start, nan_end), axis=0), (len(nan_start[0]), 2))  ## start and end of vector

print("No. sites with nan val (ia):", sum(np.isnan(r2[0])))
print("Seg sites in window with nan val (ia): ", np.take((r2[2]),np.where(np.isnan(r2[0]))) ) ##the seg sites of windows with nan values
print("Windows of nan val (ia):", nans)  ##the start position of windows with nan values

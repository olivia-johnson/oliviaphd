import os
os.chdir("/Users/olivia/oliviaphd/")

import msprime
import pyslim
import tskit
import numpy as np
import pandas as pd



genomeSize = 1e4
popnSize = 100
mutRate = 1e-6
recRate = 1e-8
l = 100
y = 2.0
d = 0.6
nWin = 10
rGen = 100
#g_s = no. summer generations
#g_w = no. winter generations

## COALESCENT BURN IN

#burnin = msprime.simulate(sample_size=2*popnSize,Ne=popnSize, length=genomeSize, mutation_rate=mutRate, recombination_rate=recRate)
# issue is having mutations in this tree. Turn mut rate to 0. ALso, makes sense, as we only need the geneaology from the burn in, and we add all muations at the end of combined coalescent + forward time.
burnin = msprime.simulate(sample_size=2*popnSize,Ne=popnSize, length=genomeSize, mutation_rate=0, recombination_rate=recRate)
burnin_ts = pyslim.annotate_defaults(burnin, model_type="WF", slim_generation=1)
burnin_ts.dump("burnin.trees")

## FORWARD SIMULATION

cmd = "slim -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d N=" + str(int(popnSize)) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=0.0 -d rr=" + str(recRate) +" ~/oliviaphd/seglift.slim"
print(cmd)
os.system(cmd)

slim_ts = pyslim.load("./treeseq.trees").simplify()

## mutation check
if (slim_ts.num_sites != l):
    print("Less than " + str(l) + "introduced mutations")
else:
    print (str(l) + " introduced mutations")
## SUMMARISE MUTATIONS
mut_met = pd.DataFrame({"mut_num": [], "mut_pos": [], "mut_type":[]})
for mut in slim_ts.mutations():
    #print(mut)
    mut_met = mut_met.append({"mut_num" : [mut.site], "mut_pos" : [(mut.position)], "mut_type":["fluctuating"]}, ignore_index=True)

mut_met = mut_met.loc[mut_met.astype(str).drop_duplicates(subset=("mut_num")).index]
    

def allele_counts(ts, sample_sets=None):
    if sample_sets is None:
       sample_sets = [ts.samples()]
    def f(x):
       return x
    return ts.sample_count_stat(sample_sets, f, len(sample_sets),
               span_normalise=False, windows='sites',
               polarised=True, mode='site', strict=False)

alco = allele_counts(slim_ts, sample_sets = None) ## ASSUMING IN RIGHT ORDER
mut_met["allele_count"] = alco
mut_met.assign(allele_freq=lambda df: mut_met.allele_count/(2*popnSize))


## SUMMARISE INDIVIDUALS
ind_met = pd.DataFrame({"id": [], "time": [], "pop": [], "season" : []})   
for ind in slim_ts.individuals():
  #  print(ind)
    if ind.time % 10 == 0:
        ind_season = "S"
    else:
        ind_season = "W"
    ind_met = ind_met.append({"id" : [ind.id], "time" : [ind.time], "pop" : [ind.population], "season" : [ind_season]}, ignore_index=True)

    
## ADD NEUTRAL MUTATIONS
mut_ts = pyslim.SlimTreeSequence(msprime.mutate(slim_ts, rate=mutRate, keep=True))


  
## SUMMARY STATISTICS
##Tajima's D - make sample_sets list of lists when sampling over time
#win = np.linspace(0, genomeSize, num=wins+1)
#win = win.astype(int)


sum_stats = pd.DataFrame({"generation": [],  "win_start": [], "win_end" :[], "tajimas_d": [], "afs":[]})
#"allele_freq":[],
#"win_num":[],


def sum_stats(ts, wins):
    #create windows
    win = np.linspace(0, ts.sequence_length, num=wins+1)
    
    # calculate Tajima's D for windows
    tajd =  ts.Tajimas_D(sample_sets=None, windows=win, mode="site")
    
    #calculate allele frequency spectrum
    fs = ts.allele_frequency_spectrum(sample_sets=None, windows = win, span_normalise=False, polarised=True)
    
    
    #Fill in df
    #sum_stats = sum_stats.append({"generation": [gen],  "win_start" : win[0:wins], "win_end":[(win[1:(wins+1)]-1)], "tajimas_d": [tajd], "afs":[fs]}, ignore_index=True))
    
    df = pd.DataFrame({"win_start" : win[0:wins], "win_end":[(win[1:(wins+1)]-1)], "tajimas_d": [tajd], "afs":[fs]}, index = pd.MultiIndex.from_tuples([(gen, (1:wins))], names=['generation','n_win']))
    

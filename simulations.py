import os
os.chdir("/Users/olivia/oliviaphd/")

import msprime
import pyslim
import tskit
import numpy as np
import pandas as pd



genomeSize = 1000
popnSize = 1000
mutRate = 1e-6
recRate = 1e-8
l = 100
y = 2.0
d = 0.6

## COALESCENT BURN IN

#burnin = msprime.simulate(sample_size=2*popnSize,Ne=popnSize, length=genomeSize, mutation_rate=mutRate, recombination_rate=recRate)
#burnin_ts = pyslim.annotate_defaults(simp, model_type="WF", slim_generation=1. reference_sequence=None)
#burnin_ts.dump("./burnin.trees")

## error in burnin if not specifiying nucleotide model of mutation.
burnin = msprime.simulate(sample_size=2*popnSize,Ne=popnSize, length=genomeSize, mutation_rate=0, recombination_rate=recRate)
mut_burnin=msprime.mutate(burnin,rate=mutRate,model=msprime.InfiniteSites(alphabet=1)) # alphabet=1 specifies nucleotide
# check mutations have been placed
for variant in mut_burnin.variants():
	print(variant.site.id, variant.site.position, variant.alleles, variant.genotypes, sep="\t")

burnin_ts = pyslim.annotate_defaults(mut_burnin, model_type="WF", slim_generation=1, reference_sequence=None)
# now i get an error here. Is something to do with the mutation type (nucleotide)
## FORWARD SIMULATION

cmd = "slim -d GenomeSize=" + str(genomeSize) + " -d L=" + str(l)+ " -d N=" + str(popnSize) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=0.0 -d rr=" + str(recRate) +" ~/oliviaphd/seglift.slim"
print(cmd)
os.system(cmd)

slim_ts = pyslim.load("./treeseq.trees").simplify()

slim_ts.num_mutations

## ADD NEUTRAL MUTATIONS
mut_ts = pyslim.SlimTreeSequence(msprime.mutate(slim_ts, rate=mutRate, keep=True))

for mut in mut_ts.mutations():
    md = pyslim.decode_mutation(mut.metadata) 
  



## SUMMARY STATISTICS

##Tajima's D - make sample_sets list of lists when sampling over time
win = np.linspace(0, genomeSize, num=11)
win = win.astype(int)
tajd = mut_ts.Tajimas_D(sample_sets=None, windows=win, mode="site")

sum_stats = pd.DataFrame([win[0:10], (win[1:11]-1), tajd], ["win_start", "win_end", "Tajimas D"])
sum_stats.pivot(columns='var', values='val')


sum_stats = pd.DataFrame({"win_start" : win[0:10], "win_end" :(win[1:11]-1), "Tajimas D": tajd}, index = pd.MultiIndex.from_tuples([(2000,1), (2000,2), (2000,3),(2000,4),(2000,5),(2000,6),(2000,6),(2000,7),(2000,8),(2000,9),], names=["Generation", "Window"]))

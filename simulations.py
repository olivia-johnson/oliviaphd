import os
import msprime
import pyslim

genomeSize = 1000
popnSize = 1000
mutRate = 1e-6
recRate = 1e-8
l = 100
y = 2.0
d = 0.6



## Call simulation

cmd = "slim -d GenomeSize=" + str(genomeSize) + " -d L=" + str(l)+ " -d N=" + str(popnSize) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=0.0 -d rr=" + str(recRate) +" seglift.slim"
print(cmd)
os.system(cmd)

sim_ts = pyslim.load("./treeseq.trees")

## burn in

#ts = msprime.simulate(sample_size = 100, Ne=popnSize, length=genomeSize, mutation_rate=mutRate, recombination_rate=recRate)
#slim_ts = pyslim.annotate_defaults(ts, model_type="WF", slim_generation=1)
#slim_ts.dump("./burnin.trees")


burnin_ts = pyslim.SlimTreeSequence(msprime.mutate(sim_ts, rate=1e-6, keep=True))


## Format Output



## Summary Statistics

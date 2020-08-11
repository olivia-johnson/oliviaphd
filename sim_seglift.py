import os
os.chdir("/Users/olivia/oliviaphd/")

import msprime
import pyslim
import tskit
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import itertools


genomeSize = int(1e6)
popnSize = int(1e4)
mutRate = 1e-6
recRate = 1e-8
l = 20
y = 2.0
d = 0.6
nWin = 20
rGen = 100
sum_gen = 8#no. summer generations
win_gen = 3#no. winter generations

## COALESCENT BURN IN

#burnin = msprime.simulate(sample_size=2*popnSize,Ne=popnSize, length=genomeSize, mutation_rate=mutRate, recombination_rate=recRate)
# issue is having mutations in this tree. Turn mut rate to 0. ALso, makes sense, as we only need the geneaology from the burn in, and we add all muations at the end of combined coalescent + forward time.
burnin = msprime.simulate(sample_size=popnSize,Ne=int(1e6), length=genomeSize, mutation_rate=0, recombination_rate=recRate)
burnin_ts = pyslim.annotate_defaults(burnin, model_type="WF", slim_generation=1)
burnin_ts.dump("burnin_seglift.trees")

## FORWARD SIMULATION
# for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
cmd = "slim -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d N=" + str(int(popnSize)) + " -d y=" + str(y) + " -d d=" + str(d) + " -d mut=0.0 -d rr=" + str(recRate) +" ~/oliviaphd/seglift.slim"
print(cmd)
os.system(cmd)

slim_ts = pyslim.load("./treeseq_seglift.trees").simplify()

## mutation check
if (slim_ts.num_sites != l):
    print("Less than " + str(l) + "introduced mutations")
else:
    print (str(l) + " introduced mutations")
    
## SUMMARISE MUTATIONS
start_time = time.time()
mut_met = pd.DataFrame({"mut_num": [], "mut_pos": []})
for mut in slim_ts.mutations():
    #print(mut)
    mut_met = mut_met.append({"mut_num" : mut.site, "mut_pos" : mut.position}, ignore_index=True)

mut_met = mut_met.loc[mut_met.astype(str).drop_duplicates(subset=("mut_num")).index]
    
print("Time for initial mut info = ", (time.time()- start_time))


def allele_counts(ts, sample_sets=None):
    if sample_sets is None:
       sample_sets = [ts.samples()]
    def f(x):
       return x
    return ts.sample_count_stat(sample_sets, f, len(sample_sets),
               span_normalise=False, windows='sites',
               polarised=True, mode='site', strict=False)


## SUMMARISE INDIVIDUALS
start_time = time.time()
rows_list = []

for ind in slim_ts.individuals():
 
  #  print(ind)
    if ind.time % 10 == 0:
        ind_season = "S"
    else:
        ind_season = "W"

    dict1 = {}
        # get input row in dictionary format
        # key = col_name
    dict1.update({"id" : ind.id})
    dict1.update({"time" : ind.time})
    dict1.update({"pop" : ind.population})
    dict1.update({"season" : ind_season})
    dict1.update({"nodes": ind.nodes}) 

    rows_list.append(dict1)

ind_met = pd.DataFrame(rows_list) 
print("Time for ind table = ", (time.time()- start_time))
## CALCULATE FITNESS
# samps = ind_met[["id", "nodes"]]


# start_time = time.time()
# counts = allele_counts(slim_ts, sample_sets = list(samps.nodes))
# print("Time for counts = " + (start_time - time.time()))
# start_time = time.time()
# if  len(samps) == len(counts[1]):
#     for i in range(len(samps)):
#         s_count = counts[:,i]
#         ns = sum(s_count == 2)
#         nhet = sum(s_count == 1)
#         nw = sum(s_count == 0)
#         z_s = ns + (d*nhet)
#         z_w = nw + (d*nhet)
#         s_fit = (1+z_s)**y
#         w_fit = (1+z_w)**y
#         if str(samps.nodes[i]) == str(ind_met.nodes[i]):
#             ind_met.loc[ind_met.id == samps.id[i], ["s_fitness", "w_fitness"]] = [s_fit, w_fit]
#         else:
#             print("Nodes do not match for sample " + str(i))
# else:
#     print("samples do not match allele counts")

# print("Time for fitness calc = ", (time.time()- start_time))

#for i in list(np.unique(ind_met.time)):

    
## ADD NEUTRAL MUTATIONS
mut_ts = pyslim.SlimTreeSequence(msprime.mutate(slim_ts, rate=mutRate, keep=True))

start_time = time.time()


rows_list2 = []
a_list = list(mut_met.mut_pos)
for mut in mut_ts.mutations():
    #print(mut)

        #print(mut)
        
    given_value = mut.position
        
    absolute_difference_function = lambda list_value : abs(list_value - given_value)

    closest_value = min(a_list, key=absolute_difference_function)
        
    dist = given_value - closest_value
    
    num = mut_met.mut_num[mut_met.mut_pos == closest_value]
    if mut.metadata == []:
        m_type = "neutral"
    else:
        m_type = "fluctuating"
        
    dict2 = {}
        # get input row in dictionary format
        # key = col_name
    dict2.update({"mut_id":mut.id})
    dict2.update({"mut_num" : mut.site})
    dict2.update({"mut_pos" : mut.position})
    dict2.update({"mut_type":m_type})
    dict2.update({"nearest_dist": abs(dist)}) 
    dict2.update({"nearest_mut_pos": closest_value})
    dict2.update({"nearest_mut_num": int(num)})

    rows_list2.append(dict2)

n_met = pd.DataFrame(rows_list2) 


print("Time for mut table = ", (time.time()- start_time))

   # n_met = n_met.append({"mut_id":mut.id, "mut_num" : mut.site, "mut_pos" : mut.position, "mut_type":type, "nearest_dist": abs(dist),"nearest_mut_pos": closest_value, "nearest_mut_num": int(num)}, ignore_index=True)

o_len = len(n_met)
s_len = len(n_met[n_met.mut_type == 'fluctuating'])

mut_ud = n_met.loc[n_met.astype(str).drop_duplicates(subset=("mut_num")).index]
n_len = len(mut_ud)

if n_len == (o_len - s_len + l):
    print("Only duplicate selected mutations removed")
else:
    print("More than duplicate selected mutations removed")
# nalco = allele_counts(mut_ts, sample_sets = None) ## ASSUMING IN RIGHT ORDER
# n_met["allele_count"] = nalco
# n_met = n_met.assign(allele_freq=lambda df: ((n_met.allele_count/(2*popnSize))))

## PLot dist by freq
   


   
## REMOVE SELECTED SITES FROM TS

list_m = mut_ud.mut_pos[mut_ud.mut_type=="fluctuating"].astype(int)
list_mu = list_m .tolist()
stat_ts = pyslim.SlimTreeSequence(mut_ts.delete_sites(list_mu, record_provenance=True))

## SUM MARY STATISTICS FUNCTION
def sum_stats(ts, wins, sample_sets=None):
    if sample_sets is None:
       sample_sets = [ts.samples()]
    rows_list3 = []
    
    #create windows
    win = np.linspace(0, ts.sequence_length, num=wins+1)
    
    # calculate Tajima's D for windows
    tajds =  ts.Tajimas_D(sample_sets=sample_sets, windows=win, mode="site")
    tajdb =  ts.Tajimas_D(sample_sets=sample_sets, windows=win, mode="branch")


    div = ts.diversity(sample_sets = sample_sets, windows = win)
    for w in range(wins):
        s_mutcount = np.sum((mut_met.mut_pos >= win[w]) & (mut_met.mut_pos < (win[w+1]-1)))
        n_mutcount = np.sum((mut_ud.mut_pos[mut_ud.mut_type=="neutral"] >= win[w]) & (mut_ud.mut_pos[mut_ud.mut_type=="neutral"] < (win[w+1]-1)))
       # print(w)
        dict3={}
        dict3.update({"n_win":w})
        dict3.update({"win_start" : win[w]})
        dict3.update({"win_end" : win[w+1]-1})
        dict3.update({"n_s_mut" : s_mutcount})
        dict3.update({"n_n_mut" : n_mutcount})
        dict3.update({"tajimas_d_site":tajds[w]})
        dict3.update({"tajimas_d_branch":tajdb[w]})
        dict3.update({"diversity": div[w]})
        #dict3.update({"afs": fs[w]}) 
       
        #print(dict3)
        rows_list3.append(dict3)
       # print(rows_list3)
### REMOVE BRACKETS
    rows_list3.append(dict3)
    s_stats = pd.DataFrame(rows_list3) 
    s_stats = s_stats.loc[s_stats.astype(str).drop_duplicates().index]

    return(s_stats)
  
    
# ## CALCULATE ALLELE COUNTS AND FREQS THROUGH TIME
# samp_list = []
# for t in np.unique(ind_met.time):
#     sample = ind_met.nodes[ind_met.time == t]
#     samples= list(itertools.chain(*sample))
#     samp_ts = stat_ts.subset(samples)
#     # samp_list.append(samples)
# nalco = allele_counts(mut_ts, sample_sets = samp_list) 
# nfreq = nalco/popnSize

# rows_list4 = []
# for i in len(nfreq[1]):
#     g_time = ind_met.time[i]
#     g_var = list(nfreq[:,i])
#     dict4 = {}
#     dict4.update({g_time.astype(str):g_var})
    
        

## CALCULATE SUMMARY STATISTICS
start_time = time.time()

for t in np.unique(ind_met.time):
    sample = ind_met.nodes[ind_met.time == t]
    samples= list(itertools.chain(*sample))
    #ac_label = "".join([str(int(t)), "_ac"])
    #af_label = "".join([str(int(t)), "_af"])
    # nalco = allele_counts(mut_ts, sample_sets = sample) ## ASSUMING IN RIGHT ORDER
    # mut_ud[allele_count] = nalco
    # mut_ud  = mut_ud.assign(allele_freq=lambda df:((mut_ud.allele_count/(2*popnSize))))

    df = sum_stats(ts = stat_ts, wins = nWin, sample_sets = [samples])
    df['Gen']=t
    if t == ind_met.time[0]:
        s_stat = df
    else:
        s_stat = pd.concat([s_stat, df])
    
# print("Time for sum stats = ", (time.time()- start_time))

# s_stat.to_string(buf = "~/oliviaphd/data/sim_s_stat.txt")

# plt.scatter(s_stat.n_win[s_stat.Gen==9800.0], s_stat.tajimas_d_branch[s_stat.Gen==9800.0], s=2)


# plt.scatter(n_met.nearest_dist[n_met.mut_type == "neutral"], n_met.allele_freq[n_met.mut_type == "neutral"], s=2)
# plt.show()
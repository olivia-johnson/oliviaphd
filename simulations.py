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


genomeSize = 1e5
popnSize = 1000
mutRate = 1e-5
recRate = 1e-8
l = 100
y = 2.0
d = 0.6
nWin = 20
rGen = 100
sum_gen = 8#no. summer generations
win_gen = 3#no. winter generations

## COALESCENT BURN IN

#burnin = msprime.simulate(sample_size=2*popnSize,Ne=popnSize, length=genomeSize, mutation_rate=mutRate, recombination_rate=recRate)
# issue is having mutations in this tree. Turn mut rate to 0. ALso, makes sense, as we only need the geneaology from the burn in, and we add all muations at the end of combined coalescent + forward time.
burnin = msprime.simulate(sample_size=2*popnSize,Ne=popnSize, length=genomeSize, mutation_rate=0, recombination_rate=recRate)
burnin_ts = pyslim.annotate_defaults(burnin, model_type="WF", slim_generation=1)
burnin_ts.dump("burnin.trees")

## FORWARD SIMULATION
# for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
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
mut_met = pd.DataFrame({"mut_num": [], "mut_pos": []})
for mut in slim_ts.mutations():
    #print(mut)
    mut_met = mut_met.append({"mut_num" : mut.site, "mut_pos" : mut.position}, ignore_index=True)

mut_met = mut_met.loc[mut_met.astype(str).drop_duplicates(subset=("mut_num")).index]
    



def allele_counts(ts, sample_sets=None):
    if sample_sets is None:
       sample_sets = [ts.samples()]
    def f(x):
       return x
    return ts.sample_count_stat(sample_sets, f, len(sample_sets),
               span_normalise=False, windows='sites',
               polarised=True, mode='site', strict=False)


## SUMMARISE INDIVIDUALS
rows_list = []
#ind_met = pd.DataFrame({"id": [], "time": [], "pop": [], "season" : [], "nodes":[], "s_fitness":[], "w_fitness":[]})   
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

    #ind_met = ind_met.append({"id" : ind.id, "time" : ind.time, "pop" : ind.population, "season" : ind_season, "nodes": ind.nodes}, ignore_index=True)

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

#n_met = pd.DataFrame({"mut_id":[],"mut_num": [], "mut_pos": [], "mut_type":[], "nearest_dist":[], "nearest_mut_pos":[], "nearest_mut_num":[]})

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
   


   
## SUMMARY STATISTICS
##Tajima's D - make sample_sets list of lists when sampling over time
#win = np.linspace(0, genomeSize, num=wins+1)
#win = win.astype(int)


#sum_stats = pd.DataFrame({"generation": [],  "win_start": [], "win_end" :[], "tajimas_d": [], "afs":[], "diversity":[]})
#"allele_freq":[],
#"win_num":[],



start_time = time.time()
def sum_stats(ts, wins, sample_sets=None):
    if sample_sets is None:
       sample_sets = [ts.samples()]
    #create windows
    rows_list3 = []
   
    #df = pd.DataFrame({"n_win":[],"win_start":[], "win_end":[], "tajimas_d": [], "afs":[], "diversity":[]})

    win = np.linspace(0, ts.sequence_length, num=wins+1)
    
    # calculate Tajima's D for windows
    tajd =  ts.Tajimas_D(sample_sets=sample_sets, windows=win, mode="site")
    
    #calculate allele frequency spectrum
    fs = ts.allele_frequency_spectrum(sample_sets=sample_sets, windows = win, span_normalise=False, polarised=True)
    
    div = ts.diversity(sample_sets = sample_sets, windows = win)
    for w in range(wins):
       # print(w)
        dict3={}
        dict3.update({"n_win":w})
        dict3.update({"win_start" : win[w]})
        dict3.update({"win_end" : win[w+1]-1})
        dict3.update({"tajimas_d":tajd[w]})
        #dict3.update({"afs": fs[w]}) 
        dict3.update({"diversity": div[w]})
        #print(dict3)
        rows_list3.append(dict3)
       # print(rows_list3)

    rows_list3.append(dict3)
    s_stats = pd.DataFrame(rows_list3) 
    s_stats = s_stats.loc[s_stats.astype(str).drop_duplicates().index]

    return(s_stats)
    ## SAMPLE BY WINDOW???  allele_counts(ts,sample_sets=None)
    #Fill in df
    #sum_stats = sum_stats.append({"generation": [gen],  "win_start" : win[0:wins], "win_end":[(win[1:(wins+1)]-1)], "tajimas_d": [tajd], "afs":[fs]}, ignore_index=True))
    
    #df = df.append({"n_win":(range(wins)[0:wins]),"win_start": win[0:wins], "win_end":(win[1:(wins+1)]-1), "tajimas_d": tajd, "afs":fs, "diversity":div}, ignore_index=True)
    
    
#pd.DataFrame({"win_start" : win[0:10], "win_end" :(win[1:11]-1), "Tajimas D": tajd}, index = pd.MultiIndex.from_tuples([(2000,1), (2000,2), (2000,3),(2000,4),(2000,5),(2000,6),(2000,6),(2000,7),(2000,8),(2000,9),], names=["Generation", "Window"]))

for t in np.unique(ind_met.time):
    sample = ind_met.nodes[ind_met.time == t]
    samples= list(itertools.chain(*sample))
    ac_label = "".join([str(int(t)), "_ac"])
    af_label = "".join([str(int(t)), "_af"])
    nalco = allele_counts(mut_ts, sample_sets = sample) ## ASSUMING IN RIGHT ORDER
    mut_ud[ac_label] = nalco
    mut_ud  = mut_ud.assign(af_label=lambda df: ((n_met.allele_count/(2*popnSize))))

    df = sum_stats(ts = mut_ts, wins = nWin, sample_sets = [samples])
    df['Gen']=t
    if t == ind_met.time[0]:
        s_stat = df
    else:
        s_stat = pd.concat([s_stat, df])
    
print("Time for sum stats = ", (time.time()- start_time))

samp_list = []
for t in np.unique(ind_met.time):
    sample = ind_met.nodes[ind_met.time == t]
    samples= list(itertools.chain(*sample))
    samp_list.append(samples)
nalco = allele_counts(mut_ts, sample_sets = samp_list) ## ASSUMING IN RIGHT ORDER
   

    mut_ud[ac_label] = nalco
    mut_ud  = mut_ud.assign(af_label=lambda df: ((n_met.allele_count/(2*popnSize))))

plt.scatter(n_met.nearest_dist[n_met.mut_type == "neutral"], n_met.allele_freq[n_met.mut_type == "neutral"], s=2)
plt.show()
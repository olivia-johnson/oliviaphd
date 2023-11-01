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
diversity=[]
simtype=[]
for i in range(10):
    
# starting the monitoring
    tracemalloc.start()
    
    # function call
    start_time = time.time()
    # msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)

    burnin_ts = msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize,)

  
    end_time=time.time()
    # displaying the memory  IN BYTES
    memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()
    burnin_time.append(end_time-start_time)
    burnin_ts = pyslim.annotate(burnin_ts, model_type="WF", tick=1,    stage="late")
    burnin_ts.dump("{0}/burnin_seglift_group_{1}_{2}.trees".format(tmpdir,group,i))
    tracemalloc.start()
    simtype.append("msprime_ts")
    muts_time = time.time()

    burnin_ts=msprime.sim_mutations(burnin_ts, rate=mutRate)
    mute_time=time.time()
# displaying the memory  IN BYTES
    memory.append(tracemalloc.get_traced_memory()[1])

# stopping the library
    tracemalloc.stop()
    burnin_time.append((end_time-start_time)+(mute_time-muts_time))
    div=burnin_ts.diversity(mode="site")
    simtype.append("msprime_mut")
    
    diversity.append(div)
    diversity.append(div)


data = {'time': burnin_time, 'memory': memory, 'diversity': diversity, 'sim_type':simtype}
bench_data = pd.DataFrame(data)
#bench_data.to_csv("{0}/burnin_benchmarking_msprime_coalescence.txt".format(tmpdir), index=False, header = True, sep = "\t")
bench_data.to_csv("{0}/burnin_benchmarking_msprime2.txt".format(tmpdir), index=False, header = True, sep = "\t")


sm_burnin_time=[]
sm_memory=[]
for i in range(10):
# starting the monitoring
    # tracemalloc.start()
    
    # function call
    start_time = time.time()
    # msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)

    os.system("slim -d sim_run={0} ~/oliviaphd/benchmarking/burnin_mut.slim ".format(i))
  
    end_time=time.time()
    # displaying the memory  IN BYTES
    # sm_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    # tracemalloc.stop()

    sm_burnin_time.append(end_time-start_time)


sts_diversity=[]
for i in range(10):
# starting the monitoring
    # tracemalloc.start()
    
    # function call
    # start_time = time.time()
    # msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)

    #os.system("slim -d sim_run={0} ~/oliviaphd/benchmarking/burnin_ts.slim".format(i))
  
    # end_time=time.time()
    
    burnin=tskit.load("{0}burnin_ts_{1}.trees".format(tmpdir,i)).simplify()
    burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
    burnin_ts=msprime.sim_mutations(burnin_ts, rate=mutRate)

    div=burnin_ts.diversity(mode="site")
    sts_diversity.append(div)
sts=pd.read_csv('{0}burnin_benchmarking_slim_ts.txt'.format(tmpdir), sep="\t", header=0)
sts.insert(3, "diversity", sts_diversity)
sts.to_csv("{0}burnin_benchmarking_slim_ts.txt".format(tmpdir), index=False, header = True, sep = "\t")

stscc_diversity=[]
for i in range(10):
# starting the monitoring
    # tracemalloc.start()
    
    # function call
    # start_time = time.time()
    # msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)

    os.system("slim -d sim_run={0} ~/oliviaphd/benchmarking/burnin_tscc.slim".format(i))
  
    # end_time=time.time()
    
    burnin=tskit.load("{0}burnin_tscc_{1}.trees".format(tmpdir,i)).simplify()
    burnin_ts = pyslim.annotate(burnin, model_type="WF", tick=1,    stage="late")
    burnin_ts=msprime.sim_mutations(burnin_ts, rate=mutRate)

    div=burnin_ts.diversity(mode="site")
    stscc_diversity.append(div)
sts=pd.read_csv('{0}burnin_benchmarking_slim_tscc.txt'.format(tmpdir), sep="\t", header=0)
sts.insert(3, "diversity", stscc_diversity)
sts.to_csv("{0}burnin_benchmarking_slim_tscc.txt".format(tmpdir), index=False, header = True, sep = "\t")

# smdata = {'time': sm_burnin_time, 'memory': m_memory}
# smbench_data = pd.DataFrame(smdata)
# smbench_data.to_csv("{0}/burnin_benchmarking_slim_mutations.txt".format(tmpdir), index=False, header = True, sep = "\t")



## FORWARD SIMS  #####

fts_time=[]
fts_sim_type=[]
fts_model=[]
for i in range(10):
# starting the monitoring
    # tracemalloc.start()
    
    # function call
    tracemalloc.start()
    
    # function call
    start_time = time.time()
    
    os.system("slim -d sim_run={0} ~/oliviaphd/benchmarking/forward_single_locus.slim".format(str(i)) )
  
    slim_ts=tskit.load("{0}/treeseq_single_group_{1}_{2}.trees".format(tmpdir,group,i)).simplify()
    
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)


    end_time=time.time()
    # displaying the memory  IN BYTES
    # fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()

    fts_time.append(end_time-start_time)
    fts_sim_type.append("ts")
    fts_model.append("single_locus")
    
for i in range(10):
# starting the monitoring
    # tracemalloc.start()
    
    # function call
    tracemalloc.start()
    
    # function call
    start_time = time.time()
      
    os.system("slim -d sim_run={0} ~/oliviaphd/benchmarking/forward_multi_locus.slim".format(str(i)))
  
    slim_ts=tskit.load("{0}/treeseq_multi_group_{1}_{2}.trees".format(tmpdir,group,i)).simplify()
    
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)


    end_time=time.time()
    # displaying the memory  IN BYTES
    # fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()

    fts_time.append(end_time-start_time)
    fts_sim_type.append("ts")
    fts_model.append("multilocus")

smdata = {'sim_type': fts_sim_type,  'model': fts_model,'time': fts_time}
smbench_data = pd.DataFrame(smdata)
smbench_data.to_csv("{0}/forward_benchmarking_time_ts.txt".format(tmpdir), index=False, header = True, sep = "\t")

fm_time=[]
fm_sim_type=[]
fm_model=[]
for i in range(10):
# starting the monitoring
    tracemalloc.start()
    
    # function call
    start_time = time.time()
    # msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)
    # model=msprime.SLiMMutationModel(type=0)
    # burnin = msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)
    # burnin_m=msprime.sim_mutations(burnin, rate=mutRate, model=model)
    # burnin_ts = pyslim.annotate(burnin_m, model_type="WF", tick=1,    stage="late")
    # burnin_ts.dump("{0}/mutburnin_seglift_group_{1}_{2}.trees".format(tmpdir,group,i))
  
    os.system("slim -d sim_run={0} ~/oliviaphd/benchmarking/forward_multilocus_mut.slim".format(str(i)) )

    end_time=time.time()
    # displaying the memory  IN BYTES
    # fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()

    fm_time.append(end_time-start_time)
    fm_sim_type.append("mut")
    fm_model.append("multilocus")


smdata = {'sim_type': fm_sim_type,  'model': fm_model,'time': fm_time}
smbench_data = pd.DataFrame(smdata)
smbench_data.to_csv("{0}/forward_mut_benchmarking_time_multi.txt".format(tmpdir), index=False, header = True, sep = "\t")

fm_time=[]
fm_sim_type=[]
fm_model=[]
for i in range(10):
# starting the monitoring
    tracemalloc.start()
    
    # function call
    start_time = time.time()
    # msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)
    # model=msprime.SLiMMutationModel(type=0)
    # burnin = msprime.sim_ancestry(samples=s_pop, population_size=s_pop, recombination_rate=recRate, sequence_length=sequenceSize)
    # burnin_m=msprime.sim_mutations(burnin, rate=mutRate, model=model)
    # burnin_ts = pyslim.annotate(burnin_m, model_type="WF", tick=1,    stage="late")
    # burnin_ts.dump("{0}/mutburnin_seglift_group_{1}_{2}.trees".format(tmpdir,group,i))
  
    os.system("slim -d sim_run={0} ~/oliviaphd/benchmarking/forward_single_locus_mut.slim".format(str(i)) )

    end_time=time.time()
    # displaying the memory  IN BYTES
    # fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()

    fm_time.append(end_time-start_time)
    fm_sim_type.append("mut")
    fm_model.append("single_locus")


smdata = {'sim_type': fm_sim_type,  'model': fm_model,'time': fm_time}
smbench_data = pd.DataFrame(smdata)
smbench_data.to_csv("{0}/forward_mut_benchmarking_time_single.txt".format(tmpdir), index=False, header = True, sep = "\t")


## Analysis  #####

calc_time=[]
total_time=[]
fts_stat_type=[]
fts_model=[]
fts_memory=[]
fts_stat=[]
dt=[]
for i in range(10):
# starting the monitoring
    # tracemalloc.start()
    slim_ts=tskit.load("{0}/treeseq_single_group_{1}_{2}.trees".format(tmpdir,group,i)).simplify()
 
    # function call
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)
    tot_time=time.time()

    tracemalloc.start()
    start_time = time.time()
    mut_ts.diversity(sample_sets=pyslim.individuals_alive_at(mut_ts, 0))
   

    end_time=time.time()
    # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()
    total_time.append(end_time-tot_time)
    calc_time.append(end_time-start_time)
    fts_stat_type.append("Tree-based")
    fts_stat.append("div")
    fts_model.append("single_locus")
    dt.append("Tree Sequence")
    
   
   # function call
    tracemalloc.start()
    tot_time=time.time()

    start_time = time.time()
    mut_ts.Tajimas_D(sample_sets=pyslim.individuals_alive_at(mut_ts, 0))
  

    end_time=time.time()
   # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
   # stopping the library
    tracemalloc.stop()
    total_time.append(end_time-tot_time)
    calc_time.append(end_time-start_time)
    fts_stat_type.append("Tree-based")
    fts_stat.append("tajimas d")
    fts_model.append("single_locus")
    dt.append("Tree Sequence")

    
    # function call
    tot_time = time.time()

    samp_ts = mut_ts.simplify(samples=pyslim.individuals_alive_at(mut_ts, 0))
    gm=samp_ts.genotype_matrix()

        # convert genotype matrix to haplotyoe array for haplotype statistics
    h= allel.HaplotypeArray(gm)
        # allele count for scikit.allel stats
    ac = h.count_alleles()
     
    
    mut_positions = [int(var.position+1) for var in samp_ts.variants()]
    tracemalloc.start()
    start_time = time.time()
    allel.sequence_diversity(mut_positions, ac)
   

    end_time=time.time()
    # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()
    total_time.append(end_time-tot_time)
    calc_time.append(end_time-start_time)
    fts_stat_type.append("Allele-based")
    fts_stat.append("div")
    fts_model.append("single_locus")
    dt.append("Tree Sequence")

   
   # function call
    tot_time = time.time()
    samp_ts = mut_ts.simplify(samples=pyslim.individuals_alive_at(mut_ts, 0))
    gm=samp_ts.genotype_matrix()


       # convert genotype matrix to haplotyoe array for haplotype statistics
    h= allel.HaplotypeArray(gm)
       # allele count for scikit.allel stats
    ac = h.count_alleles()
    
   
    mut_positions = [int(var.position+1) for var in samp_ts.variants()]
    tracemalloc.start()
    start_time = time.time()
    allel.tajima_d( ac, mut_positions)
  
    end_time=time.time()
   # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
   # stopping the library
    tracemalloc.stop()

    calc_time.append(end_time-start_time)
    total_time.append(end_time-tot_time)

    fts_stat_type.append("Allele-based")
    fts_stat.append("tajimas d")
    fts_model.append("single_locus")
    dt.append("Tree Sequence")


    tot_time = time.time()
    cmd="sed -n -e '/#OUT: 140000 140000 A/, /Individuals:/ p' {0}forward_single_mut{1}.txt > {0}final_gen.txt".format(tmpdir,i)
    os.system(cmd)
    sts=pd.read_csv("{0}final_gen.txt".format(tmpdir), sep=" ", header=1, skiprows=4, skipfooter=1, names=["No", "ID", "mutType", "Pos", "s", "d", "pop", "Tick", "Count"])
    sts=sts.sort_values('Pos')
    sts['Alt']=20000-sts['Count']
    positions=list(sts["Pos"])
    count=sts[["Count", "Alt"]].to_numpy()
    tracemalloc.start()
    start_time = time.time()
    allel.sequence_diversity(positions, count)
    end_time=time.time()
# displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])

# stopping the library
    tracemalloc.stop()

    calc_time.append(end_time-start_time)
    total_time.append(end_time-tot_time)

    fts_stat_type.append("Allele-based")
    fts_stat.append("div")
    fts_model.append("single_locus")
    dt.append("Classical")
   
    tot_time = time.time()
    cmd="sed -n -e '/#OUT: 140000 140000 A/, /Individuals:/ p' {0}forward_single_mut{1}.txt > {0}final_gen.txt".format(tmpdir,i)
    os.system(cmd)
    sts=pd.read_csv("{0}final_gen.txt".format(tmpdir), sep=" ", header=1, skiprows=4, skipfooter=1, names=["No", "ID", "mutType", "Pos", "s", "d", "pop", "Tick", "Count"])
    sts=sts.sort_values('Pos')
    sts['Alt']=20000-sts['Count']
    positions=list(sts["Pos"])
    count=sts[["Count", "Alt"]].to_numpy()
    tracemalloc.start()
    start_time = time.time()
    allel.tajima_d( count, positions)
    end_time=time.time()
   # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
   # stopping the library
    tracemalloc.stop()

    calc_time.append(end_time-start_time)
    total_time.append(end_time-tot_time)

    fts_stat_type.append("Allele-based")
    fts_stat.append("tajimas d")
    fts_model.append("single_locus")
    dt.append("Classical")
    
    
  
for i in range(10):

 
    slim_ts=tskit.load("{0}/treeseq_multi_group_{1}_{2}.trees".format(tmpdir,group,i)).simplify()

    
    # function call
    mut_ts = msprime.sim_mutations(slim_ts, rate=mutRate, model = 'infinite_alleles', keep=False)
    tracemalloc.start()
    tot_time=time.time()
    start_time = time.time()
    mut_ts.diversity(sample_sets=pyslim.individuals_alive_at(mut_ts, 0))
   

    end_time=time.time()
    # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()
    total_time.append(end_time-tot_time)
    calc_time.append(end_time-start_time)
    fts_stat_type.append("Tree-based")
    fts_stat.append("div")
    fts_model.append("multilocus")
    dt.append("Tree Sequence")

   
   # function call
    tracemalloc.start()
    tot_time=time.time()
    start_time = time.time()
    mut_ts.Tajimas_D(sample_sets=pyslim.individuals_alive_at(mut_ts, 0))
  

    end_time=time.time()
   # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
   # stopping the library
    tracemalloc.stop()
    total_time.append(end_time-tot_time)
    calc_time.append(end_time-start_time)
    fts_stat_type.append("Tree-based")
    fts_stat.append("tajimas d")
    fts_model.append("multilocus")
    dt.append("Tree Sequence")

    
    # function call
    tot_time = time.time()
    samp_ts = mut_ts.simplify(samples=pyslim.individuals_alive_at(mut_ts, 0))
    gm=samp_ts.genotype_matrix()


        # convert genotype matrix to haplotyoe array for haplotype statistics
    h= allel.HaplotypeArray(gm)
        # allele count for scikit.allel stats
    ac = h.count_alleles()
     
    
    mut_positions = [int(var.position+1) for var in samp_ts.variants()]
    tracemalloc.start()
    start_time = time.time()
    allel.sequence_diversity(mut_positions, ac)
   

    end_time=time.time()
    # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
    
    # stopping the library
    tracemalloc.stop()
    total_time.append(end_time-tot_time)
    calc_time.append(end_time-start_time)
    fts_stat_type.append("Allele-based")
    fts_stat.append("div")
    fts_model.append("multilocus")
    dt.append("Tree Sequence")

   
   # function call
    tot_time = time.time()

    samp_ts = mut_ts.simplify(samples=pyslim.individuals_alive_at(mut_ts, 0))
    gm=samp_ts.genotype_matrix()


       # convert genotype matrix to haplotyoe array for haplotype statistics
    h= allel.HaplotypeArray(gm)
       # allele count for scikit.allel stats
    ac = h.count_alleles()
    
   
    mut_positions = [int(var.position+1) for var in samp_ts.variants()]
    tracemalloc.start()
    start_time = time.time()
    allel.tajima_d( ac, mut_positions)
  

  

    end_time=time.time()
   # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
   # stopping the library
    tracemalloc.stop()

    calc_time.append(end_time-start_time)
    total_time.append(end_time-tot_time)

    fts_stat_type.append("Allele-based")
    fts_stat.append("tajimas d")
    fts_model.append("multilocus")
    dt.append("Tree Sequence")

    tot_time = time.time()
    cmd="sed -n -e '/#OUT: 140000 140000 A/, /Individuals:/ p' {0}forward_multi_mut{1}.txt > {0}final_gen.txt".format(tmpdir,i)
    os.system(cmd)
    sts=pd.read_csv("{0}final_gen.txt".format(tmpdir), sep=" ", header=1, skiprows=4, skipfooter=1, names=["No", "ID", "mutType", "Pos", "s", "d", "pop", "Tick", "Count"])
    sts=sts.sort_values('Pos')
    sts['Alt']=20000-sts['Count']
    positions=list(sts["Pos"])
    count=sts[["Count", "Alt"]].to_numpy()
    tracemalloc.start()
    start_time = time.time()
    allel.sequence_diversity(positions, count)
    end_time=time.time()
# displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])

# stopping the library
    tracemalloc.stop()

    calc_time.append(end_time-start_time)
    total_time.append(end_time-tot_time)

    fts_stat_type.append("Allele-based")
    fts_stat.append("div")
    fts_model.append("multilocus")
    dt.append("Classical")
   
    tot_time = time.time()
    cmd="sed -n -e '/#OUT: 140000 140000 A/, /Individuals:/ p' {0}forward_multi_mut{1}.txt > {0}final_gen.txt".format(tmpdir,i)
    os.system(cmd)
    sts=pd.read_csv("{0}final_gen.txt".format(tmpdir), sep=" ", header=1, skiprows=4, skipfooter=1, names=["No", "ID", "mutType", "Pos", "s", "d", "pop", "Tick", "Count"])
    sts=sts.sort_values('Pos')
    sts['Alt']=20000-sts['Count']
    positions=list(sts["Pos"])
    count=sts[["Count", "Alt"]].to_numpy()
    tracemalloc.start()
    start_time = time.time()
    allel.tajima_d(count, positions)
    end_time=time.time()
   # displaying the memory  IN BYTES
    fts_memory.append(tracemalloc.get_traced_memory()[1])
   
   # stopping the library
    tracemalloc.stop()

    calc_time.append(end_time-start_time)
    total_time.append(end_time-tot_time)

    fts_stat_type.append("Allele-based")
    fts_stat.append("tajimas d")
    fts_model.append("multilocus")
    dt.append("Classical")
  

smdata = {'stat_type': fts_stat_type, 'data_type':dt, 'model': fts_model,'calc_time': calc_time,'total_time': total_time, 'stat': fts_stat, 'memory': fts_memory}
smbench_data = pd.DataFrame(smdata)
smbench_data.to_csv("{0}/analysis_benchmarking.txt".format(tmpdir), index=False, header = True, sep = "\t")



    
    
    
    

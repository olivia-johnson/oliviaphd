import os
import sys
import msprime
import pyslim
import numpy as np
import tskit
import pandas as pd
import time
import itertools
#import allel
import statistics
#import scipy
sys.path.insert(1, '/localscratch/oliviaphd/hpc/')
import NCD


def simulate_alleles(tmpdir, group, sim_run, s_pop, w_pop, l, y, rGen, fitness_on, sum_gen, win_gen):
        genomeSize = l

        ## FORWARD SIMULATION
        # for when uneven seasons " -d g_s=" + str(sum_gen)+ " -d g_w=" + str(win_gen)
        start_time = time.time()
        tmpdir_call = "tmpdir='" + str(tmpdir)+ "'"

        cmd = 'slim -d "' +str(tmpdir_call)+ '" -d fit='+ str(fitness_on)+" -d group=" + str(group) + " -d g_s=" + str(sum_gen)+" -d g_w=" + str(win_gen)+" -d sim_run=" + str(sim_run) + " -d GenomeSize=" + str(int(genomeSize)) + " -d L=" + str(l)+ " -d n_s=" + str(int(s_pop)) + " -d n_w=" + str(int(w_pop)) + " -d y=" + str(y) + " -d rGen="+ str(rGen) +" ~/oliviaphd/hpc/chp4/witt_complex_allele.slim"

        print(cmd)
        os.system(cmd)
        # #print("Time for SLiM sim = ", (time.time()- start_time))
        print("Simulations took ", (time.time()-start_time) ,  " seconds")

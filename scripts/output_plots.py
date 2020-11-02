#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 11:49:54 2020

@author: olivia
"""
import os
os.chdir("/Users/olivia/oliviaphd/")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

## import output files from seglift_ts
    ## may need to add in parameter identifier
    
filenames = sorted(glob.glob("~/oliviaphd/data/seglift_ts/sim_s_stat_*.txt"))
output = np.loadtxt("~/oliviaphd/data/seglift_ts/sim_s_stat_*.txt", delimiter = ",")

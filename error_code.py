#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 13 11:35:44 2020

@author: olivia
"""

import msprime
import pyslim

burnin = msprime.simulate(sample_size=1000,Ne=500, length=100000, mutation_rate=1e-04, recombination_rate=1e-8)
burnin_ts = pyslim.annotate_defaults(burnin, model_type="WF", slim_generation=1)
burnin_ts.dump("burnin.trees")
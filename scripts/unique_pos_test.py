#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 19 15:28:23 2021

@author: olivia
"""




test_arr = [1,2,10,10,20]
for i, site in enumerate(test_arr):
    print (i, site)
    if site==test_arr[i+1]:
        test_arr[i+1]=test_arr[i+1]+1
        print(test_arr)
        

current_pos = 0



site_arr = np.rint(5e6*n_met['mut_pos'].values)

np.size(site_arr)

np.size(np.unique(site_arr))

for i in range(len(site_arr)-1):
    if site_arr[i]>=site_arr[i+1]:
        site_arr[i+1] = site_arr[i] + 1

np.max(site_arr)
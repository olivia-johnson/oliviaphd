#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 10:01:44 2021

@author: olivia
"""

rows_list_h=[]
for i in slim_ts.nodes():
    node_id = i.id
    hap = list()
    for p in np.unique(mut_met.mut_pos):
        mut = slim_ts.mutation_at(node_id, p,time=None)
        if mut >= 0:
            allele = 1
        else:
            allele = 0
        hap.append(allele)
    dict_h={}
    dict_h.update({"node_id": node_id})
    dict_h.update({"h": hap})
    rows_list_h.append(dict_h)
    
haplotypes =pd.DataFrame(rows_list_h) 

h1=np.array(haplotypes.haplotypes[node_id==ind.nodes[0]])
h2=np.array(haplotypes.haplotypes[node_id==ind.nodes[1]])
        
    
        
        
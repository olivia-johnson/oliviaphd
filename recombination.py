import pandas as pd
import msprime

def recombination_map(sim_type, group, l, nChrom, chromSize, recRate): 
    ## input group (parameter set identifier, number of chromosomes, chromosome size, recombination rate)
    
    ## GENERATE CHROMOSMES ##
## create recombination map for msprime and slim to simulate unlinked chromsomes
    rec_rows=[]
        # cycle through each chromosome
    for c in range(nChrom): 
        print(c)
        rec_dict = {}
            # generate start position of chromsome
        rec_dict.update({"positions": c*chromSize}) 
            # assign recombination rate of chromosome
        rec_dict.update({"rates":recRate})
        rec_rows.append(rec_dict)
        rec_dict = {}
            # assign end position of chromosome
        rec_dict.update({"positions": int((c+1)*chromSize-1)})
            # assign recombination rate of 0.5 to create break between chromosomes
        rec_dict.update({"rates":0.5})
        rec_rows.append(rec_dict)
    
    if sim_type == "seglift_l10": ## add additional seasonal loci to the end so loci contribute to fitness
            rec_dict = {}
            # assign end position of site
            rec_dict.update({"positions": (nChrom*chromSize+l-10-1)})
            # assign recombination rate of 0.5 to create break between chromosomes
            rec_dict.update({"rates":0.5})
            rec_rows.append(rec_dict)
    
    rec_data = pd.DataFrame(rec_rows)
    
    
## generate slim recombination map
        # reformat to comply with recombinate map required in slim 
    slim_rec = rec_data.positions[1:]
    slim_rec = slim_rec.reset_index(drop=True)
    slim_rec= pd.concat([slim_rec,rec_data.rates[0:-1]],axis=1)
    
    
## formulate msprime recombination map
    rec_data.positions.iloc[-1]=rec_data.positions.iloc[-1]+1
    rec_map = msprime.RecombinationMap(positions = list(rec_data.positions.astype(int)), rates= list(rec_data.rates), num_loci = int(nChrom*chromSize-1))
    
## output slim recombiation map to text file to be read into forward slim simulation
    slim_rec.to_csv("./rec_map.txt", index=False, header = False, sep = "\t")
    
    return rec_map


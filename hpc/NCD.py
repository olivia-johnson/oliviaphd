import numpy as np
import math
# NCD
def ncd(MAF, TF):
# NCD like Std Dev that measures dispersion around a target Frequency.
    ncd = math.sqrt(sum(TF-MAF)**2/len(MAF)) 
    return(ncd)

# window1_MAFs =np.array([0.1,0.25,0.5,0.15,0.35,0.45])
# window2_MAFs =np.array([0.01,0.025,0.05,0.015,0.035,0.045])
# window3_MAFs =np.array([0.35,0.45,0.5,0.45,0.5,0.45])

# target_freq =0.5 # probably set this to expect equilibrium frequency?
# ncd(window1_MAFs ,target_freq)

# ncd(window2_MAFs ,target_freq)

# ncd(window3_MAFs ,target_freq)

# # Thus, lower NCD -> more like balancing selection.
# ncd(window2_MAFs ,TF=statistics.mean(window2_MAFs))

# EEck! that shows one issue with NCD -> you have to choose sensible TF/ Values a priori consider more likely to be BalSel than neutral
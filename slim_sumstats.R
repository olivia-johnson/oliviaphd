library(popgen.tools)
#setwd("~/work/PhD/Drosophila/")

## read in SLiM output (MS)
raw <- scan(
  file = "data/slim_fs_output.txt",
  what = character(0),
  sep = "\n",
  quiet = TRUE
)

## Sep genotypes 
genotypeM = t(sapply(raw[4:length(raw)], function(x) as.numeric(strsplit(x, split = "")[[1]]), USE.NAMES = F))

#extract positions of seg sites
#pos = unlist(strsplit(raw[3], split=" "))
#p = list(pos[2:length(pos)])

# pos needs to be a vector not a list.
pos = as.numeric(strsplit(substr(raw[3], 
                                     12, 
                                     nchar(raw[3]) - 1), 
                              split =" ")[[1]])


# calculate sum_stats using popgen.tools
obs <- sim_obj(cmd = NA,
               seeds = NA,
               segsites = ncol(genotypeM),
               positions = pos,
               genome_matrix = genotypeM,
               sweep ="partial",
               select_coeff = 0.3 )

obs_ss = sum_stats(
  sim = obs,
  nwins = 2,
  ID = 1, 
  split_type = "mut",
  snp=100, #trim genome matrix
)

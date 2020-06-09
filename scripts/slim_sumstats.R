library(popgen.tools)
setwd("~/oliviaphd")

# run sim
cmd = "slim ~/oliviaphd/scripts/sim_fs.slim"
slim_out <- system(cmd, intern = TRUE)


## Sep genotypes 
out_pos = grep("positions", slim_out)
genotypeM = t(sapply(slim_out[(out_pos+1):length(slim_out)], function(x) as.numeric(strsplit(x, split = "")[[1]]), USE.NAMES = F))

#extract positions of seg sites

pos = as.numeric(strsplit(substr(slim_out[out_pos], 
                                     12, 
                                     nchar(slim_out[out_pos]) - 1), 
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

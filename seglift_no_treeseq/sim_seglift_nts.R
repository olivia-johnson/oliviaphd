 library(data.table)
 library(PopGenome)

setwd("~/oliviaphd/seglift_no_treeseq/")
 

group=3
runs = 5
genomeSize = 1e6
popnSize = 1e4
mutRate = 1e-6
recRate = 1e-8
l = 10
y = 1.5
d = 0.6
nWin = 20
sum_gen = 8 #no. summer generations
win_gen = 3 #no. winter generations

for (x in 1:runs){
# cmd = paste0("slim -d n_burnin=" , x," -d GenomeSize=" , genomeSize , " -d N=" , popnSize , " -d mut=" +mutRate+
#                " -d rr=" + recRate +" ~/oliviaphd/seglift_no_treeseq/sim_burnin.slim")
# system(cmd)

slim_cmd = paste0("slim -d n_burnin=0 -d fit=0 -d group=",group," -d sim_run=" , x," -d GenomeSize=" , genomeSize , 
                    " -d L=" , l, " -d N=" , popnSize , " -d y=" , y , " -d d=" , d , " -d mut=" , mutRate , 
                    " -d rr=" , recRate ," ~/oliviaphd/seglift_no_treeseq/seglift_nts.slim")
system(slim_cmd)
}

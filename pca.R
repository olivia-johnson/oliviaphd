# pca #
library(data.table)
require(ggbiplot)
library(ggpubr)


groups=c(1:8) ##positive groups
distance=50

sim_type="soft"
setwd(paste0("~/phd_data/Results/single_locus/", sim_type, "/"))
load(file=paste0("~/Box Sync/data/single_locus_",sim_type,".RData"))
load(file=paste0("~/Box Sync/data/af_single_locus_",sim_type,".RData"))
sum_stats[, dist:=n_win-250]

##soft and hard
fix_times = freq_data[, sum(mut_freq)==1, by=c("Gen","group", "run")][V1==TRUE, min(Gen), by=c("group", "run")]
data_s=NULL
 start_values=c(cumprod(c(20, rep.int(2, 9))), 12000, 24000, 36000,48000,50000,62000,74000,86000,980000)
values=sort(unique(freq_data$Gen))
for (i in 1:length(groups)){
  for (j in 1:20){
  fixed=fix_times[group==i & run==j, V1]
  if (fixed %in% values){
    if ((fixed+20) %in% values &length(values[values>=fixed&values<(fixed+20)])==20){
      
      sum_stat=sum_stats[group==i & run==j & Gen>=fixed & Gen<fixed+20 &edges==F &dist>=-distance & dist<=distance]
    }else{
      val=start_values[start_values>fixed]
      nearest=min(val)
sum_stat=sum_stats[group==i & run==j & Gen>=nearest & Gen<nearest+20 &edges==F &dist>=-distance & dist<=distance]
    }}else {
      val=values[values>fixed+20]
      nearest=min(val)
sum_stat=sum_stats[group==i & run==j & Gen>=nearest & Gen<nearest+20 &edges==F &dist>=-distance & dist<=distance]
      
    }
    # print(paste(fixed, min(sum_stat$Gen),length(unique(sum_stat$Gen))))
  data_s=rbind(data_s,sum_stat)
  }
}

data_s[, sim_type:="S"]
data_s[, id_pca:=paste(sim_type, group, run, Gen, sep="_"), by=c("group", "run", "Gen")]
data_s=data_s[, .(id_pca,pop_season, season, gen_year, sim_type, n_win, seg_sites, diversity, theta_w_allele,tajimas_d_branch, tajimas_d_allel, r2, haplotype_diversity, H1, H12, H123, H2H1)]

sim_type="hard"
setwd(paste0("~/phd_data/Results/single_locus/", sim_type, "/"))
load(file=paste0("~/Box Sync/data/single_locus_",sim_type,".RData"))
load(file=paste0("~/Box Sync/data/af_single_locus_",sim_type,".RData"))
sum_stats[, dist:=n_win-250]

##soft and hard
fix_times = freq_data[, sum(mut_freq)==1, by=c("Gen","group", "run")][V1==TRUE, min(Gen), by=c("group", "run")]
data_h=NULL
start_values=c(cumprod(c(20, rep.int(2, 9))), 12000, 24000, 36000,48000,50000,62000,74000,86000,980000)
values=sort(unique(freq_data$Gen))
for (i in 1:8){
  for (j in 1:20){
    fixed=fix_times[group==i & run==j, V1]
    if (fixed %in% values){
      if ((fixed+20) %in% values &length(values[values>=fixed&values<(fixed+20)])==20){
        
        sum_stat=sum_stats[group==i & run==j & Gen>=fixed & Gen<fixed+20 &edges==F &dist>=-distance & dist<=distance]
      }else{
        val=start_values[start_values>fixed]
        nearest=min(val)
        sum_stat=sum_stats[group==i & run==j & Gen>=nearest & Gen<nearest+20 &edges==F &dist>=-distance & dist<=distance]
      }}else {
        val=values[values>fixed+20]
        nearest=min(val)
        sum_stat=sum_stats[group==i & run==j & Gen>=nearest & Gen<nearest+20 &edges==F &dist>=-distance & dist<=distance]
        
      }
    # print(paste(fixed, min(sum_stat$Gen),length(unique(sum_stat$Gen))))
    data_h=rbind(data_h,sum_stat)
  }
}

data_h[, sim_type:="H"]
data_h[, id_pca:=paste(sim_type, group, run, Gen, sep="_"), by=c("group", "run", "Gen")]
data_h=data_h[, .(id_pca,pop_season, season, gen_year, sim_type, n_win, seg_sites, diversity,theta_w_allele, tajimas_d_branch, tajimas_d_allel, r2, haplotype_diversity, H1, H12, H123, H2H1)]


sim_type="wittmann_unlinked"
setwd(paste0("~/phd_data/Results/single_locus/", sim_type, "/"))
load(file=paste0("~/Box Sync/data/single_locus_",sim_type,".RData"))
# load(file=paste0("~/Box Sync/data/af_single_locus_",sim_type,".RData"))
sum_stats[, dist:=n_win-250]

data_WE=sum_stats[ Gen >=12040 & Gen<12060&edges==F &dist>=-distance & dist<=distance]
data_WE[, sim_type:="WE"]
data_WE[, id_pca:=paste(sim_type, group, run, Gen, sep="_"), by=c("group", "run", "Gen")]
data_WE=data_WE[, .(id_pca,pop_season, season, gen_year, sim_type, n_win, seg_sites, diversity,theta_w_allele, tajimas_d_branch, tajimas_d_allel, r2, haplotype_diversity, H1, H12, H123, H2H1)]

data_WL=sum_stats[ Gen >=96040 & Gen<96080 & edges==F & dist>=-distance & dist<=distance]
data_WL[, sim_type:="WL"]
data_WL[, id_pca:=paste(sim_type, group, run, Gen, sep="_"), by=c("group", "run", "Gen")]
data_WL=data_WL[, .(id_pca,pop_season, season, gen_year, sim_type, n_win, seg_sites, diversity,theta_w_allele, tajimas_d_branch, tajimas_d_allel, r2, haplotype_diversity, H1, H12, H123, H2H1)]

sim_type="balancing"
setwd(paste0("~/phd_data/Results/single_locus/", sim_type, "/"))
load(file=paste0("~/Box Sync/data/single_locus_",sim_type,".RData"))
# load(file=paste0("~/Box Sync/data/af_single_locus_",sim_type,".RData"))
sum_stats[, dist:=n_win-250]

data_BE=sum_stats[ Gen >=5120 & Gen<5140&edges==F &dist>=-distance & dist<=distance]
data_BE[, sim_type:="BE"]
data_BE[, id_pca:=paste(sim_type, group, run, Gen, sep="_"), by=c("group", "run", "Gen")]
data_BE=data_BE[, .(id_pca,pop_season, season, gen_year, sim_type, n_win, seg_sites, diversity,theta_w_allele, tajimas_d_branch, tajimas_d_allel, r2, haplotype_diversity, H1, H12, H123, H2H1)]

data_BL=sum_stats[ Gen >=96040 & Gen<96080&edges==F &dist>=-distance & dist<=distance]
data_BL[, sim_type:="BL"]
data_BL[, id_pca:=paste(sim_type, group, run, Gen, sep="_"), by=c("group", "run", "Gen")]
data_BL=data_BL[, .(id_pca,pop_season, season, gen_year, sim_type, n_win, seg_sites, diversity,theta_w_allele, tajimas_d_branch, tajimas_d_allel, r2, haplotype_diversity, H1, H12, H123, H2H1)]


data=rbind(data_s, data_h, data_WE, data_WL, data_BE, data_BL)
data[, tajimas_d_allel:=NULL]
remove(data_s, data_h, data_WE, data_WL, data_BE, data_BL, sum_stats, sum_stat, freq_data, fix_times)

# pca_data_original=data[, .(id_pca, n_win, seg_sites, diversity, tajimas_d_branch, tajimas_d_allel, r2, haplotype_diversity, H1, H12, H123, H2H1)]
pca_cols=c("seg_sites", "diversity", "tajimas_d_branch","theta_w_allele",  "r2", "haplotype_diversity", "H1", "H12", "H123", "H2H1")
pca_data=dcast(data[pop_season=="FP" & sim_type!="H" & sim_type!="S"], id_pca + sim_type+season+gen_year ~ n_win, value.var = pca_cols) ### pca file
id_cols=c("id_pca", "sim_type","season","gen_year")

pca_out = princomp(pca_data[,.SD, .SDcols = !id_cols])

# pca_data[rowSums(is.na(pca_data)) > 0]

pdf("~/phd_data/Results/single_locus/plots/pca_scree_FP_noPS.pdf", width = 6, height=5)
scree=screeplot(pca_out, type = c("lines"))
dev.off()

plot(pca_out$scores[,1:2], col=factor(pca_data$sim_type))

pca12=ggbiplot(pca_out, 1:2, var.axes = FALSE) + geom_point(aes(col=factor(pca_data$sim_type)), alpha=0.8)+coord_fixed(ratio=1)+
  scale_x_continuous(breaks=-10:10)+scale_y_continuous(breaks=-10:10)+ggtitle("FP_noPS")
ggexport(pca12, filename="~/phd_data/Results/single_locus/plots/pca_1_2_FP_noPS.pdf", width = 10, height=10)

pca23=ggbiplot(pca_out, 2:3, var.axes = FALSE) + geom_point(aes(col=factor(pca_data$sim_type)), alpha=0.8)+coord_fixed(ratio=1)+
  scale_x_continuous(breaks=-10:10)+scale_y_continuous(breaks=-10:10)+ggtitle("CP10K_noPS")
ggexport(pca23, filename="~/phd_data/Results/single_locus/plots/pca_2_3_CP10K_noPS.pdf", width = 10, height=10)

pca34=ggbiplot(pca_out, 3:4, var.axes = FALSE) + geom_point(aes(col=factor(pca_data$sim_type)), alpha=0.8)+coord_fixed(ratio=1)+
  scale_x_continuous(breaks=-10:10)+scale_y_continuous(breaks=-10:10)+ggtitle("CP10K_noPS")
ggexport(pca34, filename="~/phd_data/Results/single_locus/plots/pca_3_4_CP10K_noPS.pdf", width = 10, height=10)

pca13=ggbiplot(pca_out, c(1,3), var.axes = FALSE) + geom_point(aes(col=factor(pca_data$sim_type)), alpha=0.8)+coord_fixed(ratio=1)+
scale_x_continuous(breaks=-10:10)+scale_y_continuous(breaks=-10:10)+ggtitle("CP10K_noPS")
ggexport(pca13, filename="~/phd_data/Results/single_locus/plots/pca_1_3_CP10K_noPS.pdf", width = 10, height=10)


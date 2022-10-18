## single_locus combined
library(data.table)
library(ggplot2)
library(viridis)
library(stringr)
library(egg)
library(dplyr)
library(ggpubr)
library(psych)

sim_type="wittmann_unlinked"
setwd(paste0("~/phd_data/Results/single_locus/", sim_type, "/"))
load(file=paste0("~/Box Sync/data/single_locus_",sim_type,".RData"))
load(file=paste0("~/Box Sync/data/af_single_locus_",sim_type,".RData"))

sum_stats[gen_season == "EG", gen_year:=Gen%%20, by=c("label", "Gen")]
sum_stats[gen_season == "UG", gen_year:=Gen%%12, by=c("label", "Gen")] 
sum_stats[, dist:=n_win-250]

sums=sum_stats[Gen>60000 & edges==F, .(taj_db=mean(tajimas_d_branch), 
                            theta_w=mean(theta_w_allele),
                            div=mean(diversity), 
                            r2=mean(r2), 
                            hap_div=mean(haplotype_diversity), 
                            H1=mean(H1),
                            H12=mean(H12),
                            H2H1=mean(H2H1),
                            ehh=mean(ehh),
                            taj_da=mean(tajimas_d_allel), burninNe, exp, midpoint, block, dist, s_s), 
               by=c("group", "label", "n_win", "gen_year", "run", "season")]
genomeSize=unique(sum_stats$genomeSize)
sums=unique(sums)
sums[,linkage:=ifelse(midpoint>genomeSize, "unlinked", "linked"), by=c("group", "n_win")]
sums[, win_10k:=n_win%/%100, by=c("group", "n_win")]
sums[, mean_tajd:=mean(taj_da), by=c("group", "midpoint")]
sums[, mean_div:=mean(div), by=c("group", "midpoint")]
sums[, mean_thetaw:=mean(theta_w), by=c("group", "midpoint")]
sums[, mean_r2:=mean(r2), by=c("group", "midpoint")]
sums[, mean_hapdiv:=mean(hap_div), by=c("group", "midpoint")]
sums[, mean_ehh:=mean(ehh), by=c("group", "midpoint")]
sums[, mean_H1:=mean(H1), by=c("group", "midpoint")]
sums[, mean_H12:=mean(H12), by=c("group", "midpoint")]
sums[, mean_H2H1:=mean(H2H1), by=c("group", "midpoint")]
sums[, smean_tajd:=mean(taj_db), by=c("group", "midpoint", "season")]
sums[, smean_div:=mean(div), by=c("group", "midpoint", "season")]
sums[, smean_thetaw:=mean(theta_w), by=c("group", "midpoint", "season")]
sums[, smean_r2:=mean(r2), by=c("group", "midpoint", "season")]
sums[, smean_hapdiv:=mean(hap_div), by=c("group", "midpoint", "season")]
sums[, smean_ehh:=mean(ehh), by=c("group", "midpoint", "season")]
sums[, smean_H1:=mean(H1), by=c("group", "midpoint", "season")]
sums[, smean_H12:=mean(H12), by=c("group", "midpoint", "season")]
sums[, smean_H2H1:=mean(H2H1), by=c("group", "midpoint", "season")]
sums[, linkage_tajd:=mean(taj_da), by=c("group", "linkage")]
sums[, linkage_div:=mean(div), by=c("group", "linkage")]
sums[, linkage_thetaw:=mean(theta_w), by=c("group", "linkage")]
sums[, linkage_r2:=mean(r2), by=c("group", "linkage")]
sums[, linkage_hapdiv:=mean(hap_div), by=c("group", "linkage")]
sums[, linkage_ehh:=mean(ehh), by=c("group", "linkage")]
sums[, linkage_H1:=mean(H1), by=c("group", "linkage")]
sums[, linkage_H12:=mean(H12), by=c("group", "linkage")]
sums[, linkage_H2H1:=mean(H2H1), by=c("group", "linkage")]
sums[, rgn:=paste0(run, gen_year), by=c("group", "run", "gen_year")]

EG_gens = c(4, 9, 14, 19)

groups = c(1:6) ## neutral groups
groups= c(1:4, 6:11,14:20) ##wittmann groups
g=6

distance=20
for (g in groups){
div_plot = ggplot(sums[group==g], aes(x = dist))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = div, group=rgn),alpha = 0.2)+
  geom_hline(data = sums[group==g], aes(yintercept = exp), col = "grey2")+ 
  geom_line(data = sums[group==g],aes(x=dist, y = mean_div), col="orange")+
  geom_hline(data = subset(sums[group==g], linkage=="unlinked"), aes(yintercept = linkage_div), col = "hotpink", size=1)+
  ylab("Diversity")+
  # xlim(0, genomeSize/1000)+
 xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  theme(legend.position = "none")

taj_plot = ggplot(sums[group==g], aes(x = dist, y = taj_db))+
  geom_vline(data = subset(sums[group==g]),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = taj_db, group=rgn),alpha = 0.1)+
  geom_hline(data = sums[group==g], aes(yintercept = 0), col = "grey")+ 
  geom_line(data = sums[group==g],aes(x=dist, y = mean_tajd), col="orange")+
  geom_hline(data = subset(sums[group==g], linkage=="unlinked"), aes(yintercept = linkage_tajd), col = "hotpink", size=1)+
  ylab("Tajima's D")+
  # xlim(0, genomeSize/1000)+ 
   xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  theme(legend.position = "none")

r2_plot = ggplot(sums[group==g], aes(x = (dist), y = r2))+
  geom_vline(data = subset(sums[group==g]),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = r2, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_r2), col="orange")+ 
  geom_hline(data = subset(sums[group==g], linkage=="unlinked"), aes(yintercept = linkage_r2), col = "hotpink", size=1)+
  ylab("R^2")+
  # xlim(0, genomeSize/1000)+
 xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  theme(legend.position = "none")

thetaw_plot = ggplot(sums[group==g], aes(x = (dist), y = theta_w))+
  geom_vline(data = subset(sums[group==g]),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = theta_w, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_thetaw), col="orange")+ 
  geom_hline(data = subset(sums[group==g], linkage=="unlinked"), aes(yintercept = linkage_thetaw), col = "hotpink", size=1)+
  ylab("Theta W")+
 xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  theme(legend.position = "none")

hapdiv_plot = ggplot(sums[group==g], aes(x = dist, y = hap_div))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = hap_div, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_hapdiv), col="orange")+ 
  geom_hline(data = subset(sums[group==g], linkage=="unlinked"), aes(yintercept = linkage_hapdiv), col = "hotpink", size=1)+
  ylab("Haplotype Diversity")+
  # xlim(0, genomeSize/1000)+
 xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  theme(legend.position = "none")

ehh_plot = ggplot(sums[group==g], aes(x = dist, y = ehh))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = ehh, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_ehh), col="orange")+ 
  geom_hline(data = subset(sums[group==g], linkage=="unlinked"), aes(yintercept = linkage_ehh), col = "hotpink", size=1)+
  ylab("EHH")+
  # xlim(0, genomeSize/1000)+  
 xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  theme(legend.position = "none")

H1_plot = ggplot(sums[group==g], aes(x = (dist), y = H1))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = H1, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_H1), col="orange")+ 
  geom_hline(data = subset(sums[group==g], linkage=="unlinked"), aes(yintercept = linkage_H1), col = "hotpink", size=1)+
  ylab("H1")+
  # xlim(0, genomeSize/1000)+
 xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  theme(legend.position = "none")
H12_plot = ggplot(sums[group==g], aes(x = (dist), y = H12))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = H12, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_H12), col="orange")+ 
  geom_hline(data = subset(sums[group==g], linkage=="unlinked"), aes(yintercept = linkage_H12), col = "hotpink", size=1)+
  ylab("H12")+
  # xlim(0, genomeSize/1000)+
 xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  theme(legend.position = "none")

H2H1_plot = ggplot(sums[group==g], aes(x = dist, y = H2H1))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = H2H1, group=rgn),alpha = 0.1)+
  geom_line(data = sums[group==g],aes(x=dist, y = mean_H2H1), col="orange")+ 
  geom_hline(data = subset(sums[group==g], linkage=="unlinked"), aes(yintercept = linkage_H2H1), col = "hotpink", size=1)+
  ylab("H2/H1")+
  # xlim(0, genomeSize/1000)+
 xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  theme(legend.position = "none")

allele_plot =ggplot(data=freq_data[Gen>95000 & Gen<97000& group==g],aes(x = Gen, y= mut_freq))+
  geom_line(aes(col=label, group=block),alpha = 0.3) +
  ggtitle(unique(sums[group==g, label])) +
  xlab("Generation")+
  scale_y_continuous(limits = c(0, 1.0))+
  ylab("Allele Frequency")+
  theme(legend.position = "none")+
  labs(col= "Seasonal Mutation")

# allele_plot =ggplot(data=freq_data[group==g],aes(x = Gen, y= mut_freq))+
#   geom_line(aes(col=label, group=block),alpha = 0.3) +
#   ggtitle(unique(sums[group==g, label])) +
#   xlab("Generation")+
#   scale_y_continuous(limits = c(0, 1.0))+
#   ylab("Allele Frequency")+
#   theme(legend.position = "none")+
#   labs(col= "Seasonal Mutation")

plot=ggarrange(allele_plot, div_plot, taj_plot, thetaw_plot, r2_plot, hapdiv_plot, ehh_plot, H1_plot, H12_plot, H2H1_plot)
ggexport(plot, filename=paste0("plots/combined_rgn_", unique(sums[group==g, label]), "_", distance,".pdf"), width = 15, height=10)
}


div_plot = ggplot(data, aes(x = dist,y = div))+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = div, group=block),alpha = 0.2)+
  geom_smooth(method="loess", se=F)+
  geom_hline(data = data, aes(yintercept = exp), col = "grey2")+ 
  geom_hline(data = subset(data, linkage=="unlinked"), aes(yintercept = linkage_div), col = "hotpink", size=1)+
  facet_wrap(~gen_year, drop = TRUE)+
  ylab("Diversity")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  theme(legend.position = "none")

div_plot = ggplot(sums, aes(x = dist))+
     geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
     geom_line(aes(x = dist, y = div, group=rgn),alpha = 0.2)+
     geom_hline(data = sums, aes(yintercept = exp), col = "grey2")+ 
     geom_line(data = sums,aes(x=dist, y = mean_div), col="orange")+
     geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_div), col = "hotpink", size=1)+
     ylab("Diversity")+
     # xlim(0, genomeSize/1000)+
     xlim(-distance, distance) +
     facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
     xlab("Distance from selected site (10kb)")+
     theme(legend.position = "none")
 ggexport(div_plot, filename="plots/div_test.pdf", width = 15, height=10)

div_plot = ggplot(sums[gen_year==5|gen_year==10|gen_year==15|gen_year==0], aes(x = dist))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = div, group=rgn, col=as.factor(gen_year)),alpha = 0.5)+
  geom_hline(data = sums, aes(yintercept = exp), col = "grey2")+ 
  geom_line(data = sums,aes(x=dist, y = mean_div), col="orange")+
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_div), col = "hotpink", size=1)+
  ylab("Diversity")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  xlab("Distance from selected site (10kb)")
ggexport(div_plot, filename="plots/div_gy.pdf", width = 15, height=10)

taj_plot = ggplot(sums, aes(x = dist))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = taj_db, group=rgn),alpha = 0.2)+
  geom_hline(data = sums, aes(yintercept = 0), col = "grey2")+ 
  geom_line(data = sums,aes(x=dist, y = mean_tajd), col="orange")+
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_tajd), col = "hotpink", size=1)+
  ylab("Tajima's D")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  xlab("Distance from selected site (10kb)")+
  theme(legend.position = "none")
ggexport(taj_plot, filename="plots/tajdb.pdf", width = 15, height=10)

r2_plot = ggplot(sums, aes(x = (dist), y = r2))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = r2, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_r2), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_r2), col = "hotpink", size=1)+
  ylab("R^2")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(r2_plot, filename="plots/r2.pdf", width = 15, height=10)

thetaw_plot = ggplot(sums, aes(x = (dist), y = theta_w))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = theta_w, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_thetaw), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_thetaw), col = "hotpink", size=1)+
  ylab("Theta W")+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(thetaw_plot, filename="plots/thetaw_plot.pdf", width = 15, height=10)

hapdiv_plot = ggplot(sums, aes(x = dist, y = hap_div))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = hap_div, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_hapdiv), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_hapdiv), col = "hotpink", size=1)+
  ylab("Haplotype Diversity")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(hapdiv_plot, filename="plots/hapdiv_plot.pdf", width = 15, height=10)

ehh_plot = ggplot(sums, aes(x = dist, y = ehh))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = ehh, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_ehh), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_ehh), col = "hotpink", size=1)+
  ylab("EHH")+
  # xlim(0, genomeSize/1000)+  
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(ehh_plot, filename="plots/ehh_plot.pdf", width = 15, height=10)

H1_plot = ggplot(sums, aes(x = (dist), y = H1))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = H1, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_H1), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_H1), col = "hotpink", size=1)+
  ylab("H1")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(H1_plot, filename="plots/H1_plot.pdf", width = 15, height=10)

H12_plot = ggplot(sums, aes(x = (dist), y = H12))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = H12, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_H12), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_H12), col = "hotpink", size=1)+
  ylab("H12")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(H12_plot, filename="plots/H12_plot.pdf", width = 15, height=10)

H2H1_plot = ggplot(sums, aes(x = dist, y = H2H1))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = H2H1, group=rgn),alpha = 0.1)+
  geom_line(data = sums,aes(x=dist, y = mean_H2H1), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_H2H1), col = "hotpink", size=1)+
  ylab("H2/H1")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(H2H1_plot, filename="plots/H2H1_plot.pdf", width = 15, height=10)




data=sum_stats[Gen>60000 & h_s==0.6 & pop_season=="CP10k", .(taj_db=mean(tajimas_d_branch), 
                            theta_w=mean(theta_w_allele),
                            div=mean(diversity), 
                            r2=mean(r2), 
                            hap_div=mean(haplotype_diversity), 
                            H1=mean(H1),
                            H12=mean(H12),
                            H2H1=mean(H2H1),
                            ehh=mean(ehh),
                            taj_da=mean(tajimas_d_allel)), by=c("group", "label", "gen_year", "midpoint", "burninNe", "exp")]

div = ggplot(data, aes(x=midpoint/1e6, y=div, col=factor(gen_year)))+
  #geom_smooth(method="loess", se=F)+
  geom_line()+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red")+
  geom_hline(data = data, aes(yintercept = exp), col = "grey2")+ 
  facet_grid("label")
ggexport(div, filename="plots/div_new.pdf", width = 15, height=10)

tajda = ggplot(data, aes(x=midpoint/1e6, y=taj_da, col=factor(gen_year)))+
  geom_smooth(method="loess", se=F)+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red")+
  geom_hline(data = data, aes(yintercept = 0), col = "grey2")+ 
  facet_grid("label")
ggexport(tajda, filename="plots/taj_da_new.pdf", width = 15, height=10)

tajdb = ggplot(data, aes(x=midpoint/1e6, y=taj_db, col=factor(gen_year)))+
  geom_smooth(method="loess", se=F)+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red")+
  geom_hline(data = data, aes(yintercept = 0), col = "grey2")+ 
  facet_grid("label")
ggexport(tajdb, filename="plots/taj_db_new.pdf", width = 15, height=10)

theta_w = ggplot(data, aes(x=midpoint/1e6, y=theta_w, col=factor(gen_year)))+
  geom_smooth(method="loess", se=F)+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red")+
  facet_grid("label")
ggexport(theta_w, filename="plots/theta_w_new.pdf", width = 15, height=10)

r2 = ggplot(data, aes(x=midpoint/1e6, y=r2, col=factor(gen_year)))+
  geom_smooth(method="loess", se=F)+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red")+
  facet_grid("label")
ggexport(r2, filename="plots/r2_new.pdf", width = 15, height=10)

hap_div = ggplot(data, aes(x=midpoint/1e6, y=hap_div, col=factor(gen_year)))+
  geom_smooth(method="loess", se=F)+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red")+
  facet_grid("label")
ggexport(hap_div, filename="plots/hap_div_new.pdf", width = 15, height=10)

H1 = ggplot(data, aes(x=midpoint/1e6, y=H1, col=factor(gen_year)))+
  geom_smooth(method="loess", se=F)+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red")+
  facet_grid("label")
ggexport(H1, filename="plots/H1_new.pdf", width = 15, height=10)

H12 = ggplot(data, aes(x=midpoint/1e6, y=H12, col=factor(gen_year)))+
  geom_smooth(method="loess", se=F)+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red")+
  facet_grid("label")
ggexport(H12, filename="plots/H12_new.pdf", width = 15, height=10)

H2H1 = ggplot(data, aes(x=midpoint/1e6, y=H2H1, col=factor(gen_year)))+
  geom_smooth(method="loess", se=F)+
  geom_vline(data = subset(data),aes(xintercept = 0), col= "red")+
  facet_grid("label")
ggexport(H2H1, filename="plots/H2H1_new.pdf", width = 15, height=10)

ggplot(sum_stats[Gen >= 60000 & pop_season == "CP10k" & gen_season == "EG" & h_s == 0.6 & gen_year %in% 0:19], aes(midpoint/1e6, diversity, col=factor(gen_year))) + 
  geom_smooth(span = 0.05, n=200, method = "loess", se=FALSE) + 
  facet_grid(.~label) +
  labs(x="Mb") + theme_bw()



taj_plot = ggplot(sums, aes(x = dist))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = taj_db, group=rgn),alpha = 0.2)+
  geom_hline(data = sums, aes(yintercept = 0), col = "grey2")+ 
  geom_line(data = sums,aes(x=dist, y = mean_tajd), col="orange")+
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_tajd), col = "hotpink", size=1)+
  ylab("Tajima's D")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  xlab("Distance from selected site (10kb)")+
  theme(legend.position = "none")
ggexport(taj_plot, filename="plots/tajdb.pdf", width = 15, height=10)

r2_plot = ggplot(sums, aes(x = (dist), y = r2))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = r2, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_r2), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_r2), col = "hotpink", size=1)+
  ylab("R^2")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(r2_plot, filename="plots/r2.pdf", width = 15, height=10)

thetaw_plot = ggplot(sums, aes(x = (dist), y = theta_w))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = theta_w, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_thetaw), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_thetaw), col = "hotpink", size=1)+
  ylab("Theta W")+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(thetaw_plot, filename="plots/thetaw_plot.pdf", width = 15, height=10)

hapdiv_plot = ggplot(sums, aes(x = dist, y = hap_div))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = hap_div, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_hapdiv), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_hapdiv), col = "hotpink", size=1)+
  ylab("Haplotype Diversity")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(hapdiv_plot, filename="plots/hapdiv_plot.pdf", width = 15, height=10)

ehh_plot = ggplot(sums, aes(x = dist, y = ehh))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = ehh, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_ehh), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_ehh), col = "hotpink", size=1)+
  ylab("EHH")+
  # xlim(0, genomeSize/1000)+  
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(ehh_plot, filename="plots/ehh_plot.pdf", width = 15, height=10)

H1_plot = ggplot(sums, aes(x = (dist), y = H1))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = H1, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_H1), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_H1), col = "hotpink", size=1)+
  ylab("H1")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(H1_plot, filename="plots/H1_plot.pdf", width = 15, height=10)

H12_plot = ggplot(sums, aes(x = (dist), y = H12))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = H12, group=rgn),alpha = 0.1)+
  geom_line(aes(x=dist, y = mean_H12), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_H12), col = "hotpink", size=1)+
  ylab("H12")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(H12_plot, filename="plots/H12_plot.pdf", width = 15, height=10)

H2H1_plot = ggplot(sums, aes(x = dist, y = H2H1))+
  geom_vline(data = subset(sums),aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist, y = H2H1, group=rgn),alpha = 0.1)+
  geom_line(data = sums,aes(x=dist, y = mean_H2H1), col="orange")+ 
  geom_hline(data = subset(sums, linkage=="unlinked"), aes(yintercept = linkage_H2H1), col = "hotpink", size=1)+
  ylab("H2/H1")+
  # xlim(0, genomeSize/1000)+
  xlim(-distance, distance) +
  xlab("Distance from selected site (10kb)")+
  facet_wrap(s_s~label, drop = TRUE, scales = "free_y")+
  theme(legend.position = "none")
ggexport(H2H1_plot, filename="plots/H2H1_plot.pdf", width = 15, height=10)

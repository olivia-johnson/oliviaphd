library(data.table)
library(ggplot2)
library(viridis)
library(stringr)
library(egg)
library(tidyverse)
setwd("~/oliviaphd/seglift_substitute")

groups = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19","20", "21", "22", "23", "24", "25", "26", "27", "28")
groups=c("15", "16", "17", "18", "19","20", "21","22", "23", "24", "25", "26", "27", "28")
windows=51
chromsize=5e5

#ID = "L=10, d = 0.65, y=0.5"
##"L=10, d = 0.65, y=4"
##"No Fitness, L=10"
freq_data = NULL
sum_stats = NULL

for (g in groups){
  print(g)
  
  f_list  <- list.files(path =paste0("~/oliviaphd/seglift_substitute/group_", g, "/"),pattern ="al_freq_")
  
  s_list  <- list.files(path =paste0("~/oliviaphd/seglift_substitute/group_", g, "/"),pattern ="s_stat_")
  
  parameters <- fread(file=paste0("~/oliviaphd/seglift_substitute/group_",g, "/parameters.yml"), sep = ":")
  setkey(parameters, V1)
  
  dom = parameters["d", V2]
  epi = parameters["y", V2]
  fiton=parameters["fitness_on", V2]
  sum_gen =parameters["sum_gen", V2]
  win_gen =parameters["win_gen", V2] 
  sum_pop =parameters["s_pop", V2]
  win_pop =parameters["w_pop", V2]
  if (parameters["slim_sim", V2] == "seglift_substitute"){
    sim_type = "seg"
  } else{
    sim_type = "all"
  }
  ## collate al_freq files
  for (i in 1:(length(f_list ))){  
    filename = f_list [i]
    #print(filename)
    al_freq  = fread(file = paste0("~/oliviaphd/seglift_substitute/group_",g, "/",filename) )
    al_freq [,run := i]
    al_freq [,group:=g]
    al_freq [,d:=dom]
    al_freq [,y:=epi]
    al_freq [,fit:=fiton]
    al_freq [,s_gen:=sum_gen]
    al_freq [,w_gen:=win_gen]
    al_freq [,s_pop:=sum_pop]
    al_freq [,w_pop:=win_pop]
    al_freq[, fit_type:=sim_type]
    freq_data  = rbind(freq_data , al_freq )
  }

  ## collate sum_stat files
  for (i in 1:(length(s_list ))){  
    filename = s_list [i]
    #print(filename)
    s_stat  = fread(file = paste0("~/oliviaphd/seglift_substitute/group_",g, "/",filename))
    s_stat [,run := i]
    s_stat [,group:=g]
    s_stat [,d:=dom]
    s_stat [,y:=epi]
    s_stat [,fit:=fiton]
    s_stat [,s_gen:=sum_gen]
    s_stat [,w_gen:=win_gen]
    s_stat [,s_pop:=sum_pop]
    s_stat [,w_pop:=win_pop]
    s_stat[, fit_type:=sim_type]
    sum_stats = rbind(sum_stats, s_stat )
  }

}

freq_data [, mean_freq:=mean(mut_freq), by=c("Gen", "mut_pos", "group")]
#freq_data [, Season:= ifelse((Gen/2)%%15==0, "Winter", "Summer")]
freq_data [, block :=paste0(group, "_", run, "_",(mut_pos/5e5)+.5)]
freq_data[,symmetry:=ifelse(s_gen == w_gen, "symmetrical", "unsymmetrical")]
freq_data[,label:=ifelse(fit==0, paste0("No Fitness_", symmetry), paste(d, y, fit_type, symmetry, sep="_")),  by="group"]

sum_stats[, win_pos:=n_win%%windows]
sum_stats[, midpoint:= (win_end-win_start)/2 + win_start]
sum_stats[, block :=paste0(group, "_", run, "_",chrom)]
sum_stats[, chrom_pos:= midpoint-(chrom-1)*chromsize]
sum_stats[, dist:=(chromsize/2)-chrom_pos-1, by="n_win"]
# sum_stats[, theta_w_corr:= theta_w*1e15]
sum_stats[, mean_tajd:=mean(tajimas_d_branch), by=c("time", "chrom_pos", "group")]
sum_stats[, mean_div:=mean(diversity), by=c("time", "chrom_pos", "group")]
sum_stats[, mean_thetaw:=mean(theta_w), by=c("time", "chrom_pos", "group")]
sum_stats[,symmetry:=ifelse(s_gen == w_gen, "symmetrical", "unsymmetrical")]
# sum_stats[, mean_thetaw:=mean(theta_w_corr), by=c("time", "win_pos", "group")]
sum_stats[,label:=ifelse(fit==0, paste0("No Fitness_", symmetry), paste(d, y, fit_type, symmetry, sep="_")),  by="group"]



plot = ggplot(freq_data [group ==2 & run ==2],aes( x = Gen, y= mut_freq))+
  geom_line(aes( col = as.factor(mut_pos))) +
  ggtitle("Allele Frequency (Single Simulation, N = 1000, d=0.65)") +
  xlab("Generation")+
  ylab("Allele Frequency")

##ggsave(filename =("freq_plot_0_3.jpg"), plot = plot, width = 15, height = 10)

##group.labs <- c("1"="Rel-No Fitness", "6"= "Rel-y=4, d=0.65", "7"="Rel-y=0.5, d=0.65","5"= "Rel-y=1, d=0.65", "3"= "Rel-y=4, d=0.5","4"="Rel-y=0.5, d=0.5", "2"="Rel-y=1, d=0.5","8"="Abs-No Fitness", "13"= "Abs-y=4, d=0.65", "14"="Abs-y=0.5, d=0.65","12"= "Abs-y=1, d=0.65", "10"= "Abs-y=4, d=0.5","11"="Abs-y=0.5, d=0.5", "9"="Abs-y=1, d=0.5" )

#al_time = c( 1,13,10000,10013,20000,20013,30000,30013,40000,40013,50000,50013,60000, 60013,70000,70013,80000, 80013,90000,90013,100000, 100013, 110000, 110013,120000)
al_time=c(20000, 20007, 20013, 20014)
allele_plot =ggplot(freq_data[Gen<50000] ,aes(x = Gen, y= mut_freq, group = block))+
  geom_line(aes(col=as.factor(group)),alpha = 0.3) +
  facet_wrap(~label)+
  #ggtitle(paste0("Allele Frequency (", ID, ")")) +
  xlab("Generation")+
  scale_y_continuous(limits = c(0, 1.0))+
  ylab("Allele Frequency")+
  theme(legend.position = "none")+
  labs(col= "Seasonal Mutation")
ggsave(filename =paste0("plots/allele_freq_relvsabsfit.jpg"), plot = allele_plot , width = 15, height = 10)


freq_data[mut_freq!=1, Freq.bin:="Segregating"]
freq_data[mut_freq==1, Freq.bin:="Fixed_Summer"]
ff <- freq_data[, .N, by = c("Gen", "group", "Freq.bin")]
gg <- freq_data[, .(Freq.bin="Fixed_Winter", N=80-.N), by = c("Gen", "group")]

xx <- rbind(ff, gg)

freq_bin = ggplot(xx[Gen<50000], aes(x=Gen, y=N, col = Freq.bin))+
  geom_line()+
  facet_wrap(~factor(group, levels=c("1", "8", "2", "9", "3", "10", "4", "11", "5", "12", "6", "13", "7", "14")), labeller = as_labeller(group.labs))
ggsave(filename =paste0("plots/freq_bin_relvsabsfit.jpg"), plot = freq_bin , width = 15, height = 10)


freq_data[mut_freq!=1, Freq.bin:="Segregating"]
freq_data[mut_freq==1, Freq.bin:="Fixed_Summer"]
ff <- freq_data[, .N, by = c("Gen", "group", "run", "Freq.bin")]
gg <- freq_data[, .(Freq.bin="Fixed_Winter", N=10-.N), by = c("Gen", "group", "run")]

xx <- rbind(ff, gg)

freq_bin_run = ggplot(xx[Gen <= 50000], aes(x=Gen, y=N, col=factor(Freq.bin)))+
  geom_line()+
  facet_grid(run~group, labeller = labeller(group = group.labs))
ggsave(filename =paste0("plots/freq_bin_run.jpg"), plot = freq_bin_run , width = 15, height = 10)

ggplot(freq_data[group %in% c(6,8) & mut_freq!=1], aes(x=Gen, y=mut_freq, col=factor(mut_pos))) +
  geom_line()+
  facet_grid(run~group, labeller = labeller(group = group.labs))



times  = c(7,14,22,29, 1507, 1514, 1522, 1529, 7507, 7514, 7522, 7529, 12007, 12014, 12022, 12029, 30007, 30014, 30022, 30029, 49507, 49514, 49522, 49529)
#times  = c(8, 15, 23,10000, 10008, 10015, 10023,20000, 20008, 20015, 20023,50000,50008,50015, 50023, 80000, 80008,80015,80023, 110000, 110008, 110015, 110023)
#times  = c(6,12, 13,14,15006,15012,15013,15014,30006,30012,30013,30014,45006,45012,45013,45014,60006,60012,60013,60014,75006,75012,75013,75014,90006,90012,90013,90014,105006, 105012, 105013, 105014)
#times = c(6,12, 13, 14, 7506, 7512, 7513, 7514,22506, 22512, 22513, 22514,52506, 52512, 52513, 52514,75006, 75012, 75013, 75014,105006, 105012, 105013, 105014)
## ALL
taj_plot_all = ggplot(sum_stats[time %in% times], aes(x = chrom_pos/1000, y = tajimas_d_branch, group = block))+
  geom_vline(data = subset(sum_stats[time %in% times ]),aes(xintercept = (chromsize/2)/1000), col= "red")+
  geom_line(aes(x = chrom_pos/1000, y = tajimas_d_branch, col = group),alpha = 0.4)+
  geom_hline(yintercept = 0, col = "grey")+
  geom_line(aes(x=chrom_pos/1000, y = mean_tajd))+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance along chromosome (kb)")+
  ylab("Tajima's D")+
  theme(legend.position = "none")
ggsave(filename =paste0("plots/tajd_<20k.jpg"), plot = taj_plot_all , width = 15, height = 10)

## NON_SELECTED
taj_plot_ns = ggplot(sum_stats[time %in% times & selected==F], aes(x = chrom_pos/1000, y = tajimas_d_branch, group = block))+
  geom_line(aes(x = chrom_pos/1000, y = tajimas_d_branch, col = group),alpha = 0.4)+
  geom_hline(yintercept = 0, col = "grey")+
  geom_line(aes(x=chrom_pos/1000, y = mean_tajd))+
  facet_grid(group+selected~time, labeller = labeller(group = group.labs))+
  xlab("Distance along chromosome (kb)")+
  ylab("Tajima's D")+
  theme(legend.position = "none")

## SELECTED
taj_plot_s = ggplot(sum_stats[time %in% times & selected==T], aes(x = chrom_pos/1000, y = tajimas_d_branch, group = block))+
  geom_line(aes(x = chrom_pos/1000, y = tajimas_d_branch, col = group),alpha = 0.4)+
  geom_hline(yintercept = 0, col = "grey")+
  geom_line(aes(x=chrom_pos/1000, y = mean_tajd))+
  facet_grid(group+selected~time, labeller = labeller(group = group.labs))+
  xlab("Distance along chromosome (kb)")+
  ylab("Tajima's D")+
  theme(legend.position = "none")

dom.65 = c("1", "8", "5", "12", "6", "13", "7", "14")
dom.5 = c("1", "8", "2", "9", "3", "10", "4", "11")
## BY dominance

taj_plot_d = ggplot(sum_stats[time %in% times & group %in% dom.5], aes(x = chrom_pos/1000, y = tajimas_d_branch, group = block))+
  geom_vline(data = subset(sum_stats[time %in% times& group %in% dom.5 ]),aes(xintercept = (chromsize/2)/1000), col= "red")+
  geom_line(aes(x = chrom_pos/1000, y = tajimas_d_branch, col = group),alpha = 0.4)+
  geom_hline(yintercept = 0, col = "grey")+
  geom_line(aes(x=chrom_pos/1000, y = mean_tajd))+
  facet_grid(factor(group, levels=dom.5)~time, labeller = (group=as_labeller(group.labs)))+
  xlab("Distance along chromosome (kb)")+
  ylab("Tajima's D")+
  ggtitle("Tajima's D - d=0.5")+
  theme(legend.position = "none")


ggsave(filename =paste0("plots/tajd_relvsabsfit_dom5.jpg"), plot = taj_plot_d , width = 15, height = 10)



div_plot = ggplot(sum_stats[time %in% times], aes(x = chrom_pos/1000, y = diversity, group = block))+
  geom_vline(data = subset(sum_stats[time %in% times ]),aes(xintercept = (chromsize/2)/1000), col= "red")+
  geom_line(aes(x = chrom_pos/1000, y = diversity, col = group),alpha = 0.4)+
  geom_line(aes(x=chrom_pos/1000, y = mean_div))+
  geom_hline(yintercept = 0.00005, col = "grey", )+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  ylab("Diversity")+
  xlab("Distance along chromosome (kb)")+
  theme(legend.position = "none")

ggsave(filename =paste0("plots/div_<20k.jpg"), plot = div_plot , width = 15, height = 10)


theta_w_plot = ggplot(sum_stats[time %in% times ], aes(x = chrom_pos/1000, y = theta_w, group = block))+
  geom_vline(data = subset(sum_stats[time %in% times ],selected==T),aes(xintercept = (chromsize/2)/1000), col= "red")+
  geom_line(aes(x = chrom_pos/1000, y = theta_w, col = group),alpha = 0.4)+
  geom_line(aes(x=chrom_pos/1000, y = mean_thetaw))+
  geom_hline(yintercept = 0.0005, col = "grey", )+
  facet_grid(group+selected~time, labeller = labeller(group = group.labs))+
  xlab("Distance along chromosome (kb)")+
  ylab("Watterson's Theta")+
  theme(legend.position = "none")



ggsave(filename =paste0("plots/theta_plot_unselchrom_5.jpg"), plot = theta_w_plot , width = 15, height = 10)
ggsave(filename =paste0("plots/taj_plot_unselchrom_5.jpg"), plot = taj_plot , width = 15, height = 10)
ggsave(filename =paste0("plots/div_plot_unselchrom_5.jpg"), plot = div_plot , width = 15, height = 10)

plot = ggarrange(taj_plot_2 , div_plot , theta_w_plot , nrow = 3)
ggsave(filename =paste0("plots/summary_plot_500kb.jpg"), plot = plot, width = 25, height = 20)


h1_plot  = ggplot(sum_stats[time %in% times ], aes(x = dist/1000, y = H1, group = block))+
  geom_vline(aes(xintercept =(chromsize/2)), col= "black")+
  geom_dotplot(aes(x = dist/1000, y = H1, col = group),alpha = 0.4)+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("H1")+
  theme(legend.position = "none")


h12_plot  = ggplot(sum_stats[time %in% times  & run == 1], aes(x = dist/1000, y = H12, group = block))+
  geom_vline(aes(xintercept = 0), col= "black")+
  geom_jitter(aes(x = dist/1000, y = H12, col = group),alpha = 0.4)+
  #geom_smooth(data = sum_stats[time %in% times],aes(group = time), col = "black")+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("H1")+
  #labs(col= "Generation")
  theme(legend.position = "none")

h2h1_plot  = ggplot(sum_stats[time %in% times], aes(x = dist/1000, y = H2H1, group = block))+
  geom_vline(aes(xintercept = 0), col= "black")+
  geom_jitter(aes(x = dist/1000, y = H2H1, col = group),alpha = 0.4)+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("H1")+
  theme(legend.position = "none")





plot = ggarrange(allele_plot , taj_plot_2 , heights = c(1,2))
ggsave(filename =paste0("initial_plot_all_groups_4.jpg"), plot = plot, width = 15, height = 10)


pdf(file=paste0("~/oliviaphd/seglift_win/initial_plot_",group,".pdf"), width = 10, height = 10) 
plot = ggarrange(plot_1, plot_2, heights = c(1,2))
#plot
dev.off()
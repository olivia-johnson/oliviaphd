library(data.table)
library(ggplot2)
library(viridis)
library(stringr)
library(egg)
library(dplyr)
setwd("~/phd_data/seglift_long")

groups = c(1:8, 16:21, 43:56)
groups=c(5, 16:28)  ## equal seasons and constant population size  (CP_EG)
groups=c(1:4, 6:15) ## unequal seasons and seasonally fluctuating population size  (FP_UG)
groups=c(29:35,50:56) ## unequal season and constant population size  (CP_UG)
groups=c(36:49) ## equal seasons and seasonally fluctuating population size  (FP_EG)
windows=51
chromsize=5e5


setup="usp"

#ID = "L=10, d = 0.65, y=0.5"
  ##"L=10, d = 0.65, y=4"
  ##"No Fitness, L=10"
freq_data = NULL
sum_stats = NULL

for (g in groups){
  print(g)

  f_list  <- list.files(path =paste0("~/phd_data/seglift_long/group_", g, "/"),pattern ="al_freq_")
  
  s_list  <- list.files(path =paste0("~/phd_data/seglift_long/group_", g, "/"),pattern ="s_stat_")
  
  parameters <- fread(file=paste0("~/phd_data/seglift_long/group_",g, "/parameters.yml"), sep = ":")
  setkey(parameters, V1)
  
  fiton=parameters["fitness_on", V2]
  if (fiton==0){
    dom=0
    epi=0
  } else{
    dom = parameters["d", V2]
    epi = parameters["y", V2]
  }
  sum_gen =parameters["sum_gen", V2]
  win_gen =parameters["win_gen", V2] 
  sum_pop =parameters["s_pop", V2]
  win_pop =parameters["w_pop", V2]
  if (parameters["slim_sim", V2] == "seglift_substitute_unselchrom"){
    sim_type = "SL"  ##Segregating loci
  } else{
    sim_type = "AL"  ## All loci
  }
  ## collate al_freq files
  for (i in 1:(length(f_list ))){  
    filename = f_list [i]
    #print(filename)
    al_freq  = fread(file = paste0("~/phd_data/seglift_long/group_",g, "/",filename))
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
    s_stat  = fread(file = paste0("~/phd_data/seglift_long/group_",g, "/",filename))
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
setnames(freq_data, "Gen", "time")
freq_data [, mean_freq:=mean(mut_freq), by=c("time", "mut_pos", "group")]
#freq_data [, Season:= ifelse((Gen/2)%%15==0, "Winter", "Summer")]
freq_data[, chrom:=(mut_pos/5e5)+.5]
freq_data [, block :=paste0(group, "_", run, "_",chrom)]
freq_data[s_pop==w_pop & s_gen == w_gen, setup:="CP_EG", by = "group"] ## even season and population
freq_data[s_pop==w_pop & s_gen != w_gen, setup:="CP_UG", by = "group"] ##even population_uneven season
freq_data[s_pop!=w_pop & s_gen != w_gen, setup:="FP_UG", by="group"] ## uneven season and population
freq_data[s_pop!=w_pop & s_gen == w_gen, setup:="FP_EG", by="group"] ##uneven population_even season
freq_data[,label:=ifelse(fit==0, paste(setup, "No Fitness", sep="_"), paste(setup, d, y, sep="_")),  by="group"]
freq_data[mut_freq!=1, Freq.bin:="Segregating"]
freq_data[mut_freq==1, Freq.bin:="Fixed_Summer"]
freq_data[Freq.bin=="Segregating",n_seg:=.N, by = c("time", "label", "group","fit_type","setup")]

sum_stats[, win_pos:=n_win%%windows]
sum_stats[, midpoint:= (win_end-win_start)/2 + win_start]
sum_stats[, block :=paste0(group, "_", run, "_",chrom)]
sum_stats[, chrom_pos:= midpoint-(chrom-1)*chromsize]
sum_stats[, dist:=(chromsize/2)-chrom_pos-1, by="n_win"]
sum_stats[, selected:=ifelse(chrom%% 2==0, T, F), by="chrom"]
#sum_stats[,segregating:=F]
# sum_stats[, theta_w_corr:= theta_w*1e15]
sum_stats[s_pop==w_pop & s_gen == w_gen, setup:="CP_EG", by = "group"] ## even season and population
sum_stats[s_pop==w_pop & s_gen != w_gen, setup:="CP_UG", by = "group"] ##even population_uneven season
sum_stats[s_pop!=w_pop & s_gen != w_gen, setup:="FP_UG", by="group"] ## uneven season and population
sum_stats[s_pop!=w_pop & s_gen == w_gen, setup:="FP_EG", by="group"] ##uneven population_even season### sum_stats %>% filter(w_pop != 10000) (dplyr)
# sum_stats[, mean_thetaw:=mean(theta_w_corr), by=c("time", "win_pos", "group")]
sum_stats[,label:=ifelse(fit==0, paste("No Fitness", fit_type,setup, sep="_"), paste(d, y, fit_type, setup, sep="_")),  by="group"]

seg_data = freq_data[,.(Freq.bin), by=c("group", "run", "chrom", "time")]
sum_stat= merge(sum_stats, seg_data, by=c("group", "run", "chrom", "time"), all.x = TRUE )
sum_stat[, segregating:=ifelse(Freq.bin=="Segregating", T, F), by=c("group", "run", "chrom", "time")]
sum_stat[is.na(segregating), segregating:=F]
# sum_stat[, contributes:= ifelse(selected==T, "contribute", "neutral"), by="selected"]
# sum_stat[fit_type=="seg" & selected==T, contributes:= ifelse(Freq.bin=="Segregating", "contribute", "n_contribute"), by="Freq.bin"]
# sum_stat[is.na(contributes), contributes:="neutral"]

sum_stat[, mean_tajd:=mean(tajimas_d_branch), by=c("time", "chrom_pos", "group", "segregating")]
sum_stat[, mean_div:=mean(diversity), by=c("time", "chrom_pos", "group", "segregating")]
sum_stat[, mean_thetaw:=mean(theta_w), by=c("time", "chrom_pos", "group", "segregating")]

plot = ggplot(freq_data[group==48 & run==1 & time<10000],aes( x = time, y= mut_freq))+
  geom_line(aes( col = as.factor(mut_pos))) +
  ggtitle("Allele Frequency") +
  xlab("Generation")+
  ylab("Allele Frequency")

##ggsave(filename =("freq_plot_0_3.jpg"), plot = plot, width = 15, height = 10)
set = c("FP_UG","CP_EG")
seg_labels=unique(freq_data[Freq.bin=="Segregating" & time==50000, .(label, n_seg)])
 ##[fit_type=="AL" & d==c(0.5, 0)]

allele_plot =ggplot(data=freq_data,aes(x = time, y= mut_freq))+
  geom_line(aes(col=label, group=block),alpha = 0.3) +
  facet_wrap(~label, ncol=7)+
  #ggtitle(paste0("Allele Frequency (", ID, ")")) +
  xlab("Generation")+
  scale_y_continuous(limits = c(0, 1.0))+
  ylab("Allele Frequency")+
  theme(legend.position = "none")+
  labs(col= "Seasonal Mutation")+
  ggtitle("Allele Frequencies")

allele_plot = allele_plot + geom_text(data=seg_labels, x=45000, y=0.1, aes(label=n_seg), parse=TRUE)
ggsave(filename =paste0("new_plots/allele_freq_all.jpg"), plot = allele_plot , width = 15, height = 10)


freq_data[mut_freq!=1, Freq.bin:="Segregating"]
freq_data[mut_freq==1, Freq.bin:="Fixed_Summer"]
ff <- freq_data[, .N, by = c("time", "label", "group","fit_type","setup", "Freq.bin")]
gg <- freq_data[, .(Freq.bin="Fixed_Winter", N=80-.N), by = c("time","group","fit_type","setup", "label")]

xx <- rbind(ff, gg)

freq_bin = ggplot(xx[time<50000&fit_type=="AL"], aes(x=time, y=N, col = Freq.bin))+
  geom_line()+
  facet_wrap(~label)
ggsave(filename =paste0("plots/freq_bin_all.jpg"), plot = freq_bin , width = 15, height = 10)


freq_data[mut_freq!=1, Freq.bin:="Segregating"]
freq_data[mut_freq==1, Freq.bin:="Fixed_Summer"]
ff <- freq_data[, .N, by = c("time", "label", "group","fit_type","setup","run","Freq.bin")]
gg <- freq_data[, .(Freq.bin="Fixed_Winter", N=10-.N), by = c("time", "label", "group","fit_type","setup","run")]

xx <- rbind(ff, gg)

freq_bin_run = ggplot(xx[time<50000], aes(x=time, y=N, col=factor(Freq.bin)))+
  geom_line()+
  facet_grid(run~label)
ggsave(filename =paste0("plots/freq_bin_run_all.jpg"), plot = freq_bin_run , width = 15, height = 10)

ggplot(freq_data[group %in% c(6,8) & mut_freq!=1], aes(x=Gen, y=mut_freq, col=factor(mut_pos))) +
  geom_line()+
  facet_grid(run~group, labeller = labeller(group = group.labs))



#times  = c(8, 15, 23,50000,50008,50015, 50023, 110000, 110008, 110015, 110023)
#times  = c(8, 15, 23,10000, 10008, 10015, 10023,20000, 20008, 20015, 20023,50000,50008,50015, 50023, 80000, 80008,80015,80023, 110000, 110008, 110015, 110023)
#times  = c(6,12, 13,14,15006,15012,15013,15014,30006,30012,30013,30014,45006,45012,45013,45014,60006,60012,60013,60014,75006,75012,75013,75014,90006,90012,90013,90014,105006, 105012, 105013, 105014)
times = c(6,12,13,14,7506,7512,7513,7514,15006 ,15012,15013,15014,22506,22512,22513,22514,30006,30012,30013,30014,37506,37512,37513,37514,45006,45012,45013,45014)

#symmetrial times
times= c(7, 14,22,29, 9007, 9014, 9022, 9029, 18007, 18014, 18022, 18029, 27007, 27014, 27022, 27029,36007, 36014, 36022, 3602,45007, 45014, 45022, 45029)
## ALL

plotdata= sum_stat[time %in% times]


taj_plot_all = ggplot(plotdata, aes(x = chrom_pos/1000, y = tajimas_d_branch, group = block))+
  geom_vline(data = subset(plotdata,segregating==T),aes(xintercept = (chromsize/2)/1000), col= "red")+
  geom_line(aes(x = chrom_pos/1000, y = tajimas_d_branch, col = group),alpha = 0.4)+
  geom_hline(yintercept = 0, col = "grey")+
  geom_line(aes(x=chrom_pos/1000, y = mean_tajd))+
  facet_grid(label+segregating~time)+
  xlab("Distance along chromosome (kb)")+
  ylab("Tajima's D")+
  theme(legend.position = "none")


tajd_plot=ggplot(plotdata[y==0.5], aes(x = chrom_pos/1000, y = tajimas_d_branch, group = block))+
  geom_vline(data = subset(plotdata[y==0.5 & segregating==T],segregating==T),aes(xintercept = (chromsize/2)/1000), col= "red")+
  geom_line(aes(x = chrom_pos/1000, y = tajimas_d_branch, col = label),alpha = 0.4)+
  geom_line(aes(x=chrom_pos/1000, y = mean_tajd))+
  geom_hline(yintercept = 0, col = "grey", )+
  facet_grid(label+segregating~time)+
  ylab("Tajima's D")+
  xlab("Distance along chromosome (kb)")+
  ggtitle("uneven seasons and changing population size- y=0.5")+
  theme(legend.position = "none")

ggsave(filename =paste0("new_plots/taj_",setup,"_d5.jpg"), plot = tajd_plot , width = 20, height = 10)

times = c(6,12,13,14,3006,  3012,  3013,  3014,6006,  6012,  6013,  6014, 9006,  9012,  9013,  9014, 12006 ,12012, 12013, 12014, 15006, 15012, 15013, 15014,18006, 18012, 18013 ,18014, 21006, 21012, 21013 ,21014,24006, 24012, 24013, 24014)

div_plot = ggplot(plotdata[y==0.5], aes(x = chrom_pos/1000, y = diversity, group = block))+
  geom_vline(data = subset(plotdata[y==0.5 & segregating==T],segregating==T),aes(xintercept = (chromsize/2)/1000), col= "red")+
  geom_line(aes(x = chrom_pos/1000, y = diversity, col = label),alpha = 0.4)+
  #geom_hline(yintercept = 0.00012, col = "grey", )+ ##ep
  geom_hline(yintercept = 0.00005, col = "grey", )+ ##usp
  #geom_hline(yintercept = 0.000021, col = "grey", )+ ##es_up
  geom_line(aes(x=chrom_pos/1000, y = mean_div))+ 
  facet_grid(label+segregating~time)+
  ylab("Diversity")+
  xlab("Distance along chromosome (kb)")+
  ggtitle("uneven seasons and changing population size - y=0.5")+
  theme(legend.position = "none")

ggsave(filename =paste0("new_plots/div_",setup,"_d5.jpg"), plot = div_plot , width = 20, height = 10)

times = c(7506,7512,7513,7514,15006 ,15012,15013,15014,22506,22512,22513,22514,30006,30012,30013,30014,37506,37512,37513,37514,45006,45012,45013,45014)

theta_w_plot = ggplot(sum_stats[time %in% times& d==0.65], aes(x = chrom_pos/1000, y = theta_w, group = block))+
  geom_vline(data = subset(sum_stats[time %in% times& d==0.65 ],selected==T),aes(xintercept = (chromsize/2)/1000), col= "red")+
  geom_line(aes(x = chrom_pos/1000, y = theta_w, col = group),alpha = 0.4)+
  geom_line(aes(x=chrom_pos/1000, y = mean_thetaw))+
  geom_hline(yintercept = 0.0005, col = "grey", )+
  facet_grid(label+selected~time)+
  xlab("Distance along chromosome (kb)")+
  ylab("Watterson's Theta")+
  theme(legend.position = "none")



ggsave(filename =paste0("plots/theta_plot_unselchrom_5.jpg"), plot = theta_w_plot , width = 15, height = 10)
ggsave(filename =paste0("plots/taj_plot_unselchrom_5.jpg"), plot = taj_plot , width = 15, height = 10)
ggsave(filename =paste0("plots/div_plot_unselchrom_5.jpg"), plot = div_plot , width = 15, height = 10)

plot = ggarrange(taj_plot_2 , div_plot , theta_w_plot , nrow = 3)
ggsave(filename =paste0("plots/summary_plot_500kb.jpg"), plot = plot, width = 25, height = 20)


h1_plot  = ggplot(sum_stats[time %in% times& d==0.65], aes(x = dist/1000, y = H1, group = block))+
  geom_vline(aes(xintercept =(chromsize/2)), col= "black")+
  geom_dotplot(aes(x = dist/1000, y = H1, col = group),alpha = 0.4)+
  facet_grid(label~time,)+
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


pdf(file=paste0("~/phd_data/seglift_win/initial_plot_",group,".pdf"), width = 10, height = 10) 
plot = ggarrange(plot_1, plot_2, heights = c(1,2))
#plot
dev.off()

library(data.table)
library(ggplot2)
library(viridis)
library(stringr)
library(egg)
library(dplyr)
library(ggpubr)

setwd("~/phd_data/seglift_linked_allele")

groups = c(1:9)


alfreq_data = NULL

for (g in groups){
  print(g)
  
  f_list  <- list.files(path =paste0("~/phd_data/seglift_linked_allele/group_", g, "/"),pattern ="al_freq_")
  
  parameters <- fread(file=paste0("~/phd_data/seglift_linked_allele/group_",g, "/parameters.yml"), sep = ":")
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
  loci =parameters["l", V2]

  ## collate al_freq files
  for (i in 1:(length(f_list ))){  
    filename = f_list [i]
    #print(filename)
    al_freq  = fread(file = paste0("~/phd_data/seglift_linked_allele/group_",g, "/",filename))
    al_freq [,run := i]
    al_freq [,group:=g]
    al_freq [,d:=dom]
    al_freq [,y:=epi]
    al_freq [,l:=loci]
    al_freq [,fit:=fiton]
    al_freq [,s_gen:=sum_gen]
    al_freq [,w_gen:=win_gen]
    al_freq [,s_pop:=sum_pop]
    al_freq [,w_pop:=win_pop]
    alfreq_data  = rbind(alfreq_data , al_freq )
  }
}
setnames(alfreq_data, "Gen", "time")
alfreq_data [, mean_freq:=mean(mut_freq), by=c("time", "mut_pos", "group")]
#alfreq_data [, Season:= ifelse((Gen/2)%%15==0, "Winter", "Summer")]
alfreq_data[, chrom:=(mut_pos/5e5)+.5]
alfreq_data [, block :=paste0(group, "_", run, "_",chrom)]
alfreq_data[s_pop==w_pop & s_gen == w_gen, setup:="CP_EG", by = "group"] ## even season and population
alfreq_data[s_pop==w_pop & s_gen != w_gen, setup:="CP_UG", by = "group"] ##even population_uneven season
alfreq_data[s_pop!=w_pop & s_gen != w_gen, setup:="FP_UG", by="group"] ## uneven season and population
alfreq_data[s_pop!=w_pop & s_gen == w_gen, setup:="FP_EG", by="group"] ##uneven population_even season
alfreq_data[,label:=ifelse(fit==0, paste("No Fitness",l,setup, sep="_"), paste(l, d, y,setup, sep="_")),  by="group"]
alfreq_data[mut_freq!=1, Freq.bin:="Segregating"]
alfreq_data[mut_freq==1, Freq.bin:="Fixed_Summer"]
alfreq_data[Freq.bin=="Segregating",n_seg:=.N, by = c("time", "label", "group","setup")]

seg_labels=unique(alfreq_data[Freq.bin=="Segregating" & time==50000, .(label, n_seg)])

allele_plot =ggplot(data=alfreq_data,aes(x = time, y= mut_freq))+
  geom_line(aes(group=block, col = group),alpha = 0.05) +
  facet_wrap(~label)+
  #ggtitle(paste0("Allele Frequency (", ID, ")")) +
  xlab("Generation")+
  scale_y_continuous(limits = c(0, 1.0))+
  ylab("Allele Frequency")+
  theme(legend.position = "none")+
  labs(col= "Seasonal Mutation")
#ggtitle("Allele Frequencies")

allele_plot = allele_plot + geom_text(data=seg_labels, x=15000, y=0.1, aes(label=n_seg), parse=TRUE)

ggexport(allele_plot, filename="plots/allele_plot_linked.pdf")

alfreq_data[mut_freq!=1, Freq.bin:="Segregating"]
alfreq_data[mut_freq==1, Freq.bin:="Fixed_Summer"]
ff <- alfreq_data[, .N, by = c("time", "label", "group","setup", "Freq.bin")]
gg <- alfreq_data[, .(Freq.bin="Fixed_Winter", N=(8*l)-.N), by = c("time","group","setup", "label")]

xx <- rbind(ff, gg)

freq_bin = ggplot(xx, aes(x=time, y=N, col = Freq.bin))+
  geom_line()+
  facet_wrap(~label)
ggexport(freq_bin, filename="plots/freq_bin.pdf")



alfreq_data[mut_freq!=1, Freq.bin:="Segregating"]
alfreq_data[mut_freq==1, Freq.bin:="Fixed_Summer"]
ff <- alfreq_data[, .N, by = c("time", "label", "group","setup","run","Freq.bin")]
gg <- alfreq_data[, .(Freq.bin="Fixed_Winter", N=l-.N), by = c("time", "label", "group","setup","run")]

zz <- rbind(ff, gg)

freq_bin_run = ggplot(zz, aes(x=time, y=N, col=factor(Freq.bin)))+
  geom_line()+
  facet_grid(run~label)
ggexport(freq_bin_run, filename="plots/freq_bin_run.pdf")




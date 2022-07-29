library(data.table)
library(ggplot2)
library(viridis)
library(stringr)
library(egg)
library(dplyr)
library(ggpubr)
library(psych)
setwd("~/phd_data/Results/single_locus/")

freq_data = NULL
sum_stats = NULL

for (g in groups){
  print(g)
  
  f_list  <- list.files(path =paste0("~/phd_data/Results/single_locus/group_", g, "/"),pattern ="al_freq_")

  s_list  <- list.files(path =paste0("~/phd_data/Results/single_locus/group_", g, "/"),pattern ="sim_stat_")

  parameters <- fread(file=paste0("~/phd_data/hpc_parameters/single_locus/group_",g, ".txt"), sep = ":")
  setkey(parameters, V1)
  
  fiton=parameters["fitness_on", V2]
  if (fiton==0){
    h_s = 0
    h_w = 0
    s_s = 0
    s_w = 0
  } else{
    h_s = parameters["h_s", V2]
    h_w = parameters["h_w", V2]
    s_s = parameters["s_s", V2]
    s_w = parameters["s_w", V2]
  }
  sum_gen =parameters["sum_gen", V2]
  win_gen =parameters["win_gen", V2] 
  if (sum_gen==win_gen){
    gen_s = "EG"
  } else{
    gen_s = "UG"
  }
  sum_pop =parameters["s_pop", V2]
  win_pop =parameters["w_pop", V2]
  if (sum_pop==win_pop){
    pop_s = "CP"
  } else{
    pop_s = "FP"
  }
  burnin = parameters["burnin_Ne", V2]
  mr = parameters["mutRate", V2]
  rr = parameters["recRate", V2]
  nWin = parameters["winpChrom", V2]
  genomeSize = parameters["genomeSize", V2]
  
  g_label = ifelse(fiton==0, paste( pop_s, gen_s, "No Fitness", sep="_"), paste(pop_s, "_",gen_s,"_", "h", h_s,"s", s_s, sep=":"))
  ## collate al_freq files
  for (i in 1:(length(f_list ))){  
    filename = f_list [i]
    #print(filename)
    al_freq  = fread(file = paste0("~/phd_data/seglift_long/group_",g, "/",filename))
    # al_freq  = fread(file = paste0("~/phd_data/Results/background/group_",g, "/",filename))
    al_freq[,run:=i]
    al_freq[,label:=g_label]
    al_freq[,group:=g]
    al_freq [,h_s:=h_s]
    al_freq [,h_w:=h_w]
    al_freq [,s_s:=s_s]
    al_freq [,s_w:=s_w]
    al_freq [,fit:=fiton]
    al_freq [,s_gen:=sum_gen]
    al_freq [,w_gen:=win_gen]
    al_freq [,s_pop:=sum_pop]
    al_freq [,w_pop:=win_pop]
    al_freq[, pop_season := pop_s]
    al_freq[, gen_season := gen_s]
    freq_data  = rbind(freq_data , al_freq )
  }
  
  ## collate sum_stat files
  for (i in 1:(length(s_list ))){  
    filename = s_list [i]
    #print(filename)
    # s_stat  = fread(file = paste0("~/phd_data/Results/single_locus/group_",g, "/",filename))
    s_stat[,run:=i]
    s_stat[,label:=g_label]
    s_stat[,group:=g]
    s_stat [,h_s:=h_s]
    s_stat [,h_w:=h_w]
    s_stat [,s_s:=s_s]
    s_stat [,s_w:=s_w]
    s_stat [,fit:=fiton]
    s_stat [,s_gen:=sum_gen]
    s_stat [,w_gen:=win_gen]
    s_stat [,s_pop:=sum_pop]
    s_stat [,w_pop:=win_pop]
    s_stat[, pop_season := pop_s]
    s_stat[, gen_season := gen_s]
    s_stat[, burninNe:= burnin]
    s_stat[, mutRate:=mr]
    s_stat[, recRate:=rr]
    s_stat[, genomeSize := genomeSize]
    sum_stats = rbind(sum_stats, s_stat,fill=TRUE)
  }
}

freq_data [, mean_freq:=mean(mut_freq), by=c("Gen", "mut_pos", "group")]
freq_data [, block :=paste0(group, "_", run)]


sum_stats[, win_pos:=nWin%%windows]
sum_stats[, midpoint:= (win_end-win_start)/2 + win_start]
sum_stats[, block :=paste0(group, "_", run)]
sum_stats[, exp:=4*as.numeric(burninNe)*as.numeric(mutRate), by="group"]
setnames(sum_stats, "time", "Gen")
sum_stats[gen_season == "EG", gen_year:=Gen%%20, by=c("label", "Gen")]
sum_stats[gen_season == "UG", gen_year:=Gen%%12, by=c("label", "Gen")] 
sum_stats[gen_season == "EG", season:=ifelse((gen_year<11 & gen_year>0), "summer", "winter"), by=c("Gen")]
sum_stats[gen_season == "UG", season:=ifelse((gen_year<11 & gen_year>0), "summer", "winter"), by=c("Gen")]

sum_stats[, mean_tajd:=mean(tajimas_d_branch), by=c("Gen", "group")]
sum_stats[, mean_div:=mean(diversity), by=c("Gen", "group")]
sum_stats[, mean_thetaw:=mean(theta_w), by=c("Gen", "group")]

div_var<- sum_stat[selected==F, mean(diversity), by=c("label", "Gen", "gen_season", "pop_season", "burninNe", "mutRate")]
setnames(div_var, "V1", "gw")
div_var[gen_season == "EG", gen_year:=Gen%%20, by=c("label", "Gen")]
div_var[gen_season == "UG", gen_year:=Gen%%12, by=c("label", "Gen")] 
div_var[gen_season == "EG", season:=ifelse((gen_year<11 & gen_year>0), "summer", "winter"), by=c("Gen")]
div_var[gen_season == "UG", season:=ifelse((gen_year<11 & gen_year>0), "summer", "winter"), by=c("Gen")]
div_var[, exp:=4*as.numeric(burninNe)*as.numeric(mutRate), by="label"]
div_var[, diff_div:=gw-exp, by=c("Gen", "label", "gen_season", "pop_season")]
div_var[, div_reldiff:=gw/exp, by=c("Gen", "label", "gen_season", "pop_season")]
div_var[,reldiff_mean:=harmonic.mean(div_reldiff), by=c("label")]
div_var[,gw_hmean:=harmonic.mean(gw), by=c("label")]
div_var[,div_ne:=gw/(4*as.numeric(mutRate))]
div_var[,exp_ne:=exp/(4*as.numeric(mutRate)), by="label"]
div_var[, d_ne_diff:=div_ne/exp_ne]



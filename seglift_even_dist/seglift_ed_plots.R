library(data.table)
library(ggplot2)
library(viridis)
library(stringr)

setwd("~/oliviaphd/seglift_even_dist/")

group = 13

## list of s_stat files
f_list <- list.files(path = "~/oliviaphd/seglift_even_dist/py_out/", pattern=paste0("sim_s_stat_", group, "_"))
hap_list <- list.files(path = "~/oliviaphd/seglift_even_dist/slim_out/", pattern=paste0("ms_out_", group, "_", run, "_"))

## merge output files
s_stat = NULL
haplotypes = NULL
hap_stats = NULL
freq_data = NULL


###  For 0 selected mutations
for (i in 1:(length(f_list))){
  file_name = f_list[i]
  output = fread(file = paste0("~/oliviaphd/seglift_even_dist/py_out/",file_name))
  output[,run := i]
  s_stat = rbind(s_stat, output)
  
  # for (j in (unique(output$time))){
  #   hapout = fread(file = paste0("~/oliviaphd/seglift_even_dist/slim_out/","ms_out_", group, "_", run, "_", j, ".txt"))
  #   hapout = hapout[-(1:2),.N, by="//"]
  #   setnames(hapout, "//", "hap")
  #   setnames(hapout, "N", "hap_count")
  #   hapout[, n_sel:= str_count(hap, "1")]
  #   haps = hapout[, .N, by=n_sel]
  #   haps$run = run
  #   haps$gen = j
  #   hapout$run = run
  #   hapout$gen = j
  #   haplotypes = rbind(haplotypes, hapout)
  #   hap_stats = rbind(haps, hap_stats)
  # }
}
##For >0 selected mutations

for (i in 1:(length(f_list))){
  file_name = f_list[i]
  output = fread(file = paste0("~/oliviaphd/seglift_even_dist/py_out/",file_name))
  al_freq = fread(file = paste0("~/oliviaphd/seglift_even_dist/slim_out/al_freq_", group, "_", i, ".txt"))
  output[,run := i]
  al_freq[,run := i]
  s_stat = rbind(s_stat, output)
  freq_data = rbind(freq_data, al_freq)
  
  for (j in (unique(output$time))){
    hapout = fread(file = paste0("~/oliviaphd/seglift_even_dist/slim_out/","ms_out_", group, "_", run, "_", j, ".txt"))
    hapout = hapout[-(1:2),.N, by="//"]
    setnames(hapout, "//", "hap")
    setnames(hapout, "N", "hap_count")
    hapout[, n_sel:= str_count(hap, "1")]
    haps = hapout[, .N, by=n_sel]
    haps$run = run
    haps$gen = j
    hapout$run = run
    hapout$gen = j
    haplotypes = rbind(haplotypes, hapout)
    hap_stats = rbind(haps, hap_stats)
  }
}

s_stat[, midpoint := (s_stat$win_start + (s_stat$win_end - s_stat$win_start)/2)]

s_stat[, branch_mean := mean(tajimas_d_branch), by=c("n_win", "time")]
s_stat[, div_mean := mean(diversity), by=c("n_win", "time")]



## PLOTS ##
# ggplot(freq_data, aes(Gen, mut_freq))+
#   geom_point()+
#   facet_wrap(~factor(mut_pos))
# 
# 
# pdf(file=paste0("~/oliviaphd/seglift_even_dist/plots/t_vs_tajd_branch_", group,".pdf"), width = 10, height = 5)
# ggplot(s_stat, aes(x = time, y = tajd_branch_1)) +
#   geom_area()+
#   scale_color_viridis()+
#   xlab("Time (generations)")+
#   ylab("Tajima's D (mean)")+
#   ggtitle("t_vs_tajd_branch_1")
# dev.off()

pdf(file=paste0("~/oliviaphd/seglift_even_dist/plots/win_vs_tajd_branch_mean_", group,".pdf"), width = 10, height = 10)
ggplot(s_stat[time != 10000 ], aes(x = midpoint/1000, y = branch_mean)) +
  geom_line()+
  facet_wrap(~time, nrow=5)+
  xlab("Position (Kb)")+
  ylab("Tajima's D (mean)")+
  ggtitle("win_vs_tajd_branch_mean")+
  geom_vline(xintercept=unique(freq_data$mut_pos)/1000,color="red")
dev.off()


# pdf(file=paste0("~/oliviaphd/seglift_even_dist/plots/t_vs_tajd_branch_", group,".pdf"), width = 10, height = 5)
# ggplot(s_stat, aes(x = time, y = branch_mean)) +
#   geom_line()+
#   facet_wrap(~n_win)+
#   xlab("Time (generations)")+
#   ylab("Tajima's D (mean)") +
#   ggtitle("t_vs_tajd_branch_mean")
# dev.off()


pdf(file=paste0("~/oliviaphd/seglift_even_dist/plots/win_vs_div_", group,".pdf"), width = 10, height = 10) 
ggplot(data=s_stat, aes(midpoint/1000, div_mean))+
  geom_line()+
  facet_wrap(~time)+
  xlab("Position (Kb)")+
  ylab("Diversity (pi)") +
  geom_vline(xintercept=unique(freq_data$mut_pos)/1000,color="red", size=0.25 )+
  ggtitle("win_vs_diversity")
dev.off()


pdf(file=paste0("~/oliviaphd/seglift_even_dist/plots/win_vs_tajd_branchs_", group,".pdf"), width = 10, height = 10) 
ggplot(data=s_stat, aes(midpoint/1000, tajimas_d_branch, col = factor(run)))+
  geom_line()+
  facet_wrap(~time)+
  xlab("Position (Kb)")+
  ylab("Tajima's D") +
  #geom_vline(xintercept=unique(freq_data$mut_pos)/1000,color="red" )+
  ggtitle("win_vs_tajd_branch")
dev.off()
library(data.table)
library(ggplot2)
library(viridis)

setwd("~/oliviaphd/seglift_even_dist/")

group = 2

## list of s_stat files
f_list <- list.files(path = "~/oliviaphd/seglift_even_dist/py_out/", pattern=paste0("sim_s_stat_", group, "_"))

##sim_data = fread(file = ~/oliviaphd/seglift_even_dist/slim_out/sim_data )
## merge output files
s_stat = NULL
for (i in 1:(length(f_list))){
  file_name = f_list[i]
  output = fread(file = paste0("~/oliviaphd/seglift_even_dist/py_out/",file_name))
  output[,run := i]
  s_stat = rbind(s_stat, output)
}

s_stat[, midpoint := (s_stat$win_start + (s_stat$win_end - s_stat$win_start)/2)]

s_stat[, branch_mean := mean(tajimas_d_branch), by=c("n_win", "time")]



## PLOTS ##

pdf(file=paste0("~/oliviaphd/seglift_even_dist/plots/t_vs_tajd_branch_", group,".pdf"), width = 10, height = 5)
ggplot(s_stat, aes(x = time, y = tajd_branch_1)) +
  geom_area()+
  scale_color_viridis()+
  xlab("Time (generations)")+
  ylab("Tajima's D (mean)")+
  ggtitle("t_vs_tajd_branch_1")
dev.off()

pdf(file=paste0("~/oliviaphd/seglift_even_dist/plots/win_vs_tajd_branch_mean_", group,".pdf"), width = 10, height = 5)
ggplot(s_stat[time != 10000], aes(x = midpoint/1000, y = branch_mean)) +
  geom_line()+
  facet_wrap(~time, nrow=5)+
  xlab("Position (Kb)")+
  ylab("Tajima's D (mean)")+
  ggtitle("win_vs_tajd_branch_mean")
dev.off()


pdf(file=paste0("~/oliviaphd/seglift_even_dist/plots/t_vs_tajd_branch_", group,".pdf"), width = 10, height = 5)
ggplot(s_stat, aes(x = time, y = branch_mean)) +
  geom_line()+
  facet_wrap(~n_win)+
  xlab("Time (generations)")+
  ylab("Tajima's D (mean)") +
  ggtitle("t_vs_tajd_branch_mean")
dev.off()


#ggplot(stat_tajd, aes(x = factor(time), y = branch_mean)) +
geom_boxplot()


pdf(file=paste0("~/oliviaphd/seglift_even_dist/plots/win_vs_tajd_branch_", group,".pdf"), width = 10, height = 5) 
ggplot(data=s_stat, aes(midpoint/1000, tajimas_d_branch, col = factor(run)))+
  geom_line()+
  facet_wrap(~time)+
  xlab("Position (Kb)")+
  ylab("Tajima's D") +
  ggtitle("win_vs_tajd_branch")
dev.off()
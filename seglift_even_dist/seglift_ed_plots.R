library(data.table)
library(ggplot2)
library(viridis)

setwd("~/oliviaphd/seglift_treeseq/")

## set number of runs -- call in from python?
group = 1

## list of s_stat files
f_list <- list.files(path = "~/oliviaphd/seglift_treeseq/py_out/", pattern=paste0("sim_s_stat_", group, "_"))

## merge output files
s_stat = NULL
for (i in 1:(length(f_list))){
  file_name = f_list[i+1]
  output = fread(file = paste0("~/oliviaphd/seglift_treeseq/py_out/",file_name))
  output[,run = i]
  rbind(s_stat, output)
}

stat_tajd[, midpoint := (stat_tajd$win_start + (stat_tajd$win_end - stat_tajd$win_start)/2)]


## PLOTS ##

pdf(file="~/oliviaphd/seglift_treeseq/plots/t_vs_tajd_branch_1.pdf", width = 10, height = 5)
ggplot(stat_tajd, aes(x = time, y = tajd_branch_1)) +
  geom_area()+
  scale_color_viridis()+
  xlab("Time (generations)")+
  ylab("Tajima's D (mean)")+
  ggtitle("t_vs_tajd_branch_1")
dev.off()

pdf(file="~/oliviaphd/seglift_treeseq/plots/win_vs_tajd_branch_mean.pdf", width = 10, height = 5)
ggplot(stat_tajd, aes(x = win_mid, y = branch_mean)) +
  geom_area()+
  facet_wrap(~time)+
  xlab("Position (bp)")+
  ylab("Tajima's D (mean)")+
  ggtitle("win_vs_tajd_branch_mean")
dev.off()


pdf(file="~/oliviaphd/seglift_treeseq/plots/t_vs_tajd_branch_mean.pdf", width = 10, height = 5)
ggplot(stat_tajd, aes(x = time, y = branch_mean)) +
  geom_area()+
  facet_wrap(~n_win)+
  xlab("Time (generations)")+
  ylab("Tajima's D (mean)") +
  ggtitle("t_vs_tajd_branch_mean")
dev.off()


#ggplot(stat_tajd, aes(x = factor(time), y = branch_mean)) +
  geom_boxplot()

  
ggplot(data=output, aes(n_win, tajimas_d_branch))+
  geom_area()+
  facet_wrap(~time)

ggplot(data=stat_tajd, aes(n_win, tajimas_d_branch))+
  geom_area()+
  facet_wrap(~time)

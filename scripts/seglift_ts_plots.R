library(data.table)
library(ggplot2)
library(viridis)

setwd("~/oliviaphd/seglift_treeseq/")

## set number of runs -- call in from python?
runs = 10

## list of s_stat files
f_list <- list.files(path = "~/oliviaphd/seglift_treeseq/py_out/", pattern="sim_s_stat_")

## cycle throu 
for (i in 0:(runs-1)){
  file_name = f_list[i+1]
  if (file_name == "sim_s_stat_0"){
    stat_tajd = fread(file = paste0("~/oliviaphd/seglift_treeseq/py_out/",file_name), select = 2:5)
    stat_div = fread(file = paste0("~/oliviaphd/seglift_treeseq/py_out/",file_name), select = 2:5)
    }
  else{
    output = fread(file = paste0("~/oliviaphd/seglift_treeseq/py_out/",file_name))
    stat_tajd[, paste("sel_mut", i, sep="_"):=output$n_s_mut]
    stat_tajd[, paste("neut_mut", i, sep="_"):=output$n_n_mut]
    stat_tajd[, paste("tajd_branch", i, sep="_"):=output$tajimas_d_branch]
    stat_tajd[, paste("tajd_site", i, sep="_"):=output$tajimas_d_site]
    stat_div[, paste("sel_mut", i, sep="_"):=output$n_s_mut]
    stat_div[, paste("neut_mut", i, sep="_"):=output$n_n_mut]
    stat_div[, paste("div", i, sep="_"):=output$diversity]
  }
}

stat_tajd[, branch_mean := rowMeans(.SD), .SDcols = patterns("tajd_branch"),by = c("time", "n_win")]
stat_tajd[, site_mean := rowMeans(.SD), .SDcols = patterns("tajd_site"),by = c("time", "n_win")]
stat_tajd[, win_mid := (stat_tajd$win_start + (stat_tajd$win_end - stat_tajd$win_start)/2)]
#stat_tajd[, branch_min := pmin(.SD, na.rm=TRUE ), .SDcols = branch_idx,by = c("time", "n_win")]
#stat_tajd[, branch_max := pmax(.SD, na.rm=TRUE ), .SDcols = branch_idx,by = c("time", "n_win")]


## PLOTS ##

pdf(file="~/oliviaphd/seglift_treeseq/plots/t_vs_tajd_branch_1.pdf", width = 10, height = 5)
ggplot(stat_tajd, aes(x = time, y = tajd_branch_1, col = win_mid)) +
  geom_point()+
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
ggplot(stat_tajd, aes(x = time, y = branch_mean, col = n_win)) +
  geom_point()+
  scale_color_viridis()+
  xlab("Time (generations)")+
  ylab("Tajima's D (mean)") +
  ggtitle("t_vs_tajd_branch_mean")
dev.off()


#ggplot(stat_tajd, aes(x = factor(time), y = branch_mean)) +
  geom_boxplot()

  
ggplot(data=output, aes(n_win, tajimas_d_branch))+
  geom_area()+
  facet_wrap(~time)

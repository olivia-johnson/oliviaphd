library(data.table)
library(ggplot2)
library(viridis)
library(stringr)
library(egg)
setwd("~/oliviaphd/seglift_win")

groups = c("1", "6", "7")

#ID = "L=10, d = 0.65, y=0.5"
  ##"L=10, d = 0.65, y=4"
  ##"No Fitness, L=10"
freq_data = NULL
sum_stats = NULL

for (g in groups){
  print(g)

  f_list <- list.files(path =paste0("~/oliviaphd/seglift_win/group_", g, "/"),pattern ="al_freq_")
  
  s_list <- list.files(path =paste0("~/oliviaphd/seglift_win/group_", g, "/"),pattern ="s_stat_2sps")
  
  
  ## collate al_freq files
  for (i in 1:(length(f_list))){  
    filename = f_list[i]
    #print(filename)
    al_freq = fread(file = paste0("~/oliviaphd/seglift_win/group_",g, "/",filename))
    al_freq[,run := i]
    al_freq[,group:=g]
    freq_data = rbind(freq_data, al_freq)
  }
  
  ## collate sum_stat files
  for (i in 1:(length(s_list))){  
    filename = s_list[i]
    #print(filename)
    s_stat = fread(file = paste0("~/oliviaphd/seglift_win/group_",g, "/",filename))
    s_stat[,run := i]
    s_stat[,group:=g]
    sum_stats = rbind(sum_stats, s_stat)
  }
}

freq_data[, mean_freq:=mean(mut_freq), by=c("Gen", "mut_pos", "group")]
#freq_data[, Season:= ifelse((Gen/2)%%15==0, "Winter", "Summer")]
freq_data[, block :=paste0(group, "_", run, "_",(mut_pos/1e5)+.5)]


sum_stats[, win_pos:=n_win%%21]
sum_stats[, mean_tajd:=mean(tajimas_d_branch), by=c("time", "win_pos", "group")]
sum_stats[, mean_div:=mean(diversity), by=c("time", "win_pos", "group")]
sum_stats[, midpoint:= (win_end-win_start)/2 + win_start]
sum_stats[, block :=paste0(group, "_", run, "_",(n_win%/%21))]
sum_stats[, dist:=49999-(midpoint-(n_win%/%21)*1e5), by="n_win"]
sum_stats[, theta_w_corr:= theta_w*1e15]
setnames(sum_stats, "H2/H1", "H2H1")



plot = ggplot(freq_data[run==4],aes( x = Gen, y= mut_freq))+
  geom_line(aes( col = as.factor(mut_pos))) +
  ggtitle("Allele Frequency (Single Simulation, N = 1000, d=0.65)") +
  xlab("Generation")+
  ylab("Allele Frequency")

##ggsave(filename =("freq_plot_0_3.jpg"), plot = plot, width = 15, height = 10)

group.labs <- c("1" = "No Fitness", "6"="y=4, d=0.65","7"= "y=0.5, d=0.65")


allele_plot=ggplot(freq_data,aes( x = Gen, y= mut_freq, group = block))+
  geom_line(aes(col=as.factor(group)),alpha = 0.3) +
  facet_wrap(~group, labeller = labeller(group = group.labs))+
  #ggtitle(paste0("Allele Frequency (", ID, ")")) +
  xlab("Generation")+
  scale_y_continuous(limits = c(0, 1.0))+
  ylab("Allele Frequency")+
  theme(legend.position = "none")+
  labs(col= "Seasonal Mutation")
# 
# plot_meanfreq = ggplot(freq_data,aes( x = Gen, y= mean_freq))+
#  geom_line(alpha = 0.3) +
#   facet_wrap(~group)+
#   ggtitle(paste0("Mean Allele Frequency (", ID, ")")) +
#   xlab("Generation")+
#   scale_y_continuous(limits = c(0, 1.0))+
#   ylab("Mean Allele Frequency")+
#   labs(col= "Seasonal Mutation")

times = c(8, 15, 23,600, 608, 615, 623, 2100, 2108, 2115, 2123,4200, 4208, 4215, 4223, 6300, 6308, 6315, 6323,8100, 8108, 8115, 8123, 9900,9908, 9915, 9923)
# plot_2=ggplot(sum_stats[time %in% times], aes(x = win_end/1000, y = tajimas_d_branch, group=run)) +
#   geom_line(alpha = 0.5)+
#   #geom_line(aes(col = as.factor(run)), alpha = 0.8)+
#   geom_line(aes(y=mean_tajd), col="blue")+
#   facet_wrap(~time, nrow = 2)+
#   xlab("Position (Kb)")+
#   ylab("Tajima's D")+
#   ggtitle(paste0("Mean Tajima's D (", ID, ")"))+
#   geom_vline(xintercept=unique(freq_data$mut_pos)/1000,color="red", alpha = 0.5)+
#   theme(legend.position = "none")


# plot_meantaj=ggplot(sum_stats[time %in% times], aes(x = win_end/1000, y = mean_tajd)) +
#   geom_line()+
#   facet_wrap(~time, nrow = 3)+
#   xlab("Position (Kb)")+
#   ylab("Tajima's D (mean)")+
#   ggtitle(paste0("Mean Tajima's D (", ID, ")"))+
#   geom_vline(xintercept=unique(freq_data$mut_pos)/1000,color="red", alpha = 0.5)+
#   theme(legend.position = "none")



# taj_plot_1 = ggplot(sum_stats[time %in% times], aes(x = dist/1000, y = tajimas_d_branch, group = block))+
#   geom_line(alpha = 0.3)+
#   facet_wrap(~time, nrow = 2)+
#   xlab("Distance from selected loci (kb)")+
#   ylab("Tajima's D")+
#   ggtitle(paste0("Tajima's D (", ID, ")"))

taj_plot_2 = ggplot(sum_stats[time %in% times], aes(x = dist/1000, y = tajimas_d_branch, group = block))+
  geom_vline(aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist/1000, y = tajimas_d_branch, col = group),alpha = 0.4)+
  geom_line(aes(x=dist/1000, y = mean_tajd))+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("Tajima's D")+
  #labs(col= "Generation")
  theme(legend.position = "none")
  #ggtitle(paste0("Tajima's D (", ID, ")"))
ggsave(filename =paste0("plots/taj_plot_2sps.jpg"), plot = taj_plot_2, width = 25, height = 10)


div_plot = ggplot(sum_stats[time %in% times], aes(x = dist/1000, y = diversity, group = block))+
  geom_vline(aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist/1000, y = diversity, col = group),alpha = 0.4)+
  geom_line(aes(x=dist/1000, y = mean_div))+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("Diversity")+
  #labs(col= "Generation")
  theme(legend.position = "none")
#ggtitle(paste0("Tajima's D (", ID, ")"))
ggsave(filename =paste0("plots/div_plot_2sps.jpg"), plot = div_plot, width = 25, height = 10)


theta_w_plot = ggplot(sum_stats[time %in% times], aes(x = dist/1000, y = theta_w_corr, group = block))+
  geom_vline(aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist/1000, y = theta_w_corr, col = group),alpha = 0.4)+
  geom_line(aes(x=dist/1000, y = mean_thetaw))+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("Watterson's Theta")+
  #labs(col= "Generation")
  theme(legend.position = "none")
#ggtitle(paste0("Tajima's D (", ID, ")"))
ggsave(filename =paste0("plots/thetaw_plot_2sps.jpg"), plot = theta_w_plot, width = 25, height = 10)

plot = ggarrange(taj_plot_2, div_plot, theta_w_plot)
ggsave(filename =paste0("plots/div_plot_2sps.jpg"), plot = div_plot, width = 20, height = 10)


h1_plot = ggplot(sum_stats[time %in% times], aes(x = dist/1000, y = H1, group = block))+
  geom_vline(aes(xintercept = 0), col= "black")+
  geom_line(aes(x = dist/1000, y = H1, col = group),alpha = 0.4)+
  geom_smooth(data = sum_stats[time %in% times],aes(group = time), col = "black")+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("H1")+
  #labs(col= "Generation")
  theme(legend.position = "none")


h12_plot = ggplot(sum_stats[time %in% times & run == 1], aes(x = dist/1000, y = H12, group = block))+
  geom_vline(aes(xintercept = 0), col= "black")+
  geom_line(aes(x = dist/1000, y = H12, col = group),alpha = 0.4)+
  #geom_smooth(data = sum_stats[time %in% times],aes(group = time), col = "black")+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("H12")+
  #labs(col= "Generation")
  theme(legend.position = "none")

h2h1_plot = ggplot(sum_stats[time %in% times], aes(x = dist/1000, y = H2H1, group = block))+
  geom_vline(aes(xintercept = 0), col= "black")+
  geom_line(aes(x = dist/1000, y = H2H1, col = group),alpha = 0.4)+
  #geom_smooth(data = sum_stats[time %in% times],aes(group = time), col = "black")+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("H2/H1")+
  #labs(col= "Generation")
  theme(legend.position = "none")





plot = ggarrange(allele_plot, taj_plot_2, heights = c(1,2))
ggsave(filename =paste0("initial_plot_all_groups_4.jpg"), plot = plot, width = 15, height = 10)


pdf(file=paste0("~/oliviaphd/seglift_win/initial_plot_",group,".pdf"), width = 10, height = 10) 
plot = ggarrange(plot_1, plot_2, heights = c(1,2))
#plot
dev.off()
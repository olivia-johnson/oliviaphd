library(data.table)
library(ggplot2)
library(viridis)
library(stringr)
library(egg)
setwd("~/oliviaphd/seglift_long")

groups = c("1", "2")

#ID = "L=10, d = 0.65, y=0.5"
  ##"L=10, d = 0.65, y=4"
  ##"No Fitness, L=10"
freq_data_100kg = NULL
sum_stats_100kg = NULL

for (g in groups){
  print(g)

  f_list_100kg <- list.files(path =paste0("~/oliviaphd/seglift_long/group_", g, "/"),pattern ="al_freq_")
  
  s_list_100kg <- list.files(path =paste0("~/oliviaphd/seglift_long/group_", g, "/"),pattern ="s_stat_2sps")
  
  
  ## collate al_freq files
  for (i in 1:(length(f_list_100kg))){  
    filename = f_list_100kg[i]
    #print(filename)
    al_freq_100kg = fread(file = paste0("~/oliviaphd/seglift_long/group_",g, "/",filename))
    al_freq_100kg[,run := i]
    al_freq_100kg[,group:=g]
    freq_data_100kg = rbind(freq_data_100kg, al_freq_100kg)
  }
  
  ## collate sum_stat files
  for (i in 1:(length(s_list_100kg))){  
    filename = s_list_100kg[i]
    #print(filename)
    s_stat_100kg = fread(file = paste0("~/oliviaphd/seglift_long/group_",g, "/",filename))
    s_stat_100kg[,run := i]
    s_stat_100kg[,group:=g]
    sum_stats_100kg = rbind(sum_stats_100kg, s_stat_100kg)
  }
}

freq_data_100kg[, mean_freq:=mean(mut_freq), by=c("Gen", "mut_pos", "group")]
#freq_data_100kg[, Season:= ifelse((Gen/2)%%15==0, "Winter", "Summer")]
freq_data_100kg[, block :=paste0(group, "_", run, "_",(mut_pos/1e5)+.5)]

sum_stats_100kg[, win_pos:=n_win%%21]
sum_stats_100kg[, theta_w_corr:= theta_w*1e15]
sum_stats_100kg[, mean_tajd:=mean(tajimas_d_branch), by=c("time", "win_pos", "group")]
sum_stats_100kg[, mean_div:=mean(diversity), by=c("time", "win_pos", "group")]
sum_stats_100kg[, mean_thetaw:=mean(theta_w_corr), by=c("time", "win_pos", "group")]
sum_stats_100kg[, midpoint:= (win_end-win_start)/2 + win_start]
sum_stats_100kg[, block :=paste0(group, "_", run, "_",(n_win%/%21))]
sum_stats_100kg[, dist:=49999-(midpoint-(n_win%/%21)*1e5), by="n_win"]
setnames(sum_stats_100kg, "H2/H1", "H2H1")



plot = ggplot(freq_data_100kg[run==4],aes( x = Gen, y= mut_freq))+
  geom_line(aes( col = as.factor(mut_pos))) +
  ggtitle("Allele Frequency (Single Simulation, N = 1000, d=0.65)") +
  xlab("Generation")+
  ylab("Allele Frequency")

##ggsave(filename =("freq_plot_0_3.jpg"), plot = plot, width = 15, height = 10)

group.labs <- c("1"= "y=0.5, d=0.65", "2" = "y=4, d=0.65")


allele_plot_100kg=ggplot(freq_data_100kg,aes( x = Gen, y= mut_freq, group = block))+
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

times_100kg = c(8, 15, 23,12000, 12008, 12015, 12023,42000,42008,42015, 42023, 81000, 81008,81015,81023, 117000, 117008, 117015, 117023)
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

taj_plot_2_100kg = ggplot(sum_stats_100kg[time %in% times_100kg], aes(x = dist/1000, y = tajimas_d_branch, group = block))+
  geom_vline(aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist/1000, y = tajimas_d_branch, col = group),alpha = 0.4)+
  geom_line(aes(x=dist/1000, y = mean_tajd))+
  #geom_smooth(data = sum_stats_100kg[time %in% times_100kg],aes(group = time), col = "black")+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("Tajima's D")+
  #labs(col= "Generation")
  theme(legend.position = "none")
  #ggtitle(paste0("Tajima's D (", ID, ")"))

div_plot_100kg = ggplot(sum_stats_100kg[time %in% times_100kg], aes(x = dist/1000, y = diversity, group = block))+
  geom_vline(aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist/1000, y = diversity, col = group),alpha = 0.4)+
  geom_line(aes(x=dist/1000, y = mean_div))+
  #geom_smooth(data = sum_stats_100kg[time %in% times_100kg],aes(group = time), col = "black")+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("Diversity")+
  #labs(col= "Generation")
  theme(legend.position = "none")
#ggtitle(paste0("Tajima's D (", ID, ")"))

theta_w_plot_100kg = ggplot(sum_stats_100kg[time %in% times_100kg], aes(x = dist/1000, y = theta_w_corr, group = block))+
  geom_vline(aes(xintercept = 0), col= "red")+
  geom_line(aes(x = dist/1000, y = theta_w_corr, col = group),alpha = 0.4)+
  #geom_smooth(data = sum_stats_100kg[time %in% times_100kg],aes(group = time), col = "black")+
  geom_line(aes(x=dist/1000, y = mean_thetaw))+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("Watterson's Theta")+
  #labs(col= "Generation")
  theme(legend.position = "none")
#ggtitle(paste0("Tajima's D (", ID, ")"))
ggsave(filename =paste0("plots/theta_plot_2_100kg.jpg"), plot = theta_w_plot_100kg, width = 15, height = 10)

plot = ggarrange(taj_plot_2_100kg, div_plot_100kg, theta_w_plot_100kg, nrow = 3)
ggsave(filename =paste0("plots/summary_plot_both.jpg"), plot = plot, width = 25, height = 20)


h1_plot_100kg = ggplot(sum_stats_100kg[time %in% times & block == "1_1_6"], aes(x = dist/1000, y = H1, group = block))+
  geom_vline(aes(xintercept = 0), col= "black")+
  geom_dotplot(aes(x = dist/1000, y = H1, col = group),alpha = 0.4)+
  #geom_smooth(data = sum_stats[time %in% times],aes(group = time), col = "black")+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("H1")+
  #labs(col= "Generation")
  theme(legend.position = "none")


h12_plot_100kg = ggplot(sum_stats_100kg[time %in% times & run == 1], aes(x = dist/1000, y = H12, group = block))+
  geom_vline(aes(xintercept = 0), col= "black")+
  geom_jitter(aes(x = dist/1000, y = H12, col = group),alpha = 0.4)+
  #geom_smooth(data = sum_stats[time %in% times],aes(group = time), col = "black")+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("H1")+
  #labs(col= "Generation")
  theme(legend.position = "none")

h2h1_plot_100kg = ggplot(sum_stats_100kg[time %in% times], aes(x = dist/1000, y = H2H1, group = block))+
  geom_vline(aes(xintercept = 0), col= "black")+
  geom_jitter(aes(x = dist/1000, y = H2H1, col = group),alpha = 0.4)+
  #geom_smooth(data = sum_stats[time %in% times],aes(group = time), col = "black")+
  facet_grid(group~time, labeller = labeller(group = group.labs))+
  xlab("Distance from selected loci (kb)")+
  ylab("H1")+
  #labs(col= "Generation")
  theme(legend.position = "none")





plot = ggarrange(allele_plot_100kg, taj_plot_2_100kg, heights = c(1,2))
ggsave(filename =paste0("initial_plot_all_groups_4.jpg"), plot = plot, width = 15, height = 10)


pdf(file=paste0("~/oliviaphd/seglift_win/initial_plot_",group,".pdf"), width = 10, height = 10) 
plot = ggarrange(plot_1, plot_2, heights = c(1,2))
#plot
dev.off()
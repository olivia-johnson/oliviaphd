combined_plots
library(ggpubr)

#### CHECK times AND exp line in div plot ###

fplotdata=freq_data[d==0.65 & y==4 & fit_type=="AL"]
times=c(45007, 45014, 45022, 45029)  ## EG
#times=c(45006,45012,45013,45014)  ##UG
plotdata=sum_stat[d==0.65 & y==4 & fit_type=="AL" & time %in% times]

seg_labels=unique(fplotdata[Freq.bin=="Segregating" & time==50000, .(label, n_seg)])

allele_plot =ggplot(data=fplotdata,aes(x = time, y= mut_freq))+
  geom_line(aes(group=block), col = "seagreen",alpha = 0.3) +
  #facet_wrap(~label)+
  #ggtitle(paste0("Allele Frequency (", ID, ")")) +
  xlab("Generation")+
  scale_y_continuous(limits = c(0, 1.0),expand = c(0,0))+
  ylab("Allele Frequency")+
  theme(legend.position = "none")+
  labs(col= "Seasonal Mutation")
  #ggtitle("Allele Frequencies")

#allele_plot = allele_plot + geom_text(data=seg_labels, x=45000, y=0.1, aes(label=n_seg), parse=TRUE)


# New facet label names for dose variable
time.labs <- c("Summer\nMid", "Summer\nEnd", "Winter\nMid", "Winter\nEnd")
names(time.labs) <- c("45007", "45014", "45022", "45029")

# New facet label names for supp variable
seg.labs <- c("Neutral", "Selected")
names(seg.labs) <- c(FALSE, TRUE)

tajd_plot=ggplot(plotdata, aes(x = chrom_pos/100000, y = tajimas_d_branch, group = block))+
  geom_vline(data = subset(plotdata,segregating==T),aes(xintercept = (chromsize/2)/100000), col= "red")+
  geom_line(aes(x = chrom_pos/100000, y = tajimas_d_branch), col = "seagreen",alpha = 0.4)+
  geom_line(aes(x=chrom_pos/100000, y = mean_tajd))+
  geom_hline(yintercept = 0, col = "sienna1", )+
  facet_grid(segregating~time,labeller = labeller(time = time.labs, segregating = seg.labs))+
  theme(strip.text.y = element_text(size = 10))+
  ylab("Tajima's D")+
  xlab("Distance along chromosome (100kb)")+
  #ggtitle("Tajima's D")+
  theme(legend.position = "none")

div_plot = ggplot(plotdata, aes(x = chrom_pos/100000, y = diversity, group = block))+
  geom_vline(data = subset(plotdata,segregating==T),aes(xintercept = (chromsize/2)/100000), col= "red")+
  geom_line(aes(x = chrom_pos/100000, y = diversity), col = "seagreen",alpha = 0.4)+
  geom_hline(yintercept = 0.00012, col = "sienna1", )+ ##CP
  #geom_hline(yintercept = 0.00005, col = "grey", )+ ##FP_UG
  #geom_hline(yintercept = 0.000021, col = "grey", )+ ##FP_EG
  geom_line(aes(x=chrom_pos/100000, y = mean_div))+ 
  facet_grid(segregating~time, labeller = labeller(time = time.labs, segregating = seg.labs))+
  theme(strip.text.y = element_text(size = 10))+
  ylab("Diversity (pi)")+
  xlab("Distance along chromosome (100kb)")+
  #ggtitle("Diversity (pi)")+
  theme(legend.position = "none")

plot=ggarrange( div_plot,tajd_plot, ncol=2)

ggexport(allele_plot, plot, filename = "new_plots/decrease_diveristy.pdf",
         nrow = 2, ncol = 1)

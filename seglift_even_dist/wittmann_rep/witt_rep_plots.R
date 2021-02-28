library(data.table)
library(ggplot2)
library(viridis)
library(stringr)
library(egg)
setwd("~/oliviaphd/seglift_even_dist/wittmann_rep")

group = 1

f_list <- list.files(path =paste0("~/oliviaphd/seglift_even_dist/wittmann_rep/ed_", group, "/"),pattern ="al_freq_")

freq_data = NULL
## collate al_freq files
for (i in 1:(length(f_list))){  
  filename = f_list[i]
  #print(filename)
  al_freq = fread(file = paste0("~/oliviaphd/seglift_even_dist/wittmann_rep/ed_",group, "/",filename))
  al_freq[,run := i]
  freq_data = rbind(freq_data, al_freq)
}

freq_data_1 = freq_data
freq_data_1[, mean_freq:=mean(mut_freq), by=c("Gen", "mut_pos")]
freq_data_1[, Season:= ifelse(Gen%%30>=16, "Winter", "Summer")]
freq_data_1[, winter_al:=1-mut_freq]

sample = c(1, 10, 20, 30, 40, 50, 60, 70,80, 90,100)
plot_1 = ggplot(freq_data_1[mut_pos %in% sample & run==1],aes( x = Gen, y= mut_freq))+
  geom_line(aes( col = as.factor(mut_pos)), alpha = 0.5) +
  ggtitle("d=0.15") +
  #xlab("Generation")+
  ylab("Mean Allele Frequency")+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(limits = c(0, 1.0))+
  theme(legend.position = "none")
plot_2=ggplot(freq_data_2[mut_pos %in% sample & run==1],aes( x = Gen, y= mean_freq))+
  geom_line(aes( col = as.factor(mut_pos)), alpha = 0.5)+
  ggtitle( "d=0.5")+
  xlab("Generation")+
  scale_y_continuous(limits = c(0, 1.0))+
 # ylab("Mean Allele Frequency")+
  theme(axis.title.y = element_blank())+
  theme(legend.position = "none")
plot_3=ggplot(freq_data_3[mut_pos %in% sample & run==1],aes( x = Gen, y= mean_freq))+
  geom_line(aes( col = as.factor(mut_pos)), alpha = 0.5)+
  ggtitle( "d=0.65")+
  #xlab("Generation")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  #ylab("Mean Allele Frequency")+
  scale_y_continuous(limits = c(0, 1.0))+
  labs(col = "Selected Loci")

## FIG 7 150 gens

plot_1 = ggplot(freq_data_1[mut_pos %in% sample & run == 2 & Gen >=3000 & Gen<=3150],aes( x = Gen, y= mut_freq))+
  geom_line(aes( col = as.factor(mut_pos)), alpha = 0.8) +
  ggtitle("d=0.15") +
  scale_x_continuous(breaks = c(3000, 3015, 3030, 3045, 3060, 3075, 3090, 3105, 3120, 3135, 3150),
                     labels = c("Winter","Summer", "Winter","Summer", "Winter","Summer", "Winter","Summer", 
                                "Winter","Summer", "Winter"))+
  theme(axis.text.x=element_text(angle = 90))+
  #xlab("Generation")+
  ylab("Allele Frequency")+
  scale_y_continuous(limits = c(0, 1.0))+
  theme(axis.title.x = element_blank())+
  theme(legend.position = "none")

plot_2 = ggplot(freq_data_2[mut_pos %in% sample & run == 2 & Gen >=3000 & Gen<=3150],aes( x = Gen, y= mut_freq))+
  geom_line(aes( col = as.factor(mut_pos)), alpha = 0.8) +
  ggtitle("d=0.5") +
  scale_x_continuous(breaks = c(3000, 3015, 3030, 3045, 3060, 3075, 3090, 3105, 3120, 3135, 3150),
                     labels = c("Winter","Summer", "Winter","Summer", "Winter","Summer", "Winter","Summer", 
                                "Winter","Summer", "Winter"))+
  theme(axis.text.x=element_text(angle = 90))+
  xlab("Generation")+
  scale_y_continuous(limits = c(0, 1.0))+
  # ylab("Allele Frequency")+
  theme(axis.title.y = element_blank())+
  theme(legend.position = "none")
plot_3 = ggplot(freq_data_3[mut_pos %in% sample& run == 2 & Gen >=3000 & Gen<=3150],aes( x = Gen, y= mut_freq))+
  geom_line(aes( col = as.factor(mut_pos)), alpha = 0.8) +
  ggtitle("d=0.65") +
  scale_x_continuous(breaks = c(3000, 3015, 3030, 3045, 3060, 3075, 3090, 3105, 3120, 3135, 3150),
                     labels = c("Winter","Summer", "Winter","Summer", "Winter","Summer", "Winter","Summer", 
                                "Winter","Summer", "Winter"))+
  theme(axis.text.x=element_text(angle = 90))+
  #ylab("Allele Frequency")+
  #xlab("Generation")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  scale_y_continuous(limits = c(0, 1.0))+
  labs(col = "Selected Loci")
plot = ggarrange(plot_1, plot_2, plot_3, nrow = 1)

ggsave(filename = "fig7_exact.jpg", plot = plot, width = 15, height = 5)



##[mut_pos %in% sample & run==1]
plot_1 = ggplot(freq_data_1[run == 2 ],aes( x = Gen, y= mut_freq))+
  geom_line(aes( col = as.factor(mut_pos)), alpha = 0.2) +
  ggtitle("d=0.15") +
  xlab("Generation")+
  ylab("Allele Frequency")+
  scale_y_continuous(limits = c(0, 1.0))+
  theme(axis.title.y = element_blank())+
  theme(legend.position = "none")

plot_2 = ggplot(freq_data_2[run == 2],aes( x = Gen, y= mut_freq))+
  geom_line(aes( col = as.factor(mut_pos)), alpha = 0.2) +
  ggtitle("d=0.5") +
  #xlab("Generation")+
  scale_y_continuous(limits = c(0, 1.0))+
  # ylab("Allele Frequency")+
  theme(axis.title.y = element_blank())+
  theme(legend.position = "none")
plot_3 = ggplot(freq_data_3[run == 2],aes( x = Gen, y= mut_freq))+
  geom_line(aes( col = as.factor(mut_pos)), alpha = 0.2) +
  ggtitle("d=0.65") +
  #ylab("Allele Frequency")+
  #xlab("Generation")+
  theme(axis.title.x = element_blank())+
  theme(axis.title.y = element_blank())+
  scale_y_continuous(limits = c(0, 1.0))+
  labs(col = "Selected Loci")


ggsave(filename = "fig7_onerun2.jpg", plot = plot, width = 15, height = 5)
pdf(file="~/oliviaphd/seglift_even_dist/wittmann_rep/fig7_onerun2.pdf", width = 15, height = 5) 
plot = ggarrange(plot_1, plot_2, plot_3, nrow = 1)
dev.off()

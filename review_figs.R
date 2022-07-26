library(data.table)
library(ggpubr)
library(viridis)
library(tidyverse)
library(ggdist)
library(readxl)
library(grDevices)
setwd("~/phd_data/drosData")

data.labs = c('chrom', 'pos','AF', 'SP', 'SQ' ,'FallF', 'SprF')
bergland_snps= fread("bergland/bergland_PA.txt", col.names = data.labs)

load("machado/mel_freqdp_042016_Ne_fixed_correctBAVI.Rdata")
machado_data=as.data.table(cbind(info, freq))

m_snps<- as.data.table(read_excel("machado/elife-67577-supp1-v2.xlsx", 
                                   sheet = "Supplementaryfile1B"))
m_pops <-as.data.table(read_excel("machado/elife-67577-supp1-v2.xlsx", 
                    sheet = "Supplementaryfile1A"))


pop_keep = m_pops[Core20=="yes", .(.N, Sample, InternalName, Year, Season), by=c("Locality")]
pop_keep = pop_keep[N>2]
pops = unique(pop_keep$Locality)
for (i in pops){
  pop_keep[Locality==i & min(Year), yr:=ifelse(Year==min(Year), 1, 2)]
}
pop_keep[, time:=ifelse(Season=="spring", paste("s",yr, sep="_"), paste("f",yr, sep="_")), by=c("yr", "Season")]
pop_keep[, label:=paste(Locality, time, sep="_"), by= c("Locality", "time")]
pop_keep=pop_keep[Locality!= "MA_la" & Locality!="WI_cp"]
mpops= m_pops[Sample%in%pop_keep$Sample]
m_pops=m_pops[Locality!= "MA_la" & Locality!="WI_cp"]
cols=c("X.CHROM", "POS", pop_keep$InternalName)
col.names=c("chrom", "pos", pop_keep$label)
machado=machado_data[, ..cols]
setnames(machado, cols, col.names)

# b_pops=m_pops[Publication=="Bergland et al 2014" & Locality_code =="PA" & Year!=2012 & Sample!="PA_li_11_frost"]
# b_pops[, yr:=(Year%%2009)+1, by="Year"]
# b_pops[, label:=ifelse(Season=="spring", paste("s", yr, sep="_"), paste("f",yr, sep="_")), by=c("yr", "Season")]
# cols=c("X.CHROM", "POS", b_pops$InternalName)
# col.names=c("chrom", "pos", b_pops$label)
# bergland_data=machado_data[, ..cols]
# setnames(bergland_data, cols, col.names)
# bergland= merge(bergland_snps, bergland_data, all.x=TRUE, by=c("chrom", "pos"))
# bergland=na.omit(bergland)
# bergland[, id:=paste(chrom, pos, sep="_")]
# bergland[, fc_1:=(`s_1`-`f_1`), by=c("id")]
# bergland[, fc_2:=(`s_2`-`f_2`), by=c("id")]
# bfreq_change=bergland[, freqchange:=abs((fc_1+fc_2)/2), by=c("id")]
# # machado.af=machado.af[freqchange>0.04,.(id, m1S, m1F, m2S, m2F)]
# bergland=bergland[freqchange>0,.(id, fc_1, s_1, f_1, s_2, f_2,s_3, f_3)]
# bergland[fc_1<0, s_1:=(1-s_1), by=c("id")]
# bergland[fc_1<0, f_1:=1-f_1, by=c("id")]
# bergland[fc_1<0, s_2:=1-s_2, by=c("id")]
# bergland[fc_1<0, f_2:=1-f_2, by=c("id")]
# bergland_freq = melt(bergland, id = 1:2, variable.name = "time", value.name = "freq",variable_factor=TRUE)

snp_keep = m_snps[, .(chrom, pos)]
machado=merge(snp_keep, machado, all.x = TRUE, by=c("chrom", "pos"))
machado[, id:=paste(chrom,pos, sep="_"), by=c("chrom", "pos")]
fintersect(bergland_snps[, .(chrom, pos)],snp_keep)
machado=na.omit(machado)


fc = machado[, .(id, pa.fc=PA_li_f_1-PA_li_s_1, ca.fc=CA_es_f_1-CA_es_s_1, vi.fc=VA_ch_f_1-VA_ch_s_1)]
fc[, m_fc:=(pa.fc+ ca.fc+ vi.fc)/3, by="id"]
pol = fc[m_fc<0, id]

machado.af = melt(machado, variable.name = "locality", measure.vars = patterns("s_1", "f_1", "s_2","f_2"), value.name=c("s_1", "f_1", "s_2","f_2"), variable.factor = FALSE)
machado.af=na.omit(machado.af)
machado.af[, fc_1:=(`s_1`-`f_1`), by=c("id", "locality")]
machado.af[, fc_2:=(`s_2`-`f_2`), by=c("id", "locality")]
machado.af=machado.af[,.(id, locality, fc_1, s_1, f_1, s_2, f_2)]
# machado.af[fc_1>0, s_1:=(1-s_1), by=c("id", "locality")]
# machado.af[fc_1>0, f_1:=1-f_1, by=c("id", "locality")]
# machado.af[fc_1>0, s_2:=1-s_2, by=c("id", "locality")]
# machado.af[fc_1>0, f_2:=1-f_2, by=c("id", "locality")]
machado.af[id %in% pol, s_1:=(1-s_1), by=c("id", "locality")]
machado.af[id %in% pol, f_1:=1-f_1, by=c("id", "locality")]
machado.af[id %in% pol, s_2:=1-s_2, by=c("id", "locality")]
machado.af[id %in% pol, f_2:=1-f_2, by=c("id", "locality")]

machado.af2 = melt(machado.af[locality=="1" |locality=="2" | locality=="3"], id = 1:3, variable.name = "time", value.name = "freq",variable_factor=TRUE)
machado.af2$locality <- factor(machado.af2$locality, labels = c("Pennsylvania","California", "Virginia"))
machado.af2=distinct(machado.af2)

machado.af[, FF:=((f_1+ f_2)/2), by=c("id", "locality")]
machado.af[, SF:=((s_1+ s_2)/2), by=c("id", "locality")]
machado.af[, amp:=abs(FF-SF), by=c("id", "locality")]




plot = ggplot()+
  geom_line(data=machado.af2,aes(x=time, y=freq, group=id, col=as.factor(locality)), alpha=0.03)+
  facet_wrap("locality", labeller = label_value, nrow=1)+ 
  xlab("Season")+
  ylab("Allele Frequency")+
  scale_x_discrete(labels=c("s_1"="Spring", "f_1"="Fall", "s_2"="Spring", "f_2"="Fall"), limits=c("s_1", "f_1", "s_2", "f_2"))+
  theme(legend.position = "none")

## with bergland
bergland_freq[, locality:="Bergland et al 2014"]
plot = ggplot()+
  geom_line(data=machado.af2,aes(x=time, y=freq, group=id, col=as.factor(locality)), alpha=0.03)+
  # geom_line(data=bergland_freq, aes(x=time, y=freq, group=id, col=as.factor(locality)), alpha=0.05)+
  facet_wrap("locality", labeller = label_value, nrow=1)+ 
  xlab("Season")+
  ylab("Allele Frequency")+
  scale_x_discrete(labels=c("s_1"="Spring", "f_1"="Fall", "s_2"="Spring", "f_2"="Fall", "s_3"="Spring", "f_3"="Fall"), limits=c("s_1", "f_1", "s_2", "f_2", "s_3", "f_3"))+
  theme(legend.position = "none")
  
plot = ggplot()+
  geom_line(data=machado.af2,aes(x=time, y=freq, group=id, col=locality), alpha=0.025)+  xlab("Season")+
  facet_wrap("locality", labeller = label_value, nrow=1)+ 
  theme_bw()+
  xlab("Season")+
  ylab("Allele Frequency")+
  scale_x_discrete(breaks=c("s_1", "f_1", "s_2", "f_2"), 
  labels=c("Spring", "Fall", "Spring", "Fall"), 
  expand=c(0, 0) )+
  scale_y_continuous(limits = c(0,1), expand=c(0,0))+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0, hjust=1), panel.spacing.x = unit(1, "lines"))


ggexport(plot, filename="plots/fig_1_pol2.pdf", width = 6, height =4)
ggsave(filename =paste0("plots/fig_1_pol.jpg"), plot = plot , width = 6, height = 4)



### FIGURE 3
setwd("~/phd_data/review")
groups = c(1:7)


alfreq_data = NULL
dropcols=c("s_d", "w_d", "s_fx", "w_fx")

for (g in groups){
  print(g)
  
  f_list  <- list.files(path =paste0("~/phd_data/review/group_", g, "/"),pattern ="al_freq_")
  
  
  parameters <- fread(file=paste0("~/phd_data/review/group_",g, "/parameters.txt"), sep = ":")
  setkey(parameters, V1)
  
  fiton=parameters["fitness_on", V2]
  if (fiton==0){
    epi=0
  } else{
    epi = parameters["y", V2]
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
  loci =parameters["l", V2]
  
  g_label = ifelse(fiton==0, "No Fitness", epi)
  
  for (i in 1:(length(f_list ))){
    ## collate al_freq files
    
    filename = f_list [i]
    #print(filename)
    al_freq  = fread(file = paste0("~/phd_data/review/group_",g, "/",filename), fill=TRUE,skip="Gen")
    al_freq [,run := i]
    al_freq [,group:=g]
    al_freq [,y:=epi]
    al_freq [,l:=loci]
    al_freq [,fit:=fiton]
    al_freq [,s_gen:=sum_gen]
    al_freq [,w_gen:=win_gen]
    al_freq [,s_pop:=sum_pop]
    al_freq [,w_pop:=win_pop]
    al_freq[, pop_season := pop_s]
    al_freq[, gen_season := gen_s]
    al_freq[,label := g_label]
    alfreq_data  = rbind(alfreq_data , al_freq )
  }
}
alfreq_data[gen_season == "EG", gen_year:=Gen%%30, by=c("label", "Gen")]
alfreq_data[gen_season == "UG", gen_year:=Gen%%15, by=c("label", "Gen")] 
alfreq_data[gen_season == "EG", season:=ifelse((gen_year<16 & gen_year>0), "summer", "winter"), by=c("Gen")]
alfreq_data[gen_season == "UG", season:=ifelse((gen_year<14 & gen_year>0), "summer", "winter"), by=c("Gen")]
alfreq_data [, id :=paste0(group, "_", run, "_",mut_pos)]
alfreq_data[mut_freq!=1, Freq.bin:="Segregating"]
alfreq_data[mut_freq==1, Freq.bin:="Fixed_Summer"]
alfreq_data[Freq.bin=="Segregating",n_seg:=.N, by = c("Gen", "label","l")]
alfreq_data [, id :=paste0(group, "_", run, "_",mut_pos)]
sd_val = na.omit(alfreq_data[,.(s_d, w_d, s_fx, w_fx), by=c("id")])
rm.col=c("s_d", "w_d", "s_fx", "w_fx")
alfreq_data[, (rm.col):=NULL]
alfreq_data[, year:= ifelse(gen_season=="EG", (Gen%/%30)+1, (Gen%/%15)+1), by=c("Gen", "gen_season")]
alfreq_data[Freq.bin=="Segregating",fc_year:=max(mut_freq)-min(mut_freq), by=c("id", "year")]
alfreq_data[Freq.bin=="Segregating",fc:=mean(fc_year), by=c("id")]
alfreq_data[Gen==90015, fluctuating:=ifelse(mut_freq!=0 & mut_freq != 1, T, F), by="id"]
intermediate = alfreq_data[Gen > 20000 & Freq.bin=="Segregating"& l>10, .(ymax = max(mut_freq), ymin = min(mut_freq)), by=c("id", "year","group", "label") ]
intermediate[, yamp:=ymax-ymin, by=c("id", "year","group", "label") ]
amps = intermediate[yamp>0, .(av_max = mean(ymax), av_min = mean(ymin), av_amp = mean(yamp)), by=c("id","group", "label")]
amps[, amp := av_max-av_min, c("id", "group","label")]
amps = merge(amps, sd_val, by = "id", all.x=TRUE)

bergland_snps[SQ<0.3& chrom!= "X", amp:=abs(FallF-SprF)]
bergland_snps[SQ<0.3& chrom!= "X", id:=paste(chrom, pos, sep="_"), by=c("chrom", "pos")]
bergland_snps[SQ<0.3& chrom!= "X", label:="Bergland"]
mean_berg_amp=mean(bergland_snps$amp)
mean_machado_amp= mean(machado.af$amp)
b_sd = sd(bergland_snps$amp)

b_amp=bergland_snps[SQ<0.3& chrom!= "X", .(id, amp, label)]
m_amp = machado.af[, .(id, amp)]
m_amp[,label:="Machado"]
# m_amp=unique(m_amp[, amp:=mean(amp), by="id"])
amplitudes = rbind(amps[, .(id, label, amp)], m_amp, b_amp)
amplitudes[label=="No Fitness", label:="Drift"]

numloci=amplitudes[, .N, by="label"]
amplitudes[, data_type:=ifelse(label=="Bergland" | label=="Machado", "Empirical", "Simulation"), by="label"]

fig3A = ggplot(amplitudes[data_type=="Simulation"], aes(x = label)) + 
  geom_jitter(aes(y = amp,col = label), alpha=0.2, width = 0.38
  )  +
  geom_boxplot(aes(y = amp),alpha=0,width = 0.8,
    outlier.color = NA ## `outlier.shape = NA` works as well ## remove outliers
  ) +
  # coord_cartesian(xlim = c(1.2, NA)) +
  theme_light()+
  scale_x_discrete(limits=c("Drift", "0.5", "1", "4", "8", "12","20"), breaks=c("Drift", "0.5", "1", "4", "8", "12","20"))+
  scale_y_continuous(limits=c(-0.035, 1))+
  ylab("Amplitude of Fluctuations")+
  xlab("y")+
  annotate("text", x= 1.75, y=0.8, label ="italic(w)(italic(z)) == (1 + italic(z))^italic( y)",
           parse=TRUE, vjust = -1, size=7, family="Times")+
  geom_text(data = numloci, aes(x = label, label=paste0("n=",N), y=-0.03), size = 3.5)+
  theme(legend.position = "none", 
        axis.title.x = element_text(face = "italic", size=18, family="Times"), 
        axis.text = element_text(size=16), 
        axis.title.y = element_text(size=18))

fig3B=ggplot(amplitudes[data_type=="Empirical"], aes(x = label)) + 
  geom_jitter(aes(y = amp,col = label), alpha=0.2, width = 0.38
  )  +
  geom_boxplot(aes(y = amp),alpha=0,width = 0.8,
               outlier.color = NA ## `outlier.shape = NA` works as well ## remove outliers
  ) +
  # coord_cartesian(xlim = c(1.2, NA)) +
  theme_light()+
  scale_x_discrete(limits=c("Bergland", "Machado"), breaks=c("Bergland", "Machado"))+
  xlab("Empirical Studies")+
  scale_y_continuous(limits=c(-0.035, 1))+
  geom_text(data = numloci, aes(x = label, label=paste0("n=",N), y=-0.03), size = 3.5)+
  theme(legend.position = "none", axis.title.x = element_text(size=16),
        axis.text = element_text(size=16), axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

fig3 = ggarrange(fig3A, fig3B, widths = c(2.5, 1))

fig3 = ggplot(amplitudes, aes(x = label, y = amp, col = label)) + 
  ## add half-violin from {ggdist} package
  ggdist::stat_halfeye(
    aes(fill=label),
    adjust = .6,            ## custom bandwidth
    width = .8,             ## adjust height
    justification = -.2,    ## move geom to the right
    .width = 0              ## remove slab interval
  ) + 
  geom_boxplot(
    width = .12, 
    outlier.color = NA ## `outlier.shape = NA` works as well ## remove outliers
  ) +

  gghalves::geom_half_point(    ## add dot plots from {ggdist} package
    side = "l", ## draw jitter on the left
    range_scale = .4,## control range of jitter
    alpha = .3, ## add some transparency
    size=0.2
  )  +
  # coord_cartesian(xlim = c(1.2, NA)) +
  theme_bw()+
  scale_x_discrete(limits=c("0.5", "1", "4", "8", "12","20", "No Fitness", "Bergland", "Machado"))+
  ylab("Amplitude of Fluctuations")+
  xlab("Epistasis")+
  coord_flip()+
  theme(legend.position = "none") 






dist_compare = ggplot(data=amplitudes, aes(x=factor(label), y=amp))+
  #geom_violin(scale = "area")+
  #geom_boxplot(width=0.1, alpha = 0.5, outlier.shape = NA)+
  geom_boxplot()+
  # geom_boxplot(data=bergland_snps[SQ<0.3 & chrom!= "X"], aes(y=amp, x="Bergland"))+
  # geom_boxplot(data=machado.af[amp>0 & locality=="1"], aes(y=amp, x="Machado"))+
  geom_hline(aes(yintercept=mean_berg_amp))+
  geom_text(aes(1.5,mean_berg_amp,label = "Bergland et al", vjust = -1),size=3.25)+
  geom_hline(aes(yintercept=mean_machado_amp))+
  geom_text(aes(1.5,mean_machado_amp,label = "Machado et al", vjust = -1),size=3.25)+
  scale_x_discrete(limits=c("0.5", "1", "4", "8", "12","20", "No Fitness"))+
  scale_x_discrete(linits=c("Bergland", "Machado"))+
  ylab("Amplitude of Fluctuations")+
  xlab("Epistasis")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) 

ggexport(fig3, filename="fig3_box7.pdf", width = 11, height =6)
ggsave(filename =paste0("fig3_box7.jpg"), plot = fig3 , width = 11, height = 6)


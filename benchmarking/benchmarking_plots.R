## Benchmarking  ##
library(data.table)
library(ggplot2)
library(lubridate)
library(ggpubr)
library(writexl)
library(scales)

setwd("~/phd_data/benchmarking/")

##Burn-in
msprime_ts=fread("~/phd_data/benchmarking/burnin_benchmarking_msprime2.txt")
# msprime_ts[, sim_type:="msprime_ts"]
msprime_ts[, simulator:="msprime"]
msprime_ts[, data_type:=ifelse(sim_type=="msprime_ts", "Tree Sequence","Tree Sequence\nwith Neutral Mutations")]
msprime_ts[, final_gen:=NA]

slim_ts=fread("~/phd_data/benchmarking/burnin_benchmarking_slim_ts.txt")
slim_ts[, start_time:=period_to_seconds(hms(start_time))]
slim_ts[,end_time:=period_to_seconds(hms(end_time))]
slim_ts[, time:=ifelse(end_time>start_time, end_time-start_time,((86400-start_time)+end_time)) ]
slim_ts[, sim_type:="slim_ts"]
slim_ts[, simulator:="SLiM"]
slim_ts[, data_type:="Tree Sequence\n(10Ne)"]
slim_ts[, final_gen:=100000]


slim_tscc=fread("~/phd_data/benchmarking/burnin_benchmarking_slim_tscc.txt")
slim_tscc[, start_time:=period_to_seconds(hms(start_time))]
slim_tscc[,end_time:=period_to_seconds(hms(end_time))]
slim_tscc[, time:=ifelse(end_time>start_time, end_time-start_time,((86400-start_time)+end_time)) ]
slim_tscc[, sim_type:="slim_tscc"]
slim_tscc[, simulator:="SLiM"]
slim_tscc[, data_type:="Tree Sequence\n(Coalesced)"]


slim_mut=fread("~/phd_data/benchmarking/burnin_benchmarking_slim_mut.txt")
slim_mut[, start_time:=period_to_seconds(hms(start_time))]
slim_mut[,end_time:=period_to_seconds(hms(end_time))]
slim_mut[, time:=ifelse(end_time>start_time, end_time-start_time,((86400-start_time)+end_time))]
slim_mut[, sim_type:="slim_mut"]
slim_mut[, simulator:="SLiM"]
slim_mut[, data_type:="Classical\n(10Ne)"]
slim_mut[, final_gen:=100000]



data=rbind(msprime_ts,slim_ts[,.(time, memory, diversity,sim_type, simulator, data_type, final_gen)],slim_tscc[,.(time, memory, diversity,sim_type, simulator, data_type, final_gen)], slim_mut[,.(time, memory, diversity, sim_type, simulator, data_type,final_gen)])
data[, memory:=memory/1e6]  ## convert memory to MB
# data[simulator=="msprime", data_type:=ifelse(sim_type=="msprime_mut", "Tree Sequence with Neutral Mutation", "Tree Sequence"), by="sim_type"]

data[, `:=` (mean_time = mean(time), mean_memory=mean(memory), var_time=var(time), var_mem=var(memory), max_time=max(time), min_time=min(time), max_mem=max(memory), min_mem=min(memory), mean_div=mean(diversity), max_div=max(diversity), min_div=min(diversity),var_div=var(diversity), var_fg=var(final_gen), mean_fg=mean(final_gen), min_fg=min(final_gen), max_fg=max(final_gen)), by="sim_type"]

averages = unique(data[,.(simulator, data_type, mean_time, mean_memory, var_time, var_mem, max_mem, min_mem, max_time, min_time, mean_div, min_div, max_div, var_div, mean_fg, var_fg, min_fg, max_fg)])

## pvals
  ##msprime memory
t.test(data[sim_type=="msprime_ts", memory], data[sim_type=="msprime_mut", memory])
##msprime time
t.test(data[sim_type=="msprime_ts", time], data[sim_type=="msprime_mut", time])


##diversity p-vals
theta=4*10000*1e-7
expected_div=theta/(1+theta)
t.test(data[sim_type=="msprime_ts", diversity], mu=expected_div)
t.test(data[sim_type=="msprime_mut", diversity], mu=expected_div)
t.test(data[sim_type=="slim_ts", diversity], mu=expected_div)
t.test(data[sim_type=="slim_tscc", diversity], mu=expected_div)
t.test(data[sim_type=="slim_mut", diversity], mu=expected_div)


##plots##
mem = ggplot(data, aes(x=data_type))+
  # geom_jitter(aes(y = memory,col = sim_type), width = 0.35, size=2.2)  +
  # geom_boxplot(aes(y = memory),alpha=0,width = 0.8,outlier.color = NA )+
  geom_boxplot(aes( y=memory, fill=sim_type))+
  labs(x="Simulation Type",y="Memory (MB)") +theme_bw()+theme(legend.position = "none")+ facet_wrap("simulator", scale="free", drop=TRUE)+# scale_y_continuous(breaks=c(72:78,2500,5000,7500, 10000, 12500, 15000, 17500))+ #scale_y_log10(breaks=c(1000,3000, 5000, 10000, 15000))+
  scale_fill_manual(values = c("#C77CFF", "#F8766D", "#00BFC4", "#F8766D", "#F8766D")) #scale_fill_manual(values = c("#C77CFF", "#F8766D", "#00BFC4", "#7CAE00", "#F8766D"))

time = ggplot(data,aes(x=data_type))+
   # geom_jitter(aes(y = time,col = sim_type), width = 0.35, size=2.2)  +
  # geom_boxplot(aes(y = time),alpha=0,width = 0.8,outlier.color = NA )+
  geom_boxplot(aes(y=time, fill=sim_type))+
  labs(x="",y="Time (seconds)")+theme_bw() +theme(legend.position = "none")+ facet_wrap("simulator", scale="free", drop=TRUE)+ #scale_y_continuous(breaks=c(60:65,1250,2500,3750,5000,6250,7500, 8750 )) +#scale_y_log10()+
  scale_fill_manual(values = c("#C77CFF", "#F8766D", "#00BFC4", "#F8766D", "#F8766D")) #scale_fill_manual(values = c("#C77CFF", "#F8766D", "#00BFC4", "#7CAE00", "#F8766D"))
burnin= ggarrange(time, mem, ncol = 1, labels = c("A", "B"))
ggexport(burnin, filename="figure_3.pdf", width=7.5, height=7)
ggsave(plot=burnin, filename="burnin_oneset3.jpg", width=7.5, height=7)



div=ggplot(data, aes(x=data_type)) + 
  geom_jitter(aes(y = diversity/expected_div,col = sim_type), width = 0.35, size=2.2)  +
  geom_boxplot(aes(y = diversity/expected_div),alpha=0,width = 0.8,outlier.color = NA )+
   # geom_hline(yintercept = 1, col="grey")+ 
  labs(x="Simulation Type",y="Relative Diversity") +theme_bw()+theme(legend.position = "none")+ facet_wrap("simulator", drop=TRUE) +
  scale_fill_manual(values = c("#C77CFF", "#F8766D", "#00BFC4", "#7CAE00", "#F8766D"))
ggexport(div, filename="burnin_div_nofill.pdf")
ggsave(plot=div, filename="burnin_div.jpg", width=7.5, height=5)


div=ggplot(data, aes(x=data_type)) + 
  geom_boxplot(aes(y = diversity/expected_div, fill=sim_type))+ #coord_cartesian(y=c(0.95,1.03))+
  # geom_hline(yintercept = 1, col="grey")+ 
  labs(x="Simulation Type",y="Relative Diversity") +theme_bw()+theme(legend.position = "none")+ facet_wrap("simulator", drop=TRUE, scale='free_x') +
  scale_fill_manual(values = c("#C77CFF", "#F8766D", "#00BFC4", "#F8766D", "#F8766D")) #scale_fill_manual(values = c("#C77CFF", "#F8766D", "#00BFC4", "#7CAE00", "#F8766D"))
ggexport(div, filename="burnin_div_oneset1.pdf", height = 4, width=8)
ggsave(plot=div, filename="burnin_div_oneset1.jpg", height = 4, width=8)

div=ggplot(data[sim_type!="msprime_ts"], aes(x=data_type)) + 
  geom_boxplot(aes(y = diversity/expected_div, fill=sim_type))+ #coord_cartesian(y=c(0.95,1.03))+
  # geom_hline(yintercept = 1, col="grey")+ 
  labs(x="Simulation Type",y="Relative Diversity") +theme_bw()+theme(legend.position = "none")+ facet_grid(~simulator, drop=TRUE, scale='free_x', space='free_x',labeller=label_wrap_gen(width = 10, multi_line = TRUE)) +
  scale_fill_manual(values = c("#C77CFF", "#00BFC4", "#F8766D", "#F8766D")) #scale_fill_manual(values = c("#C77CFF", "#F8766D", "#00BFC4", "#7CAE00", "#F8766D"))
ggexport(div, filename="figure_4.pdf", height = 4, width=8)
ggsave(plot=div, filename="figure_4.jpg", height = 4, width=8)


mem_msprime = ggplot(data[simulator=="msprime"])+
  geom_boxplot(aes(x=data_type,y=memory, col=sim_type)) + labs(x="",y="Memory Usage (MB)") +theme(legend.position = "none")
# + facet_wrap("simulator", scale="free_y")

time_msprime = ggplot(data[simulator=="msprime"])+
geom_boxplot(aes(x=data_type,y=time, col=sim_type)) + labs(x="",y="Time (seconds)") +theme(legend.position = "none")
burnin_msprime= ggarrange(mem_msprime, time_msprime)
ggexport(burnin_msprime, filename="burnin_msprime.pdf")

mem_slim = ggplot(data[simulator=="SLiM"])+
  geom_boxplot(aes(x=data_type,y=memory, col=sim_type)) + labs(x="",y="Memory Usage (MB)") +theme(legend.position = "none")
# + facet_wrap("simulator", scale="free_y")

time_slim = ggplot(data[simulator=="SLiM"])+
  geom_boxplot(aes(x=data_type,y=time/60, col=sim_type)) + labs(x="",y="Time (minutes)") +theme(legend.position = "none")
burnin_slim= ggarrange(mem_slim, time_slim)
ggexport(burnin_slim, filename="burnin_slim.pdf")


ggplot(data[simulator=="SLiM"])+
  geom_boxplot(aes(x=sim_type,y=memory, col=sim_type)) + labs(y="Memory Usage (MB)") 
# + facet_wrap("simulator", scale="free_y")

mem_all = ggplot(data)+
  geom_boxplot(aes(x=sim_type,y=memory, col=sim_type)) + labs(y="Memory Usage (MB)") + scale_y_log10()
ggexport(burnin_msprime, filename="burnin_msprime.pdf")

time_all = ggplot(data)+
  geom_boxplot(aes(x=sim_type, y=time, col=sim_type))+ labs(y="Time (seconds)") + scale_y_log10()
# + facet_wrap("simulator", scale="free_y")


###  forward benchmarking

for_t_ts=fread("~/phd_data/benchmarking/forward_benchmarking_time_ts.txt")
for_m_m=fread("~/phd_data/benchmarking/forward_mut_benchmarking_time_multi.txt")
for_m_s=fread("~/phd_data/benchmarking/forward_mut_benchmarking_time_single.txt")
forward_time=rbind(for_t_ts,for_m_m,for_m_s)

forward_mem=fread("~/phd_data/benchmarking/forward_benchmarking_mem.txt")

forward_time[, `:=` (mean_time = mean(time), var_time=var(time),  max_time=max(time), min_time=min(time)), by=c("sim_type", "model")]
forward_mem[, `:=` ( mean_memory=mean(memory),  var_mem=var(memory), max_mem=max(memory), min_mem=min(memory)), by=c("sim_type", "model")]
forward_time[, sim_type:=ifelse(sim_type=="mut", "Classical", "Tree Sequence"), by="sim_type"]
forward_time[, model:=ifelse(model=="single_locus", "Single Locus", "Multilocus"), by="model"]
forward_mem[, sim_type:=ifelse(sim_type=="mut", "Classical", "Tree Sequence"), by="sim_type"]
forward_mem[, model:=ifelse(model=="single_locus", "Single Locus", "Multilocus"), by="model"]
forward=merge(unique(forward_mem[, .(sim_type, model, mean_memory, var_mem, max_mem, min_mem)]), unique(forward_time[, .(sim_type, model, mean_time, var_time, max_time, min_time)]), by=c("sim_type", "model"))




plot1 = ggplot(forward_mem)+
  geom_boxplot(aes(x=sim_type, y=memory, fill=sim_type))+
  facet_wrap(~factor(model, levels=c("Single Locus", "Multilocus")))+theme_bw()+ scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="Memory (MB)", x="Simulation Type")+theme(legend.position = "none")+scale_y_log10()
ggexport(plot, filename="forward_mem.pdf")

plot2 = ggplot(forward_time)+
  geom_boxplot(aes(x=sim_type, y=time, fill=sim_type))+theme_bw()+
  facet_wrap(~factor(model, levels=c("Single Locus", "Multilocus")))+ scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  labs(y="Time (seconds)", x="Simulation Type")+scale_y_log10()+theme(legend.position = "none", axis.title.x.bottom = element_blank())
ggexport(plot, filename="forward_time.pdf")

forwardp= ggarrange(plot2, plot1, ncol = 1, labels = c("A", "B"))
ggexport(forwardp, filename="figure_5.pdf", height = 7, width=7.5)
ggsave(plot=forwardp, filename="forwards.jpg", height = 7, width=7)



plot1 = ggplot(forward_mem, aes(x=sim_type))+
  geom_jitter(aes(y = memory,col = sim_type), width = 0.35, size=2.2)  +
  geom_boxplot(aes(y = memory),alpha=0,width = 0.8,outlier.color = NA )+
  facet_wrap("model")+theme_bw()+
  labs(y="Memory (MB)", x="Simulation Type")+theme(legend.position = "none")+scale_y_log10()

plot2 = ggplot(forward_time, aes(x=sim_type))+
  geom_jitter(aes(y = time,col = sim_type), width = 0.35, size=2.2)  +
  geom_boxplot(aes(y = time),alpha=0,width = 0.8,outlier.color = NA )+  theme_bw()+
  facet_wrap("model")+
  labs(y="Time (seconds)", x="Simulation Type")+scale_y_log10()+theme(legend.position = "none", axis.title.x.bottom = element_blank())

forwardp= ggarrange(plot2, plot1, ncol = 1, labels = c("A", "B"))
ggexport(forwardp, filename="forward_no fill.pdf")
### analysis benchmarking

analysis=fread("~/phd_data/benchmarking/analysis_benchmarking.txt")
analysis[, memory:=as.numeric(memory)/1e6] ##memory in MB
# analysis[, stat_type:=ifelse(stat_type=="ts","Tree Sequence", "Allele-based"), by="stat_type"]
analysis[, model:=ifelse(model=="single_locus", "Single Locus", "Multilocus"), by="model"]
analysis[, stat:=ifelse(stat=="div", "Nucleotide Diversity", "Tajima's D"), by="stat"]
analysis[, x_lab:=paste(data_type, "\n", stat_type), by=c("data_type", "stat_type")]

analysis[, `:=` (mean_calc_time = mean(calc_time),mean_total_time = mean(total_time), mean_memory=mean(memory), var_total_time=var(total_time),var_calc_time=var(calc_time), var_mem=var(memory), max_total_time=max(total_time), min_total_time=min(total_time),max_calc_time=max(calc_time), min_calc_time=min(calc_time), max_mem=max(memory), min_mem=min(memory)), by=c("stat_type","data_type", "model", "stat")]
analsysi_sum = unique(analysis[, .(stat_type, data_type, model, stat, mean_memory, var_mem, max_mem, min_mem, mean_calc_time, var_calc_time, min_calc_time, max_calc_time, mean_total_time, var_total_time, min_total_time, max_total_time)])
ggplot(analysis)+
  geom_boxplot(aes(x=stat, y=memory, fill=stat_type))+
  theme_bw()+
  facet_wrap(model~stat_type, scale="free_y")+
  labs(y="Memory (MB)", x="Statistic")+theme(legend.position = "none")


  
plot3.1 = ggplot(analysis[stat_type=="Allele-based" & x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_jitter(aes(y = memory), width = 0.35, size=2.2,colour="#F8766D")  +
  geom_boxplot(aes(y = memory),alpha=0,width = 0.8,outlier.color = NA )+
theme_bw()+
  facet_grid(model~stat)+
  labs(y="Memory (MB)", x="Statistic")+theme(legend.position = "none", axis.title.x.bottom = element_blank(), strip.text.y = element_blank())
plot3.2 = ggplot(analysis[stat_type=="Tree Sequence"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
geom_jitter(aes(y = memory), width = 0.35, size=2.2,colour="#00BFC4")  +
  geom_boxplot(aes(y = memory),alpha=0,width = 0.8,outlier.color = NA )+theme_bw()+
  facet_grid(model~stat)+
  labs(y="Memory (MB)", x="Statistic")+theme(legend.position = "none", axis.title.y.left = element_blank(), axis.title.x.bottom = element_blank())
plot3=ggarrange(plot3.1, plot3.2)
# ggexport(plot, filename="analysis_mem.pdf")
plot4.1 = ggplot(analysis[stat_type=="Allele-based"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_jitter(aes(y = calc_time), width = 0.35, size=2.2,colour="#F8766D")  +
  geom_boxplot(aes(y = calc_time),alpha=0,width = 0.8,outlier.color = NA )+  theme_bw()+
  facet_grid(model~stat)+
  labs(y="Time (seconds)", x="")+theme(legend.position = "none", axis.title.x.bottom = element_blank(), strip.text.y = element_blank())
plot4.2 = ggplot(analysis[stat_type=="Tree Sequence"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_jitter(aes(y = calc_time), width = 0.35, size=2.2,colour="#00BFC4")  +
  geom_boxplot(aes(y = calc_time),alpha=0,width = 0.8,outlier.color = NA )+  theme_bw()+
  facet_grid(model~stat)+
  labs(y="Time (seconds)", x="")+theme(axis.text.x=element_text(angle = -90, hjust = 0),legend.position = "none", axis.title.y.left  = element_blank(), axis.title.x.bottom = element_blank())
plot4=ggarrange(plot4.1, plot4.2)
# ggexport(plot, filename="analysis_time.pdf")

plot3.1 = ggplot(analysis[stat=="Nucleotide Diversity"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_boxplot(aes(x=stat_type, y=memory, fill=stat_type))+ scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  theme_bw()+coord_cartesian(y=c(0, 2.05))+
  facet_grid(model~stat)+
  labs(y="Memory (MB)", x="Statistic")+theme( axis.title.y.left = element_blank(),legend.position = "none", axis.title.x.bottom = element_blank(), strip.text.y = element_blank())
plot3.2 = ggplot(analysis[stat=="Tajima's D"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_boxplot(aes(x=stat_type, y=memory,fill=stat_type))+theme_bw()+ scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  facet_grid(model~stat)+coord_cartesian(y=c(0, 2))+
  labs(x="Statistic")+theme(strip.text.y = element_blank(),legend.position = "none", axis.title.y.left = element_blank(), axis.title.x.bottom = element_blank())
plot3=ggarrange(plot3.1, plot3.2, ncol=1)
plot3=annotate_figure(plot3, left = "Memory (MB)")
# ggexport(plot, filename="analysis_mem.pdf")
plot4.1 = ggplot(analysis[stat=="Nucleotide Diversity"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_boxplot(aes(x=stat_type, y=calc_time,fill=stat_type))+theme_bw()+ scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  facet_grid(model~stat)+coord_cartesian(y=c(0, 0.75))+
  labs(y="Calculation time (seconds)", x="")+theme( axis.title.y.left = element_blank(),legend.position = "none", axis.title.x.bottom = element_blank(), strip.text.y = element_blank())
plot4.2 = ggplot(analysis[stat=="Tajima's D"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_boxplot(aes(x=stat_type, y=calc_time,fill=stat_type))+theme_bw()+ scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  facet_grid(model~stat)+ coord_cartesian(y=c(0, 1.5))+
  labs(y="Calculation time (seconds)", x="")+theme(strip.text.y = element_blank(),legend.position = "none", axis.title.y.left  = element_blank(), axis.title.x.bottom = element_blank())
plot4=ggarrange(plot4.1, plot4.2, ncol=1)
plot4=annotate_figure(plot4, left = "Calculation time (seconds)")
plot5.1 = ggplot(analysis[stat=="Nucleotide Diversity"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_boxplot(aes(x=stat_type, y=total_time,fill=stat_type))+theme_bw()+ scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  facet_grid(model~stat)+coord_cartesian(y=c(0, 1.05))+
  labs(y="Total time (seconds)", x="")+theme( axis.title.y.left = element_blank(),legend.position = "none", axis.title.x.bottom = element_blank())
plot5.2 = ggplot(analysis[stat=="Tajima's D"& x_lab!="Tree Sequence \n Allele-based"], aes(x=stat_type))+
  geom_boxplot(aes(x=stat_type, y=total_time,fill=stat_type))+theme_bw()+ scale_fill_manual(values=c("#00BFC4","#F8766D"))+
  facet_grid(model~stat)+coord_cartesian(y=c(0, 1.45))+
  labs(y="Total time (seconds)", x="")+theme(legend.position = "none", axis.title.y.left  = element_blank(), axis.title.x.bottom = element_blank())
plot5=ggarrange(plot5.1, plot5.2, ncol=1)
plot5=annotate_figure(plot5, left="Total time (seconds)")

ggexport(plot, filename="analysis_time.pdf")

stat_an= ggarrange(plot3,plot4,plot5,  nrow = 1, labels = c("A", "B", "C"))
stat_an=annotate_figure(stat_an, bottom = "Calculation Type")
ggexport(stat_an, filename="figure_6.pdf", height=6, width=9)
ggsave(plot=stat_an, filename="figure6.jpg", height = 6, width=9)


 ## export to xl

write_xlsx(list( forward=forward, analysis=analsysi_sum), path="benchmark_forward_analysis.xlsx",col_names=TRUE)
write_xlsx(list( burnin=averages), path="benchmark_burnin.xlsx",col_names=TRUE)



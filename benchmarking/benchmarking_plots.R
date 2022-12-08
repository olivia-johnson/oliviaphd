## Benchmarking  ##
library(data.table)
library(ggplot2)
library(lubridate)
library(ggpubr)

setwd("~/phd_data/benchmarking/")

##Burn-in
msprime_ts=fread("~/phd_data/benchmarking/burnin_benchmarking_msprime_coalescence.txt")
msprime_ts[, sim_type:="msprime_ts"]
msprime_ts[, simulator:="msprime"]
msprime_ts[, data_type:="Tree Sequence"]

msprime_mut=fread("~/phd_data/benchmarking/burnin_benchmarking_msprime_mutations.txt")
msprime_mut[, sim_type:="msprime_mut"]
msprime_mut[, simulator:="msprime"]
msprime_mut[, data_type:="Neutral Mutations"]


slim_ts=fread("~/phd_data/benchmarking/burnin_benchmarking_slim_ts.txt")
slim_ts[, start_time:=period_to_seconds(hms(start_time))]
slim_ts[,end_time:=period_to_seconds(hms(end_time))]
slim_ts[, time:=end_time-start_time]
slim_ts[, sim_type:="slim_ts"]
slim_ts[, simulator:="SLiM"]
slim_ts[, data_type:="Tree Sequence"]


slim_mut=fread("~/phd_data/benchmarking/burnin_benchmarking_slim_mut.txt")
slim_mut[, start_time:=period_to_seconds(hms(start_time))]
slim_mut[,end_time:=period_to_seconds(hms(end_time))]
slim_mut[, time:=end_time-start_time]
slim_mut[, sim_type:="slim_mut"]
slim_mut[, simulator:="SLiM"]
slim_mut[, data_type:="Neutral Mutations"]



data=rbind(msprime_ts,msprime_mut,slim_ts[,.(time, memory, sim_type, simulator, data_type)], slim_mut[,.(time, memory, sim_type, simulator, data_type)])
data[, memory:=memory/1e6]  ## convert memory to MB

data[, `:=` (mean_time = mean(time), mean_memory=mean(memory), sd_time=sd(time), sd_mem=sd(memory), max_time=max(time), min_time=min(time), max_mem=max(memory), min_mem=min(memory)), by="sim_type"]

averages = unique(data[,.(sim_type, mean_time, mean_memory, sd_time, sd_mem, max_mem, min_mem, max_time, min_time)])


##plots##
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

for_t_m=fread("~/phd_data/benchmarking/forward_benchmarking_time_multi.txt")
for_t_s=fread("~/phd_data/benchmarking/forward_benchmarking_time_single.txt")
forward_time=rbind(for_t_m,for_t_s)

forward_mem=fread("~/phd_data/benchmarking/forward_benchmarking_mem.txt")
plot = ggplot(forward_mem)+
  geom_boxplot(aes(x=model, y=memory))+
  # facet_wrap(model, scale="free_y")+
  labs(y="Memory (MB)", x="Statistic")
ggexport(plot, filename="forward_mem.pdf")

plot = ggplot(forward_time)+
  geom_boxplot(aes(x=model, y=time))+
  # facet_wrap(model, scale="free_y")+
  labs(y="Time (s)", x="Statistic")
ggexport(plot, filename="forward_time.pdf")


### analysis benchmarking

analysis=fread("~/phd_data/benchmarking/analysis_benchmarking.txt")
analysis[, memory:=as.numeric(memory)/1e6] ##memory in MB
plot = ggplot(analysis)+
  geom_boxplot(aes(x=stat, y=memory))+
  facet_wrap(model~stat_type, scale="free_y")+
  labs(y="Memory (MB)", x="Statistic")
ggexport(plot, filename="analysis_mem.pdf")
plot = ggplot(analysis)+
  geom_boxplot(aes(x=stat, y=time))+
  facet_wrap(model~stat_type, scale="free_y")+
  labs(y="Time (s)", x="Statistic")
ggexport(plot, filename="analysis_time.pdf")


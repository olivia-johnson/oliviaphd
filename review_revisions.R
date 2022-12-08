library(data.table)
library(ggpubr)
library(viridis)
library(tidyverse)
library(ggdist)
library(readxl)
library(grDevices)
library(cowplot)
library(intervals)
setwd("~/phd_data/drosData")



rudman=as.data.table(read_excel("rudman/science.abj7484_tables_s1_s3_and_s6.xlsx", 
                                sheet = "Table S6"))
rudman[, c("chrom", "region") := tstrsplit(cluster_region, ":", fixed=TRUE)]
rudman[, c("cluster_start", "cluster_end") := tstrsplit(region, "-", fixed=TRUE)]
setnames(rudman, "marker_SNP_position", "pos")
rudman=rudman[,.(`cluster label`, cluster_start, cluster_end, chrom, pos, `time segment`)]
rudman[, snp_id:=paste0(chrom, "_", pos)]
rudman[, `:=` (cluster_start=as.numeric(cluster_start), cluster_end=as.numeric(cluster_end))]


data.labs = c('chrom', 'pos','AF', 'SP', 'SQ' ,'FallF', 'SprF')
bergland_snps= fread("bergland/bergland_PA.txt", col.names = data.labs)
bergland_snps[, snp_id:=paste0(chrom, "_", pos)]

load("machado/mel_freqdp_042016_Ne_fixed_correctBAVI.Rdata")
machado_data=as.data.table(cbind(info, freq))

m_snps<- as.data.table(read_excel("machado/elife-67577-supp1-v2.xlsx", 
                                  sheet = "Supplementaryfile1B"))
m_pops <-as.data.table(read_excel("machado/elife-67577-supp1-v2.xlsx", 
                                  sheet = "Supplementaryfile1A"))
m_snps[, snp_id:=paste0(chrom, "_", pos)]
## check overlap between rudman, machado, and bergland
intersect(rudman$snp_id, bergland_snps$snp_id)
intersect(rudman$snp_id, m_snps$snp_id)
intersect(bergland_snps$snp_id, m_snps$snp_id)


overlap_clusters=NULL
for (i in 1:length(unique(rudman$chrom))){
  chr=unique(rudman$chrom)[i]
  interval_data=data.matrix(rudman[chrom==chr, .(cluster_start, cluster_end)])
  int_data=Intervals(interval_data)
  int_overlap=interval_union(int_data)
  new_cluster=as.data.table(int_overlap)
  setnames(new_cluster, c("V1", "V2"), c("start", "end"))
  new_cluster[, chrom:=chr]
  overlap_clusters=rbind(overlap_clusters, new_cluster)
}
overlap_clusters$cluster_id <- 1:nrow(overlap_clusters)

intervals::plot(int_overlap,xlim=c(0,as.numeric(dros[chrom==chr, length])), ylim=c(0,15))



intervals::plot(int_overlap,xlim=c(0,as.numeric(dros[chrom==chr, length])), ylim=c(0,15))

 

overlap_snps=NULL
for (i in 1:length(bergland_snps$chrom)){
  chr=bergland_snps$chrom[i]
  position=bergland_snps$pos[i]
  within=overlap_clusters[chrom==chr & start<=position & end>=position]
  within[,`:=` (snp_position=position, data="bergland")]
  overlap_snps=rbind(overlap_snps, within)
}
for (i in 1:length(m_snps$chrom)){
  chr=m_snps$chrom[i]
  position=m_snps$pos[i]
  within=overlap_clusters[chrom==chr & start<=position & end>=position]
  within[,`:=` (snp_position=position, data="machado")]
  overlap_snps=rbind(overlap_snps, within)
}

m_cluster=unique(overlap_snps[data=="machado", .(cluster_id)])
b_cluster = unique(overlap_snps[data=="bergland", .(cluster_id)])

intersect(m_cluster, b_cluster)




## testing

overlap_clusters[, length:=end-start]


dros_chrom=c("X", "2L","2R","3L","3R")
dros_chromlen=c(22422827,23011544,21146708,24543557,27905053)  ## from https://www.fruitfly.org/sequence/README.RELEASE5
dros_cumm_start=c(1,22422828,45434372,66581080,91124637)
dros_cumm_end=c(22422827,45434371,66581079,91124636,119029689)
dros=data.table(dros_chrom, dros_chromlen, dros_cumm_start, dros_cumm_end)
setnames(dros, c("dros_chrom","dros_chromlen", "dros_cumm_start", "dros_cumm_end"), c("chrom", "length", "s_start", "s_end"))


for (i in 1:length(unique(dros$chrom))){
  chr=unique(dros$chrom)[i]
  overlap_clusters[chrom==chr, `:=` (c_start=start+dros[chrom==chr, s_start]-1, c_end=end+dros[chrom==chr, s_start]-1), by="cluster_id"]
  m_snps[chrom==chr, `:=` (s_pos=pos+dros[chrom==chr, s_start]-1, s_end=s_pos+1), by="snp_id"]
  bergland_snps[chrom==chr, `:=` (s_pos=pos+dros[chrom==chr, s_start]-1, s_end=s_pos+1), by="snp_id"]
  }

overlap_clusters[, `:=` (length=end-start, c_length=c_end-c_start)]
sum(overlap_clusters$c_length)
sum(overlap_clusters$length)
sum(overlap_clusters$length)/119029689

rudman_cl=overlap_clusters[,.(c_start, c_end)]
berg_snps=Intervals(data.matrix(bergland_snps[,.(s_pos, s_end)]))
mac_snps=Intervals(data.matrix(m_snps[,.(s_pos, s_end)]))


overlap=data.table(iteration=integer(), b_overlap=integer(), m_overlap=integer(), rotation=integer())
it=0
for (i in c(0, runif(10, 1, 119029689))){
  it=it+1
  rot_val=rudman_cl[, .(start=c_start+i, end=c_end+i)]
  rot_cl=Intervals(data.matrix(rot_val))
  b_ov=length(interval_difference(berg_snps, rot_cl))
  m_ov=length(interval_included(mac_snps, rot_cl))
  overlap=rbind(overlap, data.table(iteration=it, b_overlap=b_ov, m_overlap=m_ov, rotation=i))
}












dros_over=data.table(chrom=character(), start=numeric(), end=numeric(), length=numeric())
for (i in 1:length(overlap_clusters$cluster_id)){
  c.len=overlap_clusters[cluster_id==i, length]
  d.chr=sample(dros_chrom, size = 1)
  d.len=dros[chrom==d.chr, length]
  d.start=runif(1, min=1, max=d.len)
  d.end=d.start+c.len
  intervals=Intervals(data.matrix(dros_over[chrom==d.chr, .(start,end)]))
  new_int=Intervals(cbind(d.start,d.end))
  overlap=data.table(interval_intersection(intervals, new_int))
  while(d.end>d.len | length(overlap$V1)>0){
    d.chr=sample(dros_chrom, size = 1)
    d.len=dros[chrom==d.chr, length]
    d.start=runif(1, min=1, max=d.len)
    d.end=d.start+c.len
    intervals=Intervals(data.matrix(dros_over[chrom==d.chr, .(start,end)]))
    new_int=Intervals(cbind(d.start,d.end))
    overlap=data.table(interval_intersection(intervals, new_int))
  }
  dros_over=rbind(dros_over,data.table(chrom=d.chr, start=d.start,end=d.end, length=d.len))
}

check_clusters=NULL
for (i in 1:length(unique(dros_over$chrom))){
  chr=unique(dros_over$chrom)[i]
  interval_data=data.matrix(dros_over[chrom==chr, .(start, end)])
  int_data=Intervals(interval_data)
  new_cluster=as.data.table(interval_union(int_data))
  setnames(new_cluster, c("V1", "V2"), c("start", "end"))
  new_cluster[, chrom:=chr]
  check_clusters=rbind(check_clusters, new_cluster)
  
}
check_clusters

intervals::plot(Intervals(data.matrix(dros_over[chrom==chr, .(start, end)])),xlim=c(0,as.numeric(dros[chrom==chr, length])), ylim=c(0,5))


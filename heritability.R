library(dplyr)

groups = c(29:56)

ped_names=c("Gen", "Ind", "P1", "P2", "G1", "G2", "G3", "G4")

heritability=NULL

for (g in groups){
  print(g)
 
  parameters <- fread(file=paste0("~/phd_data/seglift_allele/group_",g, "/parameters.txt"), sep = ":")
  setkey(parameters, V1)

  fiton=parameters["fitness_on", V2]
  if (fiton==0){
    dom=0
    epi=0
  } else{
    dom = parameters["d", V2]
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

  g_label = ifelse(fiton==0, paste(pop_s, gen_s, "No Fitness", loci, sep="_"), paste(pop_s, gen_s,dom, epi, loci, sep="_"))

  
  
  pedigree  <- fread(file=paste0("~/phd_data/seglift_allele/group_", g,  "/pedigree_",g, "_0.txt"), col.names = ped_names)
  P1=as.data.table(table(pedigree$P1))
  P2= as.data.table(table(pedigree$P2))
  G1= as.data.table(table(pedigree$G1))
  G2= as.data.table(table(pedigree$G2))
  G3= as.data.table(table(pedigree$G3))
  G4= as.data.table(table(pedigree$G4))
  
  P=rbindlist(list(P1, P2))[, lapply(.SD, sum, na.rm = TRUE), by = V1]
  P=P[,V1:=as.integer(V1)]
  setnames(P, "V1", "Ind")
  setnames(P, "N", "Parent")
  G=rbindlist(list(G1, G2, G3, G4))[, lapply(.SD, sum, na.rm = TRUE), by = V1]
  G=G[,V1:=as.integer(V1)]
  setnames(G, "V1", "Ind")
  setnames(G, "N", "Grandparent")
  output=full_join(x = P, y=G, by="Ind")
  output[is.na(output)] <- 0
  pedigree = merge(pedigree, output, by="Ind", all=TRUE)
  pedigree[is.na(Parent),Parent:=0]
  pedigree[is.na(Grandparent),Grandparent:=0]
  pedigree[, label:=g_label]
  pedigree[, group:=g]
  pedigree[,gen_season:=gen_s]
  pedigree[,pop_season:=pop_s]
  heritability=rbind(heritability,pedigree)
}

load
heritability=na.omit(heritability, cols = c("Gen"))
heritability[, slope:=Grandparent/Parent]
heritability[,h2:=2*slope]


## assign each generation a time point in the year
heritability[gen_season == "EG", gen_year:=Gen%%30, by=c("label", "Gen")]
heritability[gen_season == "UG", gen_year:=Gen%%15, by=c("label", "Gen")] 

heritability[,m_parent:=mean(Parent), by=c("label", "gen_year")]
heritability[,m_grandparent:=mean(Grandparent), by=c("label", "gen_year")]
heritability[,g_h2:=m_grandparent/m_parent, by=c("label", "Gen")]
heritability[,m_h2:=m_grandparent/m_parent, by=c("label", "gen_year")]

heritability[gen_season == "EG", season:=ifelse((gen_year<16 & gen_year>0), "summer", "winter"), by=c("Gen")]
heritability[gen_season == "UG", season:=ifelse((gen_year<14 & gen_year>0), "summer", "winter"), by=c("Gen")]

her = unique(heritability[, .(label, m_h2, gen_year)])

ggplot(heritability[])+
  geom_line(aes(x=Gen, y=g_h2))+
  facet_wrap("label")


ggplot(heritability)+
  geom_violin(aes(x=as.factor(gen_year), y=h2), nrow=4)+
  facet_wrap("label")
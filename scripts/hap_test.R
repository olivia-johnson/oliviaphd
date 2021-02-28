library(data.table)
library(stringr)

group = 0
run = 1

hap_list <- list.files(path = "~/oliviaphd/seglift_even_dist/slim_out/", pattern=paste0("ms_out_", group, "_", run, "_"))


for (i in 1:(length(hap_list))){
  file_name = hap_list[i]
  output = fread(file = paste0("~/oliviaphd/seglift_even_dist/slim_out/",file_name))
  output = output[-(1:2),.N, by="//"]
  setnames(output, "//", "hap")
  setnames(output, "N", "hap_count")
  output[, n_sel:= str_count(hap, "1")]
  haps = output[, .N, by=n_sel]
  haps$run = run
}



for (j in (unique(output$time))){
  hapout = fread(file = paste0("~/oliviaphd/seglift_even_dist/slim_out/","ms_out_", group, "_", i, "_", j, ".txt"))
  hapout = hapout[-(1:2),.N, by="//"]
  setnames(hapout, "//", "hap")
  setnames(hapout, "N", "hap_count")
  hapout[, n_sel:= str_count(hap, "1")]
  haps = hapout[, .N, by=n_sel]
  haps$run = i
  haps$gen = j
  hapout$run = i
  hapout$gen = j
  haplotypes = rbind(haplotypes, hapout)
  hap_stats = rbind(haps, hap_stats)
}
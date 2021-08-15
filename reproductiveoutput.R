## Reproductive output
library(data.table)
library(ggplot2)

g=12 

output = fread(paste0("~/oliviaphd/reproductiveOutput_12_1.txt"), sep=",", fill=TRUE) ##fitness off, CP=10000, EG = 15

# output[is.na(output)] <- 0
# 
# colnames(output)<- c("Gen", 0:(length(output)-2))

offout=ggplot(data=output[Gen==7], aes(x=CountofReproductiveOutput))+
  geom_histogram(binwidth=1)+
  geom_line(aes(y=dpois(0:10000, 1), colour="red"))

            
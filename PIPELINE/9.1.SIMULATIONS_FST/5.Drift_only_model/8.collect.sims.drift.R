##### Evaluate FST per SNPs

library(tidyverse)
library(magrittr)
library(poolfstat)

library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(gtools)

library(foreach)
library(doMC)
library(doParallel)
registerDoMC(4)

###### user defined parameters
args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1])
print(k)
### here k is 1-100
### and we are collecting each simulation

####### Load reco_wins
reco_wins <- fread("reco_wins.master.txt")
names(reco_wins) = c("chr", "reco_level", "start", "end", "lenght" , "rho")

### collect FST
fstfiles <- system("ls fst_sims_all", intern = T)

o1 <- foreach(fil = fstfiles, .combine = "rbind")%do%{
  
  message(fil)
  tmp <- get(load(paste("fst_sims_all/",fil, sep = "")))
  tmp %>% 
    filter(ith == k)
  
}

o1$pos = as.numeric(o1$pos)

o1 %>%
  group_by(chr, pos, samp2) %>%
  arrange(pos) -> o1.arr

save(o1.arr,
     file = paste("/gpfs2/scratch/jcnunez/fst_brent/simulations/collect_FST/collect.FST.sim.drift.",k,".Rdata", sep = "")
     )

### collect counts
countfiles <- system("ls simulated_pools_counts", intern = T)

o2 <- foreach(fil = countfiles, .combine = "rbind")%do%{
  
  message(fil)
  tmp <- get(load(paste("simulated_pools_counts/",fil, sep = "")))
  tmp %>% 
    filter(ith == k)
  
}

o2 %>%
  group_by(chr, real_pos) %>%
  arrange(real_pos) -> o2.arr

save(o2.arr,
     file = paste("/gpfs2/scratch/jcnunez/fst_brent/simulations/collect_COUNTS/collect.COUNTS.sim.drift.",k,".Rdata", sep = "")
)


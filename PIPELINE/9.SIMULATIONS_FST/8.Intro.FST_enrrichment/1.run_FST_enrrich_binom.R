##### FST ENRICHMENT.
#### 

library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)

fstfiles = system("ls /gpfs2/scratch/jcnunez/fst_brent/simulations_redo/collect_FST_intro/", intern = T)

###### user defined parameters
args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1])
print(k)
   
i = fstfiles[k]

    inner=get(load(
      paste("/gpfs2/scratch/jcnunez/fst_brent/simulations_redo/collect_FST_intro/", i, sep = "")
    ))
    
    inner %>%
      filter(samp2 == "VT8") ->
      inner
    
    ###
    L = dim(inner)[1] 
    
    inner %>%
      arrange(-snp.fst) %>% 
      as.data.frame() %>%
      mutate(rank = seq(from=1,  to = L, by = 1)) %>% 
      mutate(rank_norm = rank/L) %>% 
      arrange(as.numeric(pos))->
      inner.rnf
    
    p.i=0.01
    inner.rnf %>%
      group_by(winSTA,winEND,chr) %>% 
      summarise(N = sum(rank_norm < p.i),
                Lw = n()) %>%
      mutate(p.binom = c(binom.test(N, Lw, p.i)$p.value)) %>%
      mutate(p.binom.adj = p.adjust(p.binom)) %>%
      mutate(samp1="sim",
             samp2= "VT8",
             comp="simulation"
      ) %>% mutate(ith = strsplit(i, "\\.")[[1]][5],
                   type = "simulation") ->
      inner.binomtest

save(inner.binomtest,
     file = paste("/gpfs2/scratch/jcnunez/fst_brent/simulations/collect_enrrichFSTBinom/", 
                  "FSTenrichBinom.", strsplit(i, "\\.")[[1]][5], ".Rdata", sep = "" )
     )
    
    
    
    
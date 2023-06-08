#### Plot FST
#### 

library(tidyverse)
library(magrittr)
library(foreach)

###
args = commandArgs(trailingOnly=TRUE)
samp1.s = args[1]
samp2.s = args[2]
comp_name = args[3]
chr = args[4]

###
#samp1.s="VT8_0_1"
#samp2.s="VT8F_1_1"
#comp_name="VT8_vs_reps"
#chr = "2L"
###

fst.out <- system(paste("ls ./fst_out/* | grep 'Rdata' | grep",chr, sep = " " ), intern = T)

ag.fst.df =
  foreach(i = 1:length(fst.out), .combine = "rbind")%do%{
    
    print(i)
    tmp <- get(load(fst.out[i]))
    ####
    tmp %>%
    filter(!is.na(snp.fst)) %>%
    filter(samp1 == samp1.s) %>%
    filter(samp2 == samp2.s) %>%   
    mutate(comp = comp_name) -> tmp
    
    return(tmp)

  }# open foreach i

save(ag.fst.df, 
     file = 
       paste(paste(samp1.s,samp2.s,chr,sep = "."), "ag.fst.df.Rdata", sep = ".")
     )


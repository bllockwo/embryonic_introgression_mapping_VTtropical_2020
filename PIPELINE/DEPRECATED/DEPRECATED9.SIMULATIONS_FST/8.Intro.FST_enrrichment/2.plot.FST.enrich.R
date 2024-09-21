#### Plot enrichment of FST
#### 

library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)

####
#### Load sims
sim_enr = system("ls collect_enrrichFSTBinom", intern = T)

sim_data = 
  foreach(i=sim_enr, .combine = "rbind")%do%{
    
    message(i)
    tmp <- get(load(paste("collect_enrrichFSTBinom/", i, sep = "")))
    
  }
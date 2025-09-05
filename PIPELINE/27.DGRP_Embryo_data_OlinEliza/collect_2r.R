library(tidyverse)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(adegenet)
library(reshape2)
library(FactoMineR)
require(gtools)
require(foreach)

files = 
  system("ls /gpfs2/scratch/jcnunez/fst_brent/Fall_2025_revisions/Exploration_DGRP/fix_2R/set.*.Rdata", 
         intern = T)

o = foreach(i=files,
            .combine = "rbind")%do%{
              tmp <- get(load(i))
            }

o %>%
  filter(type != "empirical_p") %>%
  ggplot(
    aes(
      x=as.numeric(pos),
      y=-log10(uci_1),
      color=type
    )
  ) + geom_vline(xintercept = 15501637) + geom_line()  ->
  testplot
ggsave(testplot, file = "testplot.pdf")

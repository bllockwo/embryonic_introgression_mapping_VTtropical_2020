## PLOT PCA
#Load libraries
library(tidyverse)
library(vcfR)
library(data.table)
library(adegenet)
library(FactoMineR)
library(magrittr)
library(reshape2)
library(zoo)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(ggrepel)
library(lubridate)
library(foreach)
library(SeqArray)
library(doMC)
registerDoMC(8)
library(fastglm)
library(lmtest )
library(gmodels)
library(lme4)
library(foreach)

fils <- system("ls pc_Corr", intern = T)

inday <- foreach(i=fils, .combine = "rbind")%do%{
  tmp <- get(load(paste("pc_Corr",i, sep = "/")))
}


inday %>%
  filter(pc %in% 1:2) %>%
  filter(type == "M") %>%
  ggplot(aes(
    x=log10(window_size),
    y=-log10(P),
    color=as.character(pc)
  )) + geom_line() +
  geom_hline(yintercept = -log10(0.05)) +
  facet_grid(pc~geno) + theme_bw() ->
plot1
ggsave(plot1, file = "plotcor.pdf", w = 6, h = 3)

###

inday %>%
  filter(pc %in% 1:2) %>%
  filter(type == "M") %>%
  group_by(pc, geno) %>%
  slice_min(P)

inday %>%
  filter(pc %in% 1:2) %>%
  filter(type == "M") %>%
  filter(window_size == 1600)


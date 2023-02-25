#### Plot FST
#### 

library(tidyverse)
library(magrittr)
library(foreach)

###

fst.out <- system("ls ./fst_out/* | grep 'Rdata'  ", intern = T)

intro.fst.df =
foreach(i = 1:length(fst.out), .combine = "rbind")%do%{
  
  print(i)
  tmp <- get(load(fst.out[i]))
  return(tmp)
  
  }# open foreach i

head(intro.fst.df)
tail(intro.fst.df)

###
###
intro.fst.df %>%
  mutate(comp.type = case_when(
    samp1 == "VT8_0_1" & samp2 %in% c("VT8F_1_1", "VT8F2_2_1", "VT8F3_3_1") ~ "VT8_group",
    samp1 == "VT10_0_1" & samp2 %in% c("VT10F_1_1", "VT10F2_2_2", "VT10F3_3_1") ~ "VT10_group",
  )) %>%
  filter(!is.na(comp.type)) %>% 
  mutate(mid.Mbp = (as.numeric(winSTA+winEND)/2)/1e6) %>% 
  ggplot(
    aes(x= mid.Mbp,
        y= FST,
        color = samp2)
  ) +
    geom_line() +
    facet_grid(comp.type*samp2~chr, scales = "free_x") +
    theme(legend.pos = "top")->
    fst.plots.comp
  ggsave(fst.plots.comp, file = "fst.plots.comp.pdf", w = 9, h = 7)
#### 
####
####
####
####

  intro.fst.df %>%
    mutate(comp.type = case_when(
      samp1 == "VT8_0_1" & samp2 %in% c("VT8F_1_1", "VT8F2_2_1", "VT8F3_3_1") ~ "VT8_vs_reps",
      samp1 == "VT10_0_1" & samp2 %in% c("VT10F_1_1", "VT10F2_2_2", "VT10F3_3_1") ~ "VT10_vs_reps",
      samp1 == "VT8_0_1" & samp2 %in% c("SKF_1_1", "SKF2_2_1", "SKF3_3_1") ~ "SK_vs_VT8",
      samp1 == "VT10_0_1" & samp2 %in% c("CHF_1_1", "CHF2_2_1", "CHF3_3_2") ~ "CH_vs_VT10",
    )) %>%
    filter(!is.na(comp.type)) %>% 
    mutate(mid.Mbp = (as.numeric(winSTA+winEND)/2)/1e6) %>% 
    ggplot(
      aes(x= mid.Mbp,
          y= FST,
          color = samp2)
    ) +
    geom_line() +
    facet_grid(comp.type~chr, scales = "free") +
    theme(legend.pos = "top")->
    fst.plots.comp
  ggsave(fst.plots.comp, file = "fst.plots.comp.pdf", w = 9, h = 7)
  
  
  ###
  intro.fst.df %>%
    filter(samp2 == "VT8_0_1" & samp1 == "VT10_0_1") %>% 
    mutate(mid.Mbp = (as.numeric(winSTA+winEND)/2)/1e6) %>% 
    ggplot(
      aes(x= mid.Mbp,
          y= FST)
    ) +
    geom_line() +
    facet_grid(.~chr, scales = "free") +
    theme(legend.pos = "top")->
    fst.plots.VT.comp
  ggsave(fst.plots.VT.comp, file = "fst.plots.VT.comp.pdf", w = 9, h = 3)
  
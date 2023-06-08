#### Plot FST
#### 

library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)

#####
#####
met = vroom("master.file.fst.txt")
#####
#####

outer=
foreach(i = 1:13, .combine = "rbind")%do%{
  message(paste(i,met$samp1[i],met$samp2[i]))
inner=
foreach(k = c("2L","2R","3L","3R","X"), .combine = "rbind")%do%{
  fil <- paste(met$samp1[i],met$samp2[i],k,"ag.fst.df.Rdata",sep =".")
  tmp <- get(load(fil)) %>% 
    separate(pos, into = c("chr.o","pos"), sep = "_") %>%
    arrange(as.numeric(pos))
  return(tmp)
}
L = dim(inner)[1] 
inner %>%
  arrange(-snp.fst) %>% 
  mutate(rank = 1:L) %>% 
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
  mutate(samp1=met$samp1[i],
         samp2= met$samp2[i],
         comp=met$comp[i]
           )->
  inner.binomtest
return(inner.binomtest)
}


#####
#####

outer %>%
  ggplot(aes(
    x=((winSTA+winEND)/2)/1e6,
    y=-log10(p.binom.adj),
    group = paste(samp1,samp2),
    color = paste(samp1,samp2)
  )) +
  geom_vline(xintercept = 15391154/1e6) +
  geom_vline(xintercept = 20276334/1e6) +
  geom_line() +
  ggtitle("FST enrrich top 1%") +
  facet_grid(paste(samp1,samp2)~chr, scales = "free_x") ->
  testP
ggsave(testP, 
       file = "fst.errich.pdf", 
       w = 8, h = 9)

###save(ag.fst.df, file = "ag.fst.df.Rdata")

head(ag.fst.df)
tail(ag.fst.df)

###
###
ag.fst.df %>%
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
  
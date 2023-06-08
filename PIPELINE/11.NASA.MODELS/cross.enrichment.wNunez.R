
rm(list = ls())

### libraries
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(magrittr)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)



left_join(tmp.0, 
          intro.fet[,c("chr","pos","sk_all_f.test_pval","ch_all_f.test_pval")]) %>%
  arrange(chr, pos) %>%
  filter(!is.na(sk_all_f.test_pval)) %>%
  filter(!is.na(ch_all_f.test_pval)) ->
  tmp.joint

tmp.joint <- tmp.joint[chr!="X"]
tmp.joint %<>% 
  mutate(L = dim(.)[1]) %>%
  mutate( GLM.rank = rank(p_lrt, ties.method = "first")/L) %>%
  mutate( SK.rank = rank(sk_all_f.test_pval, ties.method = "first")/L) %>%
  mutate( CH.rank = rank(ch_all_f.test_pval, ties.method = "first")/L) %>%
  group_by(chr) %>% arrange(pos)
setDT(tmp.joint)

#### === > Make windows
win.bp <- 1e5
step.bp <- 5e4
setkey(tmp.joint, "chr")

## prepare windows
wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  tmp <- tmp.joint %>%
                    filter(chr == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]
dim(wins)
#### === > Make windows

thrs <- c(0.01)
k = 0

cross.enrich.out = 
  foreach(pos.ith=1:dim(wins)[1], 
          .combine="rbind", 
          .errorhandling="remove")%do%{
            
            message( paste0(unlist(wins[pos.ith]), sep = "/") )
            
            tmp.joint %>% filter(chr== wins$chr[pos.ith]) %>%
              filter(pos > wins$start[pos.ith] &
                       pos < wins$end[pos.ith]) ->
              tmp.win
            
            tab <- table(tmp.win$GLM.rank  < thrs,
                         tmp.win$SK.rank <  thrs)
            
            fet <- fisher.test(tab)
            
            tmp.out <- data.table(chr=wins$chr[pos.ith], 
                                  thr=thrs , 
                                  perm=k,
                                  anchor.model=mod,
                                  win.start= wins$start[pos.ith],
                                  win.end= wins$end[pos.ith],
                                  or=fet$estimate, 
                                  p=fet$p.value, 
                                  lci=fet$conf.int[1], 
                                  uci=fet$conf.int[2],
                                  Nsnps = dim(tmp.win)[1]
            )
            
            ### return
            return(tmp.out)
            
            
          } ### close pos.ith

#} ### close models


ggplot() +
  geom_line(
    data = cross.enrich.out[chr == "2R"], aes(
      x=win.start,
      y=-log10(p)
    )) +
  # geom_point(
  #   data = cross.enrich.out[p < 0.05][chr == "2R"],
  #   aes(
  #     x=win.start,
  #     y=or)) +
  ##geom_hline(yintercept = -log10(0.05)) +
  facet_wrap(~chr, scales = "free_x")->
  tes

ggsave(tes, file ="tes.pdf")


#### Plot fst results 

#Load libraries
library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(ggrepel)
library(foreach)

### Population level FST
root <- "/gpfs2/scratch/jcnunez/fst_brent/slice_VCFs/results_fst_pops"
files = system(paste ("ls", root), intern = T)

o <- foreach(i=files, 
             .combine = "rbind",
             .errorhandling = "remove")%do%{
               
               message(i)
               info=str_split(str_split(i, "\\.")[[1]][1], "_")[[1]]
               #extract info
               stat=info[1]
               pop=str_split(info[2], "-")[[1]][1]
               set=str_split(info[2], "-")[[1]][2]
               W=info[4]
               S=info[6]
               #paste(stat,pop,set,W,S)
               
               dat<-fread(paste(root,i, sep = "/"))%>%
                 mutate(pop=pop)
             }

o %>%
  mutate(fix_fst = case_when(WEIGHTED_FST < 0 ~ 0,
                             TRUE ~ WEIGHTED_FST)) %>%
  ggplot(aes(
    x=(BIN_START+BIN_END)/2e6,
    y=fix_fst,
    color=pop
  ))    + geom_vline(xintercept = 20551633/1e6,
                     color = "black") +
  geom_line() + 
  theme_bw() +
  facet_grid(pop~.)->
  plot

ggsave(plot, file = paste("fst",".plot.pdf",sep = ""),
       w= 4, h =5) 


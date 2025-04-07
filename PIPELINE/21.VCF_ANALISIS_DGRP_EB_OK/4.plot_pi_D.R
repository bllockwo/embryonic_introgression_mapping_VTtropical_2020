#### Plot Pi/D results 

#Load libraries
library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(ggrepel)
library(foreach)

### find outputs
root <- "/gpfs2/scratch/jcnunez/fst_brent/slice_VCFs/results_pi_D_1Mb"
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
  
  dat<-fread(paste(root,i, sep = "/"))
  
  if(stat == "D"){
  data.frame(chr=dat$CHROM,
             N=dat$N_SNPS,
             value=dat$TajimaD) %>%
  mutate(stat=stat,
             pop=pop,
             set=set,
             W=W,
             S=S
             ) -> in1
    if(pop=="DGRP"){
      in1 %<>% mutate(pos_m = dat$BIN_START + 4548375)
    } else { in1 %<>% mutate(pos_m = dat$BIN_START )}
  }
  
  if(stat == "Pi"){
    data.frame(chr=dat$CHROM,
               N=dat$N_VARIANTS,
               value=dat$PI) %>%
      mutate(stat=stat,
             pop=pop,
             set=set,
             W=W,
             S=S
      ) -> in1
    if(pop=="DGRP"){
      #pos_mean = (dat$BIN_START+dat$BIN_END)/2
      in1 %<>% mutate(pos_m = dat$BIN_START + 4548375)
    } else { in1 %<>% mutate(#pos_m = (dat$BIN_START+dat$BIN_END)/2 
                              pos_m = dat$BIN_START 
                             )}
  }
  return(in1)
  }



for(i in c("Pi","D")){
  o %>%
    filter(stat == i) %>%
    ggplot(aes(
      x=pop,
      y=value,
      color=set
    )) + geom_boxplot() +
    scale_color_manual(values = c("steelblue", "grey2", "firebrick"))+
    theme_bw()->
    plot
  
  if(i == "D"){
    plot + geom_hline(yintercept = 0, 
                      linetype = "dashed", linewidth = 0.52) -> plot
  }
  if(i == "Pi"){
    plot + geom_hline(yintercept = 0.005, 
                      linetype = "dashed", linewidth = 0.52) -> plot
  }
  
  ggsave(plot, file = paste(i,"box.plot.pdf",sep = ""),
         w= 3, h =3) 
}



for(i in c("Pi","D")){
  o %>%
    filter(pop != "DGRP") %>%
    filter(stat == i) %>%
    ggplot(aes(
      x=pos_m/1e6,
      y=value,
      color=set
    )) + geom_line() + 
    geom_vline(xintercept = 20551633/1e6,
               color = "red") +
    scale_color_manual(values = c("steelblue", "grey2", "firebrick"))+
    facet_grid(pop~stat) ->
    plot
  
  if(i == "D"){
    plot + geom_hline(yintercept = 0, 
                      linetype = "dashed", linewidth = 0.52) -> plot
  }
  if(i == "Pi"){
    plot + geom_hline(yintercept = 0.005, 
                      linetype = "dashed", linewidth = 0.52) -> plot
  }
  
  ggsave(plot, file = paste(i,".plot.pdf",sep = ""),
         w= 3, h =5) 
}


library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)

fstfiles = system("ls /gpfs2/scratch/jcnunez/fst_brent/simulations/collect_FST_intro/", intern = T)

#args = commandArgs(trailingOnly=TRUE)
#k=as.numeric(args[1])
#print(k)


fst.wins = 
foreach(k=1:100, .combine = "rbind")%do%{
  i = fstfiles[k]
  
  message(k)
  inner=get(load(
    paste("/gpfs2/scratch/jcnunez/fst_brent/simulations/collect_FST_intro/", i, sep = "")
  ))
  
  inner %>% 
    filter(samp2 == "VT8") %>%
    group_by(chr, winSTA, winEND) %>%
    summarise(winSTA = mean(winSTA),
              winEND = mean(winEND),
              fstm = mean(FST)
    ) %>%
    mutate(ith = k, type = "sim") -> o
  
  return(o)
  
}

save(fst.wins, file = "fst.wins.mean.simsIntro.Rdata")


###  
met = vroom("/gpfs2/scratch/jcnunez/fst_brent/master.file.fst.txt")
met[grep("VT8",met$samp1),] -> met

real_fst=
  foreach(i = 1:dim(met)[1], .combine = "rbind")%do%{
    message(paste(i,met$samp1[i],met$samp2[i]))
    inner=
      foreach(k = c("2L","2R","3L","3R","X"), .combine = "rbind")%do%{
        fil <- paste(met$samp1[i],met$samp2[i],k,"ag.fst.df.Rdata",sep =".")
        tmp <- get(load(paste("../",fil, sep = "") )) %>% 
          separate(pos, into = c("chr.o","pos"), sep = "_") %>%
          arrange(as.numeric(pos)) %>%
          group_by(chr, winSTA, winEND) %>%
          summarise(winSTA = mean(winSTA),
                    winEND = mean(winEND),
                    fstm = mean(FST)
          ) %>%
          mutate(ith = 0, type = "real")
        
        return(tmp)
      }
    }
save(real_fst, file = "real_fst.mean.Rdata")


load("fst.wins.mean.simsIntro.Rdata")
load("real_fst.mean.Rdata")


fst.wins %>%
  group_by(chr,winSTA,winEND) %>%
  summarize(
    sdfst = sd(fstm, na.rm = T),
    fst95 = quantile(fstm, 0.90),
    fst05 = quantile(fstm, 0.10),
    fstm = mean(fstm, na.rm = T),
            ) %>%
rbind(., real_fst) %>%
  ggplot(aes(
    x=( winSTA   + winEND)/2e6,
    y=fstm,
    ymin=fstm-2.576*(sdfst/sqrt(100)),
    ymax=fstm+2.576*(sdfst/sqrt(100)),
    #ymin=fst05,
    #ymax=fst95,
    color = type,
    #group=ith
  )) +
  geom_point() +
  geom_ribbon(alpha = 0.8) +
  geom_smooth() +
  theme_bw() +
  facet_wrap(~chr, nrow = 1, scales = "free_x") ->
  sim.fst.plot

ggsave(sim.fst.plot, file = "sim.fst.plot.intro.pdf", h = 3, w = 9)

#####
foreach(i = 1:dim(real_fst)[1], .combine = "rbind")%do%{
  message(i)
  
  tmp = real_fst[i,]
  tmp.sim = fst.wins %>% filter(chr == tmp$chr) %>%
                         filter(winSTA > tmp$winSTA & winEND < tmp$winEND)
  
  100-length(tmp.sim$fstm > tmp$fstm)
                                

}


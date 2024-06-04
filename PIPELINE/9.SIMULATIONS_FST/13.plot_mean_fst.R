library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)

fstfiles = system("ls /gpfs2/scratch/jcnunez/fst_brent/simulations/collect_FST/", intern = T)

args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1])
print(k)


fst.wins = 
foreach(k=1:100, .combine = "rbind")%do%{
  i = fstfiles[k]
  
  message(k)
  inner=get(load(
    paste("/gpfs2/scratch/jcnunez/fst_brent/simulations/collect_FST/", i, sep = "")
  ))
  
  inner %>% 
    filter(samp2 == "VT8") %>%
    group_by(chr, winSTA, winEND) %>%
    summarise(winSTA = mean(winSTA),
              winEND = mean(winEND),
              fstm = mean(FST)
    ) %>%
    mutate(ith = k, type = "sim")
  
}

save(fst.wins, file = "fst.wins.mean.sims.Rdata")


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

    
rbind(real_fst, fst.wins) %>%
  ggplot(aes(
    x=( winSTA   + winEND)/2,
    y=fstm,
    color = type,
    group=ith
  )) +
  geom_line() +
  facet_wrap(~chr) ->
  sim.fst.plot

ggsave(sim.fst.plot, file = "sim.fst.plot.pdf")

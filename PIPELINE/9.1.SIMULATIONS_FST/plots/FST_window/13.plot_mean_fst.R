library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)
library(data.table)

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
load("fst.wins.mean.simsIntro.Rdata")

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
          group_by(chr, winSTA, winEND, samp1, samp2 ) %>%
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
    fst95 = quantile(fstm, 0.95),
    fst05 = quantile(fstm, 0.05),
    fstm = mean(fstm, na.rm = T),
            ) ->
  sim.summaries

real_fst$fstm[real_fst$fstm < 0] = 0

ggplot()+
  geom_line(data=real_fst,
            aes(
              x=( winSTA   + winEND)/2e6,
              y=fstm,
              color = samp2
            ) 
  ) +
  geom_point(data=real_fst,
             aes(
               x=( winSTA   + winEND)/2e6,
               y=fstm,
               color = samp2
             ) 
  ) +
  geom_ribbon(data=sim.summaries,
              alpha = 0.8,linewidth = 1.2,
              aes(
                x=( winSTA   + winEND)/2e6,
                y=fstm,
                ymin=0,
                ymax=fstm+2.576*(sdfst/sqrt(100)),
              )
  )  +
  facet_wrap(.~chr, nrow = 1, scales = "free_x") +
  scale_color_brewer(palette = "Accent") +
  theme_minimal() + theme(legend.position = "bottom")->
  sim.fst.plot

ggsave(sim.fst.plot, file = "sim.fst.plot.intro.pdf", h = 4, w = 9)

###### ---> Characterize the patterns of FST
real_fst %>%
  dplyr::select(chr, winSTA, winEND, samp2, fstm ) ->
  real_fst_snapshot
setDT(real_fst_snapshot)

sim.summaries %>%
  dplyr::select(chr,winSTA,winEND, fst95) %>%
  group_by(chr) %>%
  arrange(chr, winSTA) ->
  sim_fst_snapshot
setDT(sim_fst_snapshot)

setkey(sim_fst_snapshot, chr, winSTA,  winEND)

foverlaps(real_fst_snapshot, 
          sim_fst_snapshot, 
          by.x=c("chr", "winSTA", "winEND"),
          type="any") ->
  overlap_dt
overlap_dt$fst95[overlap_dt$fst95 < 0] = 0

overlap_dt %>%
  group_by(chr, i.winSTA, i.winEND, samp2 ) %>%
  summarise(mean.sim = mean(fst95),
            mean.fstm = mean(fstm )
  ) %>%
  ungroup() %>%
  group_by(chr, i.winSTA, i.winEND) %>%
  summarise(sig.test = sum(mean.fstm > mean.sim+0.01)) %>%
  filter(sig.test == 6) ->
  allwins.beat

dim(allwins.beat)
allwins.beat$chr %>% table

allwins.beat %>%
   ggplot(aes(
    x=(i.winSTA+i.winEND)/2,
    y=sig.test
  )) +
  geom_point(shape = 15, size = 3) +  
  theme_minimal() +  
  facet_wrap(.~chr, nrow = 1, scales = "free_x") ->
  count.fst.beatsim

ggsave(count.fst.beatsim, file = "count.fst.beatsim.pdf", h = 4, w = 9)



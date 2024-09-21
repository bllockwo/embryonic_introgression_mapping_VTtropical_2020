#### Plot enrichment of FST
#### 

library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)

####
#### Load sims
sim_enr = system("ls collect_enrrichFSTBinom", intern = T)

sim_data = 
  foreach(i=sim_enr, .combine = "rbind")%do%{
    
    message(i)
    tmp <- get(load(paste("collect_enrrichFSTBinom/", i, sep = "")))
    
  }

save("sim_data", file = "enrich.sim.data.Rdata")

sim_data %>%
  group_by(winSTA,  winEND, chr) %>%
  summarise(sim95 = quantile(p.binom.adj, 0.95)) ->
  sim.95.envelope

real_enr = get(load("/gpfs2/scratch/jcnunez/fst_brent/fst.binomial.test.VT8.Rdata"))

ggplot() +
  geom_line(
    data=real_enr,
    aes(
    x=(winSTA+winEND)/2,
    y=-log10(p.binom.adj),
    color=samp2
  )) +
  geom_line(
    data=sim.95.envelope,
    aes(
      x=(winSTA+winEND)/2,
      y=-log10(sim95)
    ), color = "grey", linewidth = 1.9) + 
  theme_minimal() +  
  facet_grid(samp2~chr, scales = "free_x") ->
  real_enrich.plot

ggsave(real_enrich.plot, 
       file = "real_enrich.plot.pdf", h = 6, w = 9)





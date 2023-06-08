### Select Best Models of 2R

library(ggplot2)
library(data.table)
library(patchwork)
library(ggrepel)
library(tidyverse)
library(forcats)

load("/netfiles/nunezlab/Drosophila_resources/Datasets/Nunez2023_Inversion.Seasonal/o2.globalOmnibus.Rdata")

sets <- data.table(mod=c(-1, 0, 1:11),
                   label=LETTERS[1:13],
                   year=c(NA, rep(1, 12)),
                   start=c(NA, NA, 0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(NA, NA, 7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))
####
####

left_join(o2.ag, sets) %>% 
  mutate(model.name = paste(var, start, end, sep = "_")) ->
  o2.ag.annot

#####

o2.ag.annot %>%
  group_by(cluster) %>%
  filter(cluster %in% c(#"1.Europe_W", 
                        #"3.Europe_E", 
                        "5.Cville") ) %>%
  filter(chr == "2R")  ->
  best.mods.2R.cvil

####
####
####
best.mods.2R.cvil %>% arrange(-rr) %>% filter(!var %in% c("precip.var","precip.ave")) %>% filter(sig == "TRUE" & inv == "Inside Inversion") -> Enrrich_models_2R
####
####
  ggplot() +
    geom_linerange(data = Enrrich_models_2R, 
                   aes(x=fct_reorder(model.name, rr),
                       ymin=(prop.real/prop.perm.uci), 
                       ymax=(prop.real/prop.perm.lci),
                       color = cluster), position = position_dodge(width = 0.5), 
                   alpha = 0.9, linewidth = 1.2) +
    geom_point(data = Enrrich_models_2R, 
               aes(x=fct_reorder(model.name, rr),
                   y = rr,
                   fill = cluster),
               #position=position_dodge2(width=.5), 
               size=3, alpha = 1.0, shape = 21,
               position = position_dodge(width = 0.5)) +
    facet_grid(inv~.) +
    geom_hline(yintercept=1, linetype = "dashed") +
    coord_flip() ->
    t1
  ggsave(t1, file = "t1.pdf", w = 6, h = 3.5)
   
    

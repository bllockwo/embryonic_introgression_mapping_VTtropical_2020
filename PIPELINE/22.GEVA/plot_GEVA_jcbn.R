### GEVA age estimate

library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(ggrepel)
library(foreach)
library(forcats)

markers = fread("/netfiles/nunezlab/D_melanogaster_resources/Inversions/inversion_makers_kapun.dm5.txt")
names(markers)[2] = "chr"

markers %>%
  filter(inversion == "In(2R)Ns") ->
  markers.2R
  


geva_dat <- fread("/netfiles/nunezlab/Shared_Resources/in_transit/lproud/VA.cm_GEVA_complete_named.txt")

geva_dat %>%
  separate(id,
           remove = F,
           into = c("chr", "start", "end" ),
           sep = "\\.") ->
  geva_dat.ed
setDT(geva_dat.ed)

geva_dat.ed %>%
  filter(chr == "2R") %>%
  filter(position %in% markers.2R$position) %>%
  summarise(m.TMRCA = median(TMRCA)*0.06666667,
            )

geva_dat.ed %>%
  filter(chr == "2R") %>%
  filter(position %in% 20551633) %>%
  summarise(m.TMRCA = median(TMRCA)*0.06666667,
  )

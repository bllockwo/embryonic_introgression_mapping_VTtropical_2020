library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)
library(foreach)
library(data.table)
library(gmodels)

bestMods="/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.Nunez_et_al_Supergene_paper/Revision_Best_Models/top10models.Rdata"

darMods = get(load(bestMods))

darMods %>%
  filter(chr == "2R") %>% as.data.frame() %>% 
  arrange(prop.rr) %>%
  filter(cluster == "5.Cville") ->
  bets_models2R

### Extract models 


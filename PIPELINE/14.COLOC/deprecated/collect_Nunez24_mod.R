#### Collect models
library(foreach)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
mod.i= as.numeric(args[1])
variable.i=args[2]
cluster.i=args[3]

#mod.i=8
#variable.i="temp.max" 
#cluster.i="5.Cville"

files_loc <- "/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.Nunez_et_al_Supergene_paper/Raw_Output_GLMv2"

files <- system(paste("ls", files_loc), intern = T) 

model_collected =
foreach(i=files, .combine = "rbind",
        .errorhandling = "remove")%do%{
  tmp <- get(load(paste(files_loc, i, sep = "/") ))
  
  message(i)
  
  tmp %>%
    filter(
      cluster == cluster.i) %>%
    filter(
      mod == mod.i) %>%
    filter(
      variable == variable.i) ->
    tmp.flt
}

mod.name = paste(variable.i,mod.i,cluster.i, sep = "_" )

save(model_collected, file = paste(mod.name, "Rdata", sep = "."))


#### Cross-Enrichment 
#### 

library(data.table)
library(tidyverse)
library(magrittr)
library(foreach)
library(doMC)
registerDoMC(4)
library(tidyr)

####

models=
  c(
    "temp.var;5;5.Cville",
    "temp.var;1;5.Cville",
    "temp.propmin;10;5.Cville",
    "temp.propMax;10;5.Cville",
    "temp.propMax;6;5.Cville",
    "temp.min;6;5.Cville",
    "temp.min;3;5.Cville"
  )



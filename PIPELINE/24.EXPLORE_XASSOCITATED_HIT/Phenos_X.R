#### Explore X-associated SNP

library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)
library(forcats)
library(data.table)
library(gmodels)


#### Phenos ....
#### Phenos ....#### Phenos ....
#### Phenos ....
#### Phenos ....
#### Phenos ....
#### Phenos ....
xhist <- fread("Xsnp_info.txt")

wolb = fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/wolbachia.status.txt")
names(wolb)[1] = "sampleid"
E_O_data <- fread("../slice_VCFs/Eliza_Olins_data.txt")
names(E_O_data)[1] = "sampleid"

load("../slice_VCFs/homozyg_SP70.Rdata")
phenos <- readRDS("../slice_VCFs/wideform.fixed.phenotable.RDS")
#names(phenos)

MCT = "HighThermalToleranceExtreme_VaryingWithTemperature_F"
FCT = "HighThermalToleranceExtreme_VaryingWithTemperature_M"

data.frame(
  sampleid = paste("line",phenos$ral_id, sep = "_"),
  CT_F = phenos[,"HighThermalToleranceExtreme_VaryingWithTemperature_F"],
  CT_M = phenos[,"HighThermalToleranceExtreme_VaryingWithTemperature_M"]
) -> phenos_ct
names(phenos_ct) = c("sampleid", "CTF", "CTM")

left_join(phenos_ct, homozyg_SP70) %>%
  left_join(wolb) %>%
  full_join(E_O_data[,c("sampleid","Proportion")]) %>%
  left_join(xhist)->
  ALL_PHENO_DATA_DGRP_SP70_X

ALL_PHENO_DATA_DGRP_SP70_X %>%
  filter(Infection_Status == "n") %>%
  filter(!is.na(SP70)) %>%
  filter(Xsnp != "./.") %>%
  ggplot(aes(
    y=Proportion,
    x=CTF
  )) + geom_point() + geom_smooth(method = "lm") +
  facet_wrap(~Xsnp) ->
  X2R_box
ggsave(X2R_box, file = "X2R_box.pdf")


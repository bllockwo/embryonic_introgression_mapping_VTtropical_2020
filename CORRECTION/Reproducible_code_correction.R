## reproducible code
library(plyr); library(dplyr)
library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)

#setwd(".../GitHub/embryonic_introgression_mapping_VTtropical_2020/CORRECTION/")

dat <- fread("DGRP.ProperlyRelabeled.data.txt")
## Generate a joint genotype label for Sp70 and sog; called "geno"
dat %<>%
  mutate(geno = paste("SP70=",SP70_allele,"-","sog=",SOG_allele, 
                      sep =""))

#####
# Comnbined Genotypes
dat %>%
  ggplot(aes(
    x=geno,
    y=Proportion
  )) + geom_boxplot() +
  ggtitle("Embryonic Heat Tolerance in the DGRP", 
          subtitle = "SP70 and sog genotype combinations") ->
  results1
ggsave(results1, file = "Embryonic.Means.Corrected.pdf")

# Just SP70
dat %>%
  ggplot(aes(
    x=SP70_allele,
    y=Proportion
  )) + geom_boxplot() +
  ggtitle("Embryonic Heat Tolerance in the DGRP", 
          subtitle = "SP70 only") ->
  results2
ggsave(results2, file = "Embryonic.Means.Corrected.Sp70.pdf")

# Just sog
dat %>%
  ggplot(aes(
    x=SOG_allele,
    y=Proportion
  )) + geom_boxplot() +
  ggtitle("Embryonic Heat Tolerance in the DGRP", 
          subtitle = "sog only") ->
  results3
ggsave(results3, file = "Embryonic.Means.Corrected.Sog.pdf")

### Adult Phenotypes
### Load data from Nunez 24
adult_phens <- readRDS("wideform.fixed.phenotable.RDS")
### Extract the CTmax data 
adult_phens %>%
  select("ral_id", 
         "HighThermalToleranceExtreme_VaryingWithTemperature_F",                            
         "HighThermalToleranceExtreme_VaryingWithTemperature_M") %>%
  mutate(line = paste("DGRP", ral_id, sep = "_"))  %>%
  reshape2::melt(id = c("line", "ral_id") ) %>%
  separate(variable, into = c("pheno","variance","sex"), sep = "_") ->
  ctmax_phens
## ... and merge with embryonic data
left_join(dat, ctmax_phens) -> dat_w_adults

## plot both genes
dat_w_adults %>%
  ggplot(aes(
    x=geno,
    y=value,
    fill = sex
  )) + geom_boxplot() +
  ggtitle("Adult Heat Tolerance in the DGRP", 
          subtitle = "SP70 and sog genotype combinations") ->
  results4
ggsave(results4, file = "Adult.Means.Corrected.pdf")

## plot just sp70
dat_w_adults %>%
  ggplot(aes(
    x=SP70_allele,
    y=value,
    fill = sex
  )) + geom_boxplot() +
  ggtitle("Adult Heat Tolerance in the DGRP", 
          subtitle = "SP70 genotype") + theme_bw() ->
  results4
ggsave(results4, file = "Adult.Means.Corrected.Sp70.pdf")

## plot just sog
dat_w_adults %>%
  ggplot(aes(
    x=SOG_allele,
    y=value,
    fill = sex
  )) + geom_boxplot() +
  ggtitle("Adult Heat Tolerance in the DGRP", 
          subtitle = "sog genotype") + theme_bw() ->
  results5
ggsave(results5, file = "Adult.Means.Corrected.sog.pdf")


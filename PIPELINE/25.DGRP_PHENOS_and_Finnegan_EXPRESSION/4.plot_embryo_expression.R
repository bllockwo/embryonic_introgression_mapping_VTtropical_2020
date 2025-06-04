library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)
library(forcats)
library(data.table)
library(gmodels)

###load data
dat <- fread("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/25.DGRP_PHENOS_and_Finnegan_EXPRESSION/Top_genes_expression_summary.txt")

dat %>%
  ggplot(aes(
    x=treatment, 
    y=count,
    ymin=count-se,
    ymax=count+se,
    color=genome
  )) +
  geom_errorbar(width = 0.1,position=position_dodge(width=0.5)) + 
  geom_point(position=position_dodge(width=0.5)) +
  facet_wrap(~factor(gene_name , levels = c("SP70", "Cpr57A", "sog")), scales = "free") +
  scale_color_brewer(palette = "Dark2") + theme_bw()
## save 3*9


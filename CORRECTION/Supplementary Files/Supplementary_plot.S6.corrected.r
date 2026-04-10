library(tidyverse)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(adegenet)
library(reshape2)
library(FactoMineR)
require(gtools)
require(foreach)

genofile <- seqOpen("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/dgrp2.gds.gz")
seqResetFilter(genofile)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      alleles=seqGetData(genofile, "allele"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile)) %>%
  mutate(snp_id = paste(chr,pos, sep = "_"))

snps.dt %>% 
  filter(snp_id %in%
           c(
             "2R_16439138",
             "X_15496974",
             "X_15501637",
             "X_15741847")) -> snps.sel

seqSetFilter(genofile,  
             variant.id=snps.sel$variant.id)                      

allelic.data <- data.table(
  chr=seqGetData(genofile, "chromosome"),
  pos=seqGetData(genofile, "position"),
  variant.id=seqGetData(genofile, "variant.id"),
  nAlleles=seqNumAllele(genofile),
  missing=seqMissing(genofile),
  alleles=seqGetData(genofile, "allele")
) %>%
  mutate(snp_id = paste(chr,pos, sep = "_"))


snpgdsGetGeno(genofile, snp.id=snps.sel$variant.id,
              verbose=TRUE, with.id=TRUE) ->
  GENO.matrix


GENO.matrix$genotype %>%
  as.data.frame() ->
  GENO.matrix.dt

colnames(GENO.matrix.dt) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), sep="_")
rownames(GENO.matrix.dt) <- seqGetData(genofile, "sample.id")
### load data
egg_data <- fread("DGRP.ProperlyRelabeled.data.txt")
egg_data %<>% separate(line, into = c("DGRP","line_id")) %>%
  mutate(Line = paste("line", line_id, sep = "_"))

GENO.matrix.dt %>%
  mutate(Line = rownames(.)) %>%
  full_join(egg_data) %>%
  filter(!is.na(Proportion)) %>%
  select(`2R_16439138`, "X_15496974", "X_15501637", "X_15741847",
         Proportion) %>%
  .[complete.cases(.),] %>%
  reshape2::melt(id = c("2R_16439138", "Proportion")) %>%
  mutate(Joint_Geno = paste(`2R_16439138`, value, sep = "_")) %>%
  ggplot(aes(
    x=Joint_Geno,
    y=Proportion
  )) + geom_boxplot() +
  facet_grid(~variable) + theme_classic() +
  ylim(0, 0.8)->
  Supp.Fig.S6.Corrected

ggsave(Supp.Fig.S6.Corrected, file = "Supp.Fig.S6.Corrected.pdf",
       h=2.5, w = 6)

#### Just Genotypes  
GENO.matrix.dt %>%
  mutate(Line = rownames(.)) %>%
  full_join(egg_data) %>%
  filter(!is.na(Proportion)) %>%
  select(`2R_16439138`, "X_15496974", "X_15501637", "X_15741847",
         Proportion) %>%
  .[complete.cases(.),] %>%
  reshape2::melt(id = c("2R_16439138", "Proportion")) %>%
  mutate(Joint_Geno = paste(`2R_16439138`, value, sep = "_")) %>%
  ggplot(aes(
    x=as.character(value),
    y=Proportion
  )) + geom_boxplot() +
  facet_grid(~variable) + theme_classic() +
  ylim(0, 0.8)->
  Supp.Fig.S6.Corrected.ABC

ggsave(Supp.Fig.S6.Corrected.ABC, file = "Supp.Fig.S6.Corrected.ABC.pdf",
       h=2.5, w = 6)

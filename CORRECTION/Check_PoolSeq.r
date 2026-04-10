### Check of pool-seq

library(tidyverse)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(adegenet)
library(reshape2)
library(FactoMineR)
library(gtools)
library(foreach)
library(ape)

wolbachia_st <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/wolbachia.status.txt")
phenos <- readRDS("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2024.Nunez_et_al_Genetics/Phenotyping/wideform.fixed.phenotable.RDS")
DGRP <- "/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/dgrp2.gds.gz"
invs <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/Inversion.status.txt")
POOLS <- "/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds"

genofile <- seqOpen(POOLS, allow.duplicate=T)
DGRP <- seqOpen(DGRP, allow.duplicate=T)


### RESRET
### RESRET
seqResetFilter(genofile)
### RESRET
### RESRET
### 
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      alleles=seqGetData(genofile, "allele"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))

snps.dt %<>% mutate(SNP_id = paste(chr,pos, sep = "_"))

#SNP_id %in% c(
#"2R_16439138",
#"X_15501637"
#"2R_20551633",
#"X_15607604"

snps.dt %>%
  filter(chr == "2R") %>%
  filter(nAlleles == 2) %>%
  filter(pos == 20552016) -> flt_2R
  #filter(pos > 20551633-100 & pos < 20551633+100) 

snps.dt %>%
  filter(chr == "X") %>%
  filter(nAlleles == 2) %>%
  filter(pos > 15607604-100 & pos < 15607604+100) -> flt_X

#flt2L
rbind(flt_2R
      #,flt_X
      ) -> filt.pools

seqSetFilter(genofile, 
             variant.id=filt.pools$variant.id,
             sample.id= c("SK_0_1", "SKF_1_1", "SKF3_3_1", "SKF2_2_1",
  "VT8_0_1", "VT8F_1_1", "VT8F2_2_1", "VT8F3_3_1"))

### get allele frequency data
ad <- seqGetData(genofile, "annotation/format/AD")$data
dp <- seqGetData(genofile, "annotation/format/DP")

### divide the matrices to generate allele frequency calls
dat <- ad/dp

### check the dimensions of the allele frequency matrix
dim(dat)  

####
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")
#colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")
#colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")

rownames(dat) <- seqGetData(genofile, "sample.id")
#rownames(ad) <- seqGetData(genofile, "sample.id")
#rownames(dp) <- seqGetData(genofile, "sample.id")
pdf("nj.pdf")
plot(nj(dist(dat)))
dev.off()

### Geno Step
snpgdsGetGeno(genofile, 
              verbose=TRUE, with.id=TRUE,
              snp.id =flt_i$variant.id) ->
  GENO.matrix

MATRIX <- GENO.matrix$genotype

#colnames(MATRIX)  <- c()
rownames(MATRIX)  <- GENO.matrix$sample.id

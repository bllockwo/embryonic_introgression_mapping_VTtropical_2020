### EXTR. DATA
### 

library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(data.table)
library(matrixStats)
library(foreach)
library(doParallel)
library(doMC)
registerDoMC(4)

library(SeqArray)
library(gdsfmt)
library(SNPRelate)

library(DescTools)

###
samps <- fread("/netfiles02/lockwood_lab/IntrogressionProject/Population_files/samp.guide.files.introg.txt")
samps %>% 
  filter(Pop %in% c("SK","VT8",
                    "SKF","VT8F",
                    "SKF2","VT8F2",
                    "SKF3","VT8F3",
                    "ME_temp","PAN_trop"
  )) -> samps.flt


### Load the genotype locality
genofile.path <- "/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds"

### load the genofile into memory
genofile <- seqOpen(genofile.path)

#### Load the filtering SNP object -- which JCBN created
filtering.dt <- get(load("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/SNP.filtering.guide.Rdata"))

filtering.dt %>%
  filter(is.na(libs)) %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) ->
  chosen.snps
####

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))


snps.dt %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) %>%
  right_join(chosen.snps, by = c("chr", "pos", "SNP_id") ) ->
  snps.dt

###
###
###

snps.dt <- snps.dt[nAlleles==2]

seqSetFilter(genofile, sample.id = samps.flt$Sample_id, 
             variant.id=snps.dt$variant.id)

snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

seqSetFilter(genofile,
             snps.dt[chr%in%c("2L", "2R", "3L", "3R", "X")][missing<.05]$variant.id)


#### extract matrices of cov and calls
### get allele frequency data
ad <- seqGetData(genofile, "annotation/format/AD")$data
dp <- seqGetData(genofile, "annotation/format/DP")

### divide the matrices to generate allele frequency calls
dat <- ad/dp

### check the dimensions of the allele frequency matrix
dim(dat)  

####
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")
colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")

rownames(dat) <- seqGetData(genofile, "sample.id")
rownames(ad) <- seqGetData(genofile, "sample.id")
rownames(dp) <- seqGetData(genofile, "sample.id")

####
dat %>%
  t() %>% 
  as.data.frame -> dat_t

### row means
rowMeans(dat_t) -> row_menas_cov
dat_t[which(row_menas_cov > 0),] -> dat_t.nofix

dat_t.nofix %>%
  mutate(SNP_id = rownames(.)) %>%
  separate(SNP_id, into = c("chr", "pos", "feature"), sep = "_") ->
  dat_t.nofix.meta

dat_t.nofix.meta %>% .$chr %>% table


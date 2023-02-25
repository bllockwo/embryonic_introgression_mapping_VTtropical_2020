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

seqSetFilter(genofile, variant.id=snps.dt$variant.id)

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

ad %>%
  t() %>% 
  as.data.frame -> ad_t

dp %>%
  t() %>% 
  as.data.frame -> dp_t

#### you are done and now may want to save these objects for future use
####

##### Example script
high.bulk = "CHF2_2_1"
low.bulk = "VT8_0_1"

ad_t %>%
  .[,c(high.bulk,low.bulk)] %>%
  mutate(SNP_id = rownames(.)) ->
  ad.sel
names(ad.sel) = c("hb_ad", "lb_ad", "SNP_id")

dp_t %>%
  .[,c(high.bulk,low.bulk)] %>%
  mutate(SNP_id = rownames(.)) ->
  dp.sel
names(dp.sel) = c("hb_dp", "lb_dp", "SNP_id")

left_join(ad.sel, dp.sel) %>% 
  dplyr::select(SNP_id, hb_ad, lb_ad, hb_dp, lb_dp) %>%
  mutate(hb_ref = hb_dp-hb_ad,
         lb_ref = lb_dp-lb_ad
         ) ->
  tmp.dt

tmp.dt %>% head %>% dplyr::select(!c(hb_dp,lb_dp) ) -> tmp.dt.h

tmp.matrx <-
  #matrix(c(3, 1, 1, 3),
  matrix( c(tmp.dt$hb_ref[i], tmp.dt$hb_ad[i], tmp.dt$lb_ref[i],  tmp.dt$lb_ad[i])  ,
         nrow = 2
         #dimnames = list(High = c("REF", "ALT"),
         #                 Low =  c("REF", "ALT"))
         )

GTest(tmp.matrx)



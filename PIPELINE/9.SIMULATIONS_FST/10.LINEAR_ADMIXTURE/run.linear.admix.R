### Linear admixture

### libraries
library(SeqArray)
library(data.table)
library(foreach)
library(tidyverse)
library(magrittr)
library(vroom)

library(SeqArray)
library(gdsfmt)
library(SNPRelate)


###

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

dat_t.nofix %<>% 
  mutate(SNP_id = rownames(.)) %>% 
  separate(SNP_id, into =c("chr","pos","feature"), sep = "_")

####
children <- c( "SKF_1_1","VT8F_1_1",
               "SKF2_2_1","VT8F2_2_1",
               "SKF3_3_1","VT8F3_3_1")

reco_wins = fread("/gpfs2/scratch/jcnunez/fst_brent/simulations/reco_wins.master.txt")
names(reco_wins) = c("chr", "cM","start","end", "L", "rho")

global_admixture = 
  foreach(i = 1:dim(reco_wins)[1], 
          .combine = "rbind", 
          .errorhandling = "remove")%do%{
  oin=          
  foreach(au.i=children, 
               .combine = "rbind", 
               .errorhandling = "remove")%do%{
                 
                 message(au.i)
                 test_pop = au.i 
                 #"SK","VT8"
                 dat_t.nofix %>% 
                   filter(chr == reco_wins$chr[i]) %>%
                   filter(pos > reco_wins$star[i] &  pos < reco_wins$end[i]) -> dat_t.nofix.ch
                 
                 dat_t.nofix.ch[sample(dim(dat_t.nofix.ch)[1],1000),
                             c("SK_0_1","VT8_0_1", test_pop)] -> 
                   dat.tmp
                 
                 setnames(dat.tmp, names(dat.tmp), c("SK", "VT", "TEST"))
                 
                 summary(lm(TEST~0+SK+VT, dat.tmp)) -> mod.out
                 
                 as.data.frame(mod.out$coefficients) %>%
                   mutate(sampleId = test_pop,
                          chr = reco_wins$chr[i],
                          star = reco_wins$star[i],
                          end = reco_wins$end[i]
                   ) %>%
                   mutate(source_pop = rownames(.)) ->
                   o
                 
                 return(o)
               }# inner loop
          return(oin)
          }#outer loop
save(global_admixture, file = "global_admixture.linear.Rdata")

global_admixture$star = as.numeric(global_admixture$star)
global_admixture$end = as.numeric(global_admixture$end)

global_admixture %>%
  filter(source_pop == "SK") %>%
  ggplot(aes(
    x=(star+end)/2,
    y=Estimate,
    ymin=Estimate-1.96*`Std. Error`,
    ymax=Estimate+1.96*`Std. Error`,
    color=sampleId
  )) +
  geom_ribbon(alpha = 0.1) + geom_line() +
  ylab("Estimated Genomic Ancestry") +
  facet_wrap(~chr, ncol = 1, scales = "free_x") -> admixlin

ggsave(admixlin, file = "admixlin.pdf", h = 9, w = 5)



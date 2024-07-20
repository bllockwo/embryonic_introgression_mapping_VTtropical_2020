##### Build PCA with simulations

library(plyr)
library(tidyverse)
library(magrittr)
library(reshape2)
library(SeqArray)
library(FactoMineR)
library(factoextra)
library(data.table)
library(foreach)

#####
genofile.path <- "/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds"
genofile <- seqOpen(genofile.path)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))

snps.dt <- snps.dt[nAlleles==2]
seqSetFilter(genofile, variant.id=snps.dt$variant.id)
snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

seqSetFilter(genofile,
             snps.dt[chr%in%c("2L", "2R", "3L", "3R", "X")][missing<.05]$variant.id)

### get allele frequency data
ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")

dat <- ad$data/dp
dim(dat)  

colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")

#### Now bring in simulated data
intro_data = system("ls /gpfs2/scratch/jcnunez/fst_brent/simulations/AFS_intro_objects/*", intern = T)
Drift_data = system("ls /gpfs2/scratch/jcnunez/fst_brent/simulations/AFS_drift_objects/*", intern = T)

guide_set <- get(load(intro_data[1]))
guide_set_drif <- get(load(Drift_data[1]))

intro_snps <- colnames(guide_set)
drift_snps <- colnames(guide_set_drif)

shared_snps <- intro_snps[which(intro_snps %in% drift_snps)]

snps_valid <- shared_snps[which(shared_snps %in% colnames(dat))]

intro_list = list()
for(i in 1:length(intro_data)){
  message(i)
  tmp = get(load(intro_data[i]))
  intro_list[[i]] = tmp[,snps_valid]
}

do.call(rbind, intro_list) -> Intro.AFS.valid

save(Intro.AFS.valid, file = "intro.AF.data.Rdata")

###

drift_list = list()
for(i in 1:length(Drift_data)){
  message(i)
  tmp = get(load(Drift_data[i]))
  drift_list[[i]] = tmp[,snps_valid]
}

do.call(rbind, drift_list) -> Drift.AFS.valid

save(Drift.AFS.valid, file = "Drift.AFS.valid.Rdata")

##### Evaluate FST per SNPs

library(tidyverse)
library(magrittr)
library(poolfstat)

library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(gtools)

library(foreach)
library(doMC)
library(doParallel)
registerDoMC(4)

###### user defined parameters
args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1])
print(k)

####### Load reco_wins
reco_wins <- fread("reco_wins.master.txt")
names(reco_wins) = c("chr", "reco_level", "start", "end", "lenght" , "rho")

###
set <- paste(reco_wins$chr[k], reco_wins$start[k], reco_wins$end[k], "drift", sep = "_")

### load sims 
fileset = paste(set,".fst.all.Rdata", sep ="")

sim_fst <- get(load( paste("/gpfs2/scratch/jcnunez/fst_brent/simulations/fst_sims_all", 
                          fileset, 
                          sep = "/") ))

#### load real data
##### objects needed for the analyses
samps = fread("/netfiles02/lockwood_lab/IntrogressionProject/Population_files/samp.guide.files.introg.txt")
genofile.path <- "/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds"
filtering.dt <- get(load("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/SNP.filtering.guide.Rdata"))

### load the genofile into memory
genofile <- seqOpen(genofile.path)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile),
                      allele=seqGetData(genofile, "allele")) 
snps.dt <- snps.dt[nAlleles==2]

####
snps.dt %>%
  filter(chr == reco_wins$chr[k] & pos > reco_wins$start[k] & pos < reco_wins$end[k] ) %>% 
  mutate(reco_rho = reco_wins$reco_level[k]) ->
  snps.dt.win


###
seqResetFilter(genofile)
seqSetFilter(genofile,
             sample.id=samps$Sample_id[c(grep( "VT8" , samps$Sample_id), 
                                         grep( "SK" , samps$Sample_id))],
             variant.id=
               snps.dt.win$variant.id)

###
print("Create ad and dp objects")
ad <- seqGetData(genofile, "annotation/format/AD")
ad <- ad$data
dp <- seqGetData(genofile, "annotation/format/DP")

colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")
####


##### create function
fst.test <- function(sa1, sa2){
  
  samps_to_compare = c(sa1, sa2)
  
  #sa1="SK_0_1"
  #sa2="SKF_1_1"
  
  pool_sizes = c(samps$nFlies[which(samps$Sample_id == sa1)],
                 samps$nFlies[which(samps$Sample_id == sa2)])
  
  ad.matrix = ad[which(rownames(ad) %in% samps_to_compare),]
  rd.matrix = dp[which(rownames(dp) %in% samps_to_compare ),]
  
  pool <- new("pooldata",
              npools=dim(ad.matrix)[1], #### Rows = Number of pools
              nsnp=dim(ad.matrix)[2], ### Columns = Number of SNPs
              refallele.readcount=t(ad.matrix),
              readcoverage=t(rd.matrix),
              poolsizes=pool_sizes * 2,
              poolnames = samps_to_compare)
  
  fst.out <- computeFST(pool, method = "Anova")
  
  data.frame(snp.fst = fst.out$snp.FST) %>% 
    mutate(
      pos = rownames(.),
      chr = reco_wins$chr[k],
      samp1 = sa1,
      samp2 = sa2,
      FST = fst.out$FST,
      winSTA = reco_wins$start[k],
      winEND = reco_wins$end[k]
    ) -> ou
  
  ######
  sim_fst
  
} #### close function









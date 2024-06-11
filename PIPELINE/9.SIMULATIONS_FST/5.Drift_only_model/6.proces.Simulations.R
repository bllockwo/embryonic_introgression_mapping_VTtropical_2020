##### Processing simulation outputs

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
filtering.dt <- get(load("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/SNP.filtering.guide.Rdata"))

##### objects needed for the analyses
samps = fread("/netfiles02/lockwood_lab/IntrogressionProject/Population_files/samp.guide.files.introg.txt")
genofile.path <- "/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds"

### load the genofile into memory
genofile <- seqOpen(genofile.path)

filtering.dt %>%
  filter(is.na(libs)) %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) ->
  chosen.snps
####

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile),
                      allele=seqGetData(genofile, "allele")) %>%
  separate(allele, into = c("ref_allele","alt_allele"), sep = ",")
snps.dt %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) %>%
  right_join(chosen.snps, by = c("chr", "pos", "SNP_id") ) ->
  snps.dt
#### extract 
samps %<>%
  filter(Pop %in% c(
    "SK",
    "VT8"))
setDT(snps.dt)
snps.dt <- snps.dt[nAlleles==2]

snps.dt %>%
  filter(chr == reco_wins$chr[k] & pos > reco_wins$start[k] & pos < reco_wins$end[k] ) ->
  snps.dt.win
#### ==
seqResetFilter(genofile)
seqSetFilter(genofile,
             sample.id=samps$Sample_id,
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

ad %>% t %>% as.data.frame() %>% mutate(SNP_id = rownames(.)) -> ad.real
names(ad.real)[1:2] = c("SK_0_1_ad",  "VT8_0_1_ad")
dp %>% t %>% as.data.frame() %>% mutate(SNP_id = rownames(.)) -> dp.real
names(dp.real)[1:2] = c("SK_0_1_dp",  "VT8_0_1_dp")
left_join(ad.real, dp.real) -> real.dat.slice

#################
### Load sim data
root <- "/gpfs2/scratch/jcnunez/fst_brent/simulations"
set <- paste(reco_wins$chr[k], reco_wins$start[k], reco_wins$end[k], "drift", sep = "_")
filesL <- system(paste("ls ", root, "/", "output", "/",set, sep = "" ), intern = T)

#### --> load base data
info = strsplit(set, "_")
condenced_Set = paste(info[[1]][1],info[[1]][2],info[[1]][3], sep = "_")
meta = fread(paste(root, "/", "win_meta", "/", condenced_Set, ".meta.txt" , sep = "" ) )
names(meta) = c("chr", "level", "start",  "end", "Lenght", "rate_rho")

AFO = fread(paste(root, "/", "win_data", "/", condenced_Set, ".data.txt" , sep = "" ) )
names(AFO) = c("chr", "real_pos", "AF",  "slim_pos")
####

datLoad = 
foreach(fil = filesL, 
        .combine = "rbind", 
        .errorhandling = "remove")%do%{
  
  tmp <- fread(paste(root, "/", "output", "/", set, "/", fil , sep = "" ))
  names(tmp) = c("slim_pos", "AFsim", "chr", "ith" , "recoC", "start", "end", "recoR")
  tmp <- arrange(tmp, slim_pos)
  left_join(AFO, tmp) -> joint.tmp
  
  joint.tmp %<>% mutate(ith = unique(joint.tmp$ith[!is.na(joint.tmp$ith)])) 
  joint.tmp %<>% mutate(recoC = unique(joint.tmp$recoC[!is.na(joint.tmp$recoC)])) 
  joint.tmp %<>% mutate(start = unique(joint.tmp$start[!is.na(joint.tmp$start)])) 
  joint.tmp %<>% mutate(end = unique(joint.tmp$end[!is.na(joint.tmp$end)])) 
  joint.tmp %<>% mutate(recoR = unique(joint.tmp$recoR[!is.na(joint.tmp$recoR)])) 
  joint.tmp$AFsim[is.na(joint.tmp$AFsim)] = 0
  
  joint.tmp %>%
    mutate(SNP_id = paste(chr, real_pos, sep = "_")) %>%
    group_by(chr, real_pos) %>%
    mutate(AFsim_evo = rbinom(1, 60, AFsim)/60) %>%
    mutate(Counts_sim_pool = rbinom(1, 85, AFsim_evo)) %>%
    mutate(Cov_sim_pool = 85) ->
    tmp.processed
  
  left_join(tmp.processed, real.dat.slice) -> tmp.final
  
  return(tmp.final)
  
}

##### save simualtion object
outfol <- "/gpfs2/scratch/jcnunez/fst_brent/simulations/simulated_pools_counts"  
save(datLoad, file = paste(outfol, paste(set,".simPoolCounts.all.Rdata", sep = ""), sep ="/"))


#### calculate FST
fst.df.SK = foreach(i = 1:100, .combine = "rbind" )%dopar%{
  
  message(i)
  datLoad %>%
    filter(ith == i) ->
    dat.tmp.ith
  
  pool_sizes = c(30, 30)
  
  ad.matrix = rbind(dat.tmp.ith$Counts_sim_pool, dat.tmp.ith$SK_0_1_ad)
  rd.matrix = rbind(dat.tmp.ith$Cov_sim_pool, dat.tmp.ith$SK_0_1_dp)

  pool.sk <- new("pooldata",
              npools=dim(ad.matrix)[1], #### Rows = Number of pools
              nsnp=dim(ad.matrix)[2], ### Columns = Number of SNPs
              refallele.readcount=t(ad.matrix),
              readcoverage=t(rd.matrix),
              poolsizes=pool_sizes * 2,
              poolnames = c("sim","SK")
  )
  
  fst.out <- computeFST(pool.sk, method = "Anova")

  data.frame(snp.fst = fst.out$snp.FST) %>% 
    mutate(
      pos_slim = rownames(.),
      pos = dat.tmp.ith$real_pos,
      SNP_id = dat.tmp.ith$SNP_id,
      chr = reco_wins$chr[k],
      samp1 = "sim",
      samp2 = "SK",
      FST = fst.out$FST,
      winSTA = reco_wins$start[k],
      winEND = reco_wins$end[k],
      ith = i
    ) 
  
}## i   

fst.df.VT8 = foreach(i = 1:100, .combine = "rbind" )%dopar%{
  
  message(i)
  datLoad %>%
    filter(ith == i) ->
    dat.tmp.ith
  
  pool_sizes = c(30, 30)
  
  ad.matrix = rbind(dat.tmp.ith$Counts_sim_pool, dat.tmp.ith$VT8_0_1_ad)
  rd.matrix = rbind(dat.tmp.ith$Cov_sim_pool, dat.tmp.ith$VT8_0_1_dp)
  
  pool.vt <- new("pooldata",
                 npools=dim(ad.matrix)[1], #### Rows = Number of pools
                 nsnp=dim(ad.matrix)[2], ### Columns = Number of SNPs
                 refallele.readcount=t(ad.matrix),
                 readcoverage=t(rd.matrix),
                 poolsizes=pool_sizes * 2,
                 poolnames = c("sim","VT8")
  )
  
  fst.out <- computeFST(pool.vt, method = "Anova")
  
  data.frame(snp.fst = fst.out$snp.FST) %>% 
    mutate(
      pos_slim = rownames(.),
      pos = dat.tmp.ith$real_pos,
      SNP_id = dat.tmp.ith$SNP_id,
      chr = reco_wins$chr[k],
      samp1 = "sim",
      samp2 = "VT8",
      FST = fst.out$FST,
      winSTA = reco_wins$start[k],
      winEND = reco_wins$end[k],
      ith = i
    ) 
  
}## i   

####
fst.all.df.joint = rbind(fst.df.VT8,fst.df.SK)

####
outfol <- "/gpfs2/scratch/jcnunez/fst_brent/simulations/fst_sims_all"  
save(fst.all.df.joint, file = paste(outfol, paste(set,".fst.all.Rdata", sep = ""), sep ="/"))

####

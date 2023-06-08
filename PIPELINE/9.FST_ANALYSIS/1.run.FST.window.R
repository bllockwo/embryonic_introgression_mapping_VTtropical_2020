##### Calculate FST in Introgression Lines

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
##array 1-231
print(k)

##### objects needed for the analyses

# samp files
samps = fread("/netfiles02/lockwood_lab/IntrogressionProject/Population_files/samp.guide.files.introg.txt")
removed.samps = c("CH_0_2","VT10_0_2","VT10F2_2_1","CHF3_3_1")

samps %<>%
  filter(!Sample_id %in% removed.samps)

#####
#Generate outfile object
L = dim(samps)[1]

message("Create combination vector")
comp_vector = combinations(
  L,
  2, 
  v=1:L,
  set=TRUE, 
  repeats.allowed=FALSE)

###
message("Loop to esrtimate time difference")
data.frame(comp_vector,
           samp1 = rep(NA, dim(comp_vector)[1]),
           samp2 = rep(NA, dim(comp_vector)[1])
) -> comp_vector.t

comp.pairs =
foreach(i = 1:dim(comp_vector)[1], 
        .combine = "rbind", 
        .errorhandling = "remove")%do%{
          
  message(i)
  id1 = comp_vector.t[i,1]
  id2 = comp_vector.t[i,2]
  
  samp1 = samps$Sample_id[id1]
  samp2 = samps$Sample_id[id2]
  
  order.vect <- sort(c(samp1, samp2))
  
  tmp <- data.frame(samp1, samp2, ordered.vec = paste(order.vect[1],order.vect[2] , sep = "|") )
  
  return(tmp)
}

#####
genofile.path <- "/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds"

### load the genofile into memory
genofile <- seqOpen(genofile.path)

######
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

#####
snps.dt %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) %>%
  right_join(chosen.snps, by = c("chr", "pos", "SNP_id") ) ->
  snps.dt
####
#snp.dt <- snp.dt[cm_mb>0 & !is.na(cm_mb)]
snps.dt <- snps.dt[nAlleles==2]

#### ==
seqSetFilter(genofile,
             snps.dt[chr%in%c("2L", "2R", "3L", "3R", "X")][missing<.05]$variant.id)

####
win.bp = 1e6
step.bp = 0.5e6


wins <- foreach(chr.i=c("2L","2R","3L","3R","X"),
                .combine="rbind", 
                .errorhandling="remove")%dopar%{
                  
                  tmp <- snps.dt %>%
                    filter(chr == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]
####
dim(wins)
message(paste(wins$chr[k], wins$start[k], wins$end[k], sep = "-"))

###
print("Create ad and dp objects")

ad <- seqGetData(genofile, "annotation/format/AD")
ad <- ad$data
dp <- seqGetData(genofile, "annotation/format/DP")

print("Create dat object")

#Add metadata ad
colnames(ad) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(ad) <- seqGetData(genofile, "sample.id")

#Add metadata dp
colnames(dp) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dp) <- seqGetData(genofile, "sample.id")

###
snps.dt %>%
  filter(chr == wins$chr[k] & pos > wins$start[k] & pos < wins$end[k] ) %>%
  .$SNP_id ->
  snps_of_win

###
fst.df = foreach(i = 1:dim(comp.pairs)[1], .combine = "rbind" )%dopar%{
  
  print(i/dim(comp.pairs)[1] * 100)
  
  samps_to_compare = c(comp.pairs$samp1[i], comp.pairs$samp2[i])
  
  pool_sizes = c(samps$nFlies[which(samps$Sample_id == comp.pairs$samp1[i])],
                 samps$nFlies[which(samps$Sample_id == comp.pairs$samp2[i])])
  
  ad.matrix = ad[which(rownames(ad) %in% samps_to_compare), 
                 which(colnames(ad) %in% snps_of_win)]

  rd.matrix = dp[which(rownames(dp) %in% samps_to_compare ),
                 which(colnames(dp) %in% snps_of_win)]
  
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
    chr = wins$chr[k],
    samp1 = comp.pairs$samp1[i],
    samp2 = comp.pairs$samp2[i],
    FST = fst.out$FST,
    winSTA = wins$start[k],
    winEND = wins$end[k],
  ) 
  
}## i   


###save
Filename = paste(wins$chr[k], wins$start[k], wins$end[k], sep = "-")
outfile <- "~/scratch/fst_brent/fst_out"

save(fst.df, 
     file = paste(outfile,"/Fst.",Filename, ".Rdata" , sep = "" ) )


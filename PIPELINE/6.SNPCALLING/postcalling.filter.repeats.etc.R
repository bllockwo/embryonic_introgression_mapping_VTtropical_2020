#### Filter SNPs

## module load intel/18.0 intelmpi/18.0 R/3.6.3; R

### libraries
library(data.table)
library(gdata)
library(lubridate)
library(foreach)
library(SeqArray)
library(tidyr)
library(doMC)
registerDoMC(5)
#  library(patchwork)
#  library(ggplot2)

### open GDS
genofile <- seqOpen("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds", allow.duplicate=T)

### load BED files
beds.fn <- list.files("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/RepeatFilterFiles/", ".bed", full.names=T)
beds.fn <- beds.fn[!grepl("combine", beds.fn)]

rep.bed <- foreach(i=beds.fn, .combine="rbind")%do%{
  # i<- beds.fn[5]
  message(i)
  tmp <- fread(i, skip=1)
  tmp[,chr:=gsub("chr", "", V1)]
  setnames(tmp, c("V2", "V3"), c("start", "end"))
  ### tmp[chr=="2R"][start<=6481518][end>=6481518]
  tmp <- tmp[,c("chr", "start", "end"), with=F]
  tmp[,repLib:=tstrsplit(i, "/") %>% last %>% gsub(".bed", "", .)]
  tmp # tmp[chr=="2R"][start<=6481518][end>=6481518]
}
rep.bed[chr=="2R"][start<=6481518][end>=6481518]

### load basic SNP table
snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     nAlleles=seqGetData(genofile, "$num_allele"),
                     id=seqGetData(genofile, "variant.id"),
                     DP=seqGetData(genofile, "annotation/info/DP"))


snp.dt <- snp.dt[nAlleles==2]
setkey(snp.dt)
snp.dt[,start:=pos]
snp.dt[,end:=pos]

### merge repetitive regions from BED files with SNP table
setkey(snp.dt, start, end)
setkey(rep.bed, start, end)

snp.dt <- foreach(chr.i=c("2L", "2R", "3L", "3R", "X"))%dopar%{
  #chr.i<-"2L"
  print(chr.i)
  tmp <- foverlaps(snp.dt[chr==chr.i], rep.bed[chr==chr.i])
  tmp.ag <- tmp[,list(.N, libs=paste(repLib, collapse=";"), DP=mean(DP)), list(chr=i.chr, pos, id)]
  tmp.ag[libs=="NA", libs:=NA]
  tmp.ag[is.na(libs), N:=0]
  tmp.ag
}

snp.dt <- rbindlist(snp.dt)
snp.dt[,start:=pos]
snp.dt[,end:=pos]
setkey(snp.dt, start, end)

### load in recombination rate file
recRate <- fread("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/RecRates-All-Chromosomes-100kb.r6.bed")
setnames(recRate, names(recRate), c("chr", "start", "end", "cm_mb"))
recRate[,chr:=gsub("chr", "", chr)]
recRate
setkey(recRate, start, end)

snp.dt <- foreach(chr.i=c("2L", "2R", "3L", "3R", "X"))%dopar%{
  #chr.i<-"2L"
  print(chr.i)
  tmp <- foverlaps(snp.dt[chr==chr.i], recRate[chr==chr.i])
  tmp[,c("chr", "pos", "id", "N", "libs", "cm_mb"), with=F]
}
snp.dt <- rbindlist(snp.dt)
snp.dt[,start:=pos]
snp.dt[,end:=pos]
setkey(snp.dt, start, end)

### add in Inversion identifier
inv.dt <- fread("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/InversionsMap_hglft_v6_inv_startStop.txt")
setnames(inv.dt, "chrom", "chr")
setnames(inv.dt, "stop", "end")
setkey(inv.dt, start, end)

snp.dt <- foreach(chr.i=c("2L", "2R", "3L", "3R", "X"))%dopar%{
  #chr.i<-"3R"
  print(chr.i)
  tmp <- foverlaps(snp.dt[chr==chr.i], inv.dt[chr==chr.i])
  tmp <- tmp[,c("i.chr", "pos", "id", "N", "libs", "cm_mb", "invName"), with=F]
  dim(tmp)
  tmp2 <- tmp[,list(invName=paste(invName, collapse=";")), list(chr=i.chr, pos, id, N, libs, cm_mb)]
  tmp2[invName=="NA", invName:="none"]
  
  tmp2
}
snp.dt <- rbindlist(snp.dt)

save(snp.dt, file="/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/SNP.filtering.guide.Rdata")


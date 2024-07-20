##### PCA analyses with simulation and real data
library(tidyverse)
library(magrittr)
library(reshape2)
library(SeqArray)
library(FactoMineR)
library(factoextra)
library(data.table)
library(foreach)

out="AFS_drift_objects"
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
             snps.dt[chr%in%c("2L", "2R", "3L", "3R")][missing<.05]$variant.id)

### get allele frequency data
ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")

dat <- ad$data/dp
dim(dat)  

colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")

###### Bring in simulations 

tmp_fil = system("ls /gpfs2/scratch/jcnunez/fst_brent/simulations/DrifOnly_collect_COUNTS", intern = T)
root = "/gpfs2/scratch/jcnunez/fst_brent/simulations/DrifOnly_collect_COUNTS"
  
foreach(i=1:length(tmp_fil),
                    .combine = "rbind",
                    .errorhandling = "remove")%do%{
  
    message(i)
    info_tmp = tmp_fil[i]
    ith =  strsplit(info_tmp, "\\.")[[1]][5]
                                      
    tmp_dat = get(load(paste(root, info_tmp, sep = "/")))
    
    tmp_dat %>%
      group_by(chr, real_pos) %>%
      arrange(chr, real_pos) -> tmp_dat.sort
    
    sim_af = data.frame(t(tmp_dat$Counts_sim_pool/tmp_dat$Cov_sim_pool))
    names(sim_af) = tmp_dat.sort$SNP_id
    rownames(sim_af) = paste("drift", ith, sep = "_")
    #return(sim_af)
    save(sim_af, 
         file = paste(out,
                      paste("drift", ith, "Rdata" , sep = "."), 
                      sep ="/"))
}

###



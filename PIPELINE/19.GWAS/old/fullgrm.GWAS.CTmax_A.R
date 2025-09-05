### libraries
library(GMMAT)
library(data.table)
library(magrittr)
library(tidyverse)
library(foreach)
library(SeqArray)
library(SeqVarTools)

POPSET = "AA_samps_DGRP.txt"

args = commandArgs(trailingOnly=TRUE)
#create job id from array
ith = as.numeric(args[1])


# saveRDS(task.id, paste0(INPUT, "newphenotaskidchart"))
INPUT = "/gpfs2/scratch/jcnunez/fst_brent/GWAS_EB_OK/"#path of working directory

#load grm
grm = readRDS(paste0(INPUT, "GRM.Aug"))
#load in our phenotype file 
phenotable = readRDS(paste0(INPUT, "wideform.fixed.phenotable.RDS"))

wolbachia = fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/wolbachia.status.txt")
names(wolbachia)[1] = "ral_id"
#fix ral ids
dt = phenotable %>% 
  mutate(ral_id = gsub("^", "line_", ral_id)) %>%
  dplyr::select(ral_id, HighThermalToleranceExtreme_VaryingWithTemperature_F) %>%
  left_join(wolbachia) %>%
  filter(!is.na(HighThermalToleranceExtreme_VaryingWithTemperature_F))

samps_ids <- fread(POPSET, header = F)

dt %<>%
  filter(ral_id %in% samps_ids$V1)

###filter stuff .. grm
grm_Set = grm[dt$ral_id,dt$ral_id]

####
ingds = paste0(INPUT,"dgrp.gds")
genofile <- seqOpen(ingds)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))

###
win.bp = 1e6
step.bp = 1e6+1


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

win.i = wins[ith]

snps.dt.flt = snps.dt %>% 
  filter(chr == win.i$chr) %>%
  filter(pos >= win.i$start) %>%
  filter(pos <= win.i$end) 

#GMMAT creates a null model that uses wolbachia status as fixed effect on phenotype, and the identity matrix as a (nonexistent) random effect. family refers to the method used to calculate the model
#fit a GLMM to the data
modelqtl <- glmmkin(fixed =  HighThermalToleranceExtreme_VaryingWithTemperature_F ~ Infection_Status, 
                    data = dt, 
                    kins = grm_Set, id = "ral_id",
                    family = gaussian(link = "identity"))

#run GMMAT
#the glmm.score function scores the impact of each snp agains the null model. minor allele frequency and missing data cutoff can be specified to reduce unwanted computation. compututation can also be sped up by specifying the number of snps calculated per batch, and how many cores to use
seqSetFilter(genofile, 
             sample.id = dt$ral_id,
             variant.id = snps.dt.flt$variant.id
)

####
system("mkdir GWAS_out_AA")
glmm.score(modelqtl, infile = genofile, 
           outfile = paste("./GWAS_out_AA/batch", win.i$i, "txt" , sep = "."), 
           MAF.range = c(0.05,0.95), 
           miss.cutoff = 0.15, 
           nperbatch = 400, 
           ncores = 5)

print ("done")

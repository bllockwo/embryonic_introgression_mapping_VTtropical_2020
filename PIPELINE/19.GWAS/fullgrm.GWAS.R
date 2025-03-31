### libraries
library(GMMAT)
library(data.table)
library(tidyverse)
library(foreach)
library(SeqArray)
library(SeqVarTools)

args = commandArgs(trailingOnly=TRUE)
#create job id from array
ith = as.numeric(args[1])


# saveRDS(task.id, paste0(INPUT, "newphenotaskidchart"))
INPUT = "/gpfs2/scratch/jcnunez/fst_brent/GWAS_EB_OK/"#path of working directory

#load grm
grm = readRDS(paste0(INPUT, "GRM.Aug"))
#load in our phenotype file 
phenotable = fread(paste0(INPUT, "Mar28.data.txt"))
#fix ral ids
dt = phenotable %>% 
  dplyr::rename("ral_id" = Line) %>% 
  mutate(ral_id = gsub("DGRP", "line", ral_id)) %>%
  mutate(WolbNum = case_when(Wolbachia_Status == "Yes" ~ 1,
                             Wolbachia_Status == "No" ~ 0
  ))

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
modelqtl <- glmmkin(fixed =  Proportion ~ Wolbachia_Status, 
                    data = dt, 
                    kins = grm_Set, id = "ral_id",
                    family = binomial(link = "logit"))

#run GMMAT
#the glmm.score function scores the impact of each snp agains the null model. minor allele frequency and missing data cutoff can be specified to reduce unwanted computation. compututation can also be sped up by specifying the number of snps calculated per batch, and how many cores to use
seqSetFilter(genofile, 
             sample.id = dt$ral_id,
             variant.id = snps.dt.flt$variant.id
)
####
system("mkdir GWAS_out")
glmm.score(modelqtl, infile = genofile, 
           outfile = paste("./GWAS_out/batch", win.i$i, "txt" , sep = "."), 
           MAF.range = c(0.05,0.95), 
           miss.cutoff = 0.15, 
           nperbatch = 400, 
           ncores = 5)

print ("done")



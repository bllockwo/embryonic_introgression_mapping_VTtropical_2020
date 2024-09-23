library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)

### gwas loc
loc <- "/gpfs2/scratch/jcnunez/fst_brent/phenotype_analysis/lpgwas/"
fils <- system(paste("ls ", loc), intern = T)

topsnp_chr = "2R"
topsnp_pos = "16439138"
  
top.snps =
  foreach(i=1:length(fils), .combine = "rbind")%do%{
    message(i)
    
### load data
tam <- fread(paste(loc,fils[i], sep = "/")) %>%
          mutate(pheno = fils[i])
p.adj = p.adjust(tam$PVAL, "fdr")

tam %>%
  mutate(p.adj = p.adj) %>%
  filter(CHR == topsnp_chr) %>%
  filter(POS == topsnp_pos) ->
  tmp.snp

return(tmp.snp)
  }

top.snps %>% filter(PVAL < 0.1) ->
  signif_hits

### Load phenos
wolbachia_st <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/wolbachia.status.txt")
C_lines <- fread("/gpfs2/scratch/jcnunez/fst_brent/DGRP_analysis/snp0.txt", header = F) %>%
  mutate(snp_stat = "C")
A_lines <- fread("/gpfs2/scratch/jcnunez/fst_brent/DGRP_analysis/snp2.txt", header = F) %>%
  mutate(snp_stat = "A")
invs <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/Inversion.status.txt")
names(invs)[1] = "line_id"

rbind(C_lines, A_lines) ->
  lines_stats
names(lines_stats)[1] = "line_id"

phenos <- readRDS("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2024.Nunez_et_al_Genetics/Phenotyping/wideform.fixed.phenotable.RDS")
phenos %>%
  mutate(line_id = paste("line",ral_id, sep = "_")) %>%
  left_join(lines_stats) ->
  pheno_embryo_SNP

phenos.targ = which(names(pheno_embryo_SNP) %in% gsub(".txt","",signif_hits$pheno))

pheno_embryo_SNP = as.data.frame(pheno_embryo_SNP)
pheno_embryo_SNP[,c(phenos.targ,
                    which(names(pheno_embryo_SNP) %in% 
                            c("line_id", "snp_stat")))] -> flt.dat
                       
flt.dat %>%
  filter(!is.na(snp_stat)) %>%
  melt(id = c("line_id", "snp_stat")) %>%
  left_join(wolbachia_st) %>%
  left_join(invs)->
  flt.dat.melt

linear_anovas = 
foreach(i=unique(flt.dat.melt$variable),
        .combine = "rbind")%do%{
          
          lm(value~snp_stat+Infection_Status+`In(2R)NS`, 
                 data = filter(flt.dat.melt, variable == i )) %>%
            anova ->
            test.out
          
          data.frame(
            pheno = i,
            p_snp = test.out$`Pr(>F)`[1],
            p_wolb = test.out$`Pr(>F)`[2],
            p_inv = test.out$`Pr(>F)`[3]
            
          )
        }

linear_anovas %>%
  filter(p_snp < 0.05)

flt.dat.melt %>%
  filter(variable == 
           "HighThermalToleranceExtreme_VaryingWithTemperature_M") %>%
  ggplot(aes(
    x= snp_stat,
    y= value,
    fill =Infection_Status
  )) + geom_boxplot() ->
  phenos.plot

ggsave(phenos.plot, file = "phenos.plot.pdf",
       h = 3.5, w = 4)

filter(flt.dat.melt,
       variable == 
         "HighThermalToleranceExtreme_VaryingWithTemperature_M" 
 ) -> trait.frT.data

t.test(value~snp_stat,
       data = trait.frT.data) 

t.test(value~Infection_Status,
       data = filter(trait.frT.data,
                     (snp_stat == "A" & Infection_Status == "n") |
                     (snp_stat == "A" & Infection_Status == "y") 
                     )) 

t.test(value~Infection_Status,
       data = filter(trait.frT.data,
                     (snp_stat == "C" & Infection_Status == "n") |
                       (snp_stat == "C" & Infection_Status == "y") 
       )) 

t.test(value~Infection_Status,
       data = filter(trait.frT.data,
                     (snp_stat == "C" & Infection_Status == "y") |
                       (snp_stat == "A" & Infection_Status == "n") 
       )) 


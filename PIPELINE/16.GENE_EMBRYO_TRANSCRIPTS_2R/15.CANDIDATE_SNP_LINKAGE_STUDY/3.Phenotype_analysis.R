library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)

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

embryos <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/Embryonic_Thermal_Tolerance/Sep9_data.txt")
names(embryos)[1] = "line_id"
embryos$line_id = gsub("DGRP", "line", embryos$line_id)

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

##extract phenos

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


####### Include embryos
####### Include embryos
####### Include embryos
####### Include embryos

pheno_embryo_SNP %>%
  as.data.frame() %>%
  dplyr::select(HighThermalToleranceExtreme_VaryingWithTemperature_M,
                line_id, snp_stat) %>%
  left_join(wolbachia_st) ->
  CTmax.data

CTmax.data %<>%
  mutate(mean.phen = mean(CTmax.data$`HighThermalToleranceExtreme_VaryingWithTemperature_M`, na.rm = T))
embryos %<>%
  mutate(emb.phen = mean(embryos$Proportion, na.rm = T))

  
CTmax.data %>%
  full_join(embryos) ->
  joint.data

joint.data %>%
  ggplot(aes(
    x=HighThermalToleranceExtreme_VaryingWithTemperature_M,
    y=Proportion,
    color = Infection_Status
  )) + geom_point() +
  geom_smooth(method = "lm", color = "black")  ->
  ant.plei

ggsave(ant.plei, file = "ant.plei.pdf")

CTmax.data %>%
  left_join(embryos) %>%
  ggplot(aes(
    x=HighThermalToleranceExtreme_VaryingWithTemperature_M,
    y=Proportion,
    color = snp_stat,
    shape = Infection_Status
  )) + geom_point() +
  geom_smooth(method = "lm", color = "black")  ->
  ant.plei

ggsave(ant.plei, file = "ant.plei.pdf")

### Stock
stocks <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/Embryonic_Thermal_Tolerance/stock_Assessment.sep23.txt")
stocks %>%
  dplyr::select(Stock=Stock,
                line_id=Genotype) ->
  stocks.slic

joint.data$line_id = gsub("DGRP","line", stocks.slic$line_id)

joint.data %>%
  left_join(stocks.slic) %>%
  mutate(completion = case_when(
    !is.na(HighThermalToleranceExtreme_VaryingWithTemperature_M) &
    !is.na(Proportion) ~ "complete",
    
    is.na(HighThermalToleranceExtreme_VaryingWithTemperature_M) &
    !is.na(Proportion) ~ "miss.CTmax",
    
    !is.na(HighThermalToleranceExtreme_VaryingWithTemperature_M) &
    is.na(Proportion) ~ "miss.embryo",  
    
    is.na(HighThermalToleranceExtreme_VaryingWithTemperature_M) &
    is.na(Proportion) ~ "miss.both",  
  )) ->
  completion.table

write.table(completion.table, 
            file = "completion.table.Sep23.2024.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

completion.table %>%
  filter(completion == "miss.embryo")
  

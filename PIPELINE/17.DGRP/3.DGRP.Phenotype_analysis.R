library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)

### Load phenos
wolbachia_st <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/wolbachia.status.txt")
invs <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/Inversion.status.txt")
names(invs)[1] = "line_id"

SNPX <- fread("/gpfs2/scratch/jcnunez/fst_brent/X_hits/X_snp_info.txt", header = T)
names(SNPX)[1] = "line_id"

SNP2R <- fread("/gpfs2/scratch/jcnunez/fst_brent/X_hits/2R_snp_info.txt", header = T)
names(SNP2R)[1] = "line_id"

embryos <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/Embryonic_Thermal_Tolerance/Final_data.txt")
names(embryos)[1] = "line_id"
embryos$line_id = gsub("DGRP", "line", embryos$line_id)


phenos <- readRDS("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2024.Nunez_et_al_Genetics/Phenotyping/wideform.fixed.phenotable.RDS")
phenos %>%
  mutate(line_id = paste("line",ral_id, sep = "_")) %>%
  full_join(SNPX) %>%
  full_join(SNP2R) %>%
  full_join(embryos) ->
  pheno_embryo_SNP



### gwas loc
loc <- "/gpfs2/scratch/jcnunez/fst_brent/phenotype_analysis/lpgwas/"
fils <- system(paste("ls ", loc), intern = T)

topsnp_chr1 = "2R"
topsnp_pos1 = "16439138"
topsnp_chr2 = "X"
topsnp_pos2 = "15496974"

top.snps =
  foreach(i=1:length(fils), .combine = "rbind")%do%{
    message(i)
    
    ### load data
    tam <- fread(paste(loc,fils[i], sep = "/")) %>%
      mutate(pheno = fils[i])
    #pheno="HighThermalToleranceExtreme_VaryingWithTemperature_M.txt"; tam <- fread(paste(loc,pheno, sep = "/")) %>% mutate(pheno = pheno)
    
    p.adj = p.adjust(tam$PVAL, "fdr")
    
    tam %>%
      mutate(p.adj = p.adj) %>%
      filter(CHR == topsnp_chr1) %>%
      filter(POS == topsnp_pos1) ->
      tmp.snp1
    
    tam %>%
      mutate(p.adj = p.adj) %>%
      filter(CHR == topsnp_chr2) %>%
      filter(POS == topsnp_pos2) ->
      tmp.snp2
    
    o <- rbind(tmp.snp1, tmp.snp2)
    
    return(o)
  }

save(top.snps, file = "top.hits.X.2R.GWAS.rdata")

top.snps %>% filter(PVAL < 0.1) ->
  signif_hits

signif_hits %>%
  filter(CHR == "2R")
signif_hits %>%
  filter(CHR == "X")



##extract phenos

#phenos.targ = which(names(pheno_embryo_SNP) %in% gsub(".txt","",signif_hits$pheno))

phenos.targ = c("HighThermalToleranceExtreme_VaryingWithTemperature_F",                            
                "HighThermalToleranceExtreme_VaryingWithTemperature_M", "Proportion") 

pheno_embryo_SNP = as.data.frame(pheno_embryo_SNP)
pheno_embryo_SNP$Proportion %>% is.na %>% table

pheno_embryo_SNP[,c(
  which(names(pheno_embryo_SNP) %in% 
          c("line_id", "snp_stat",phenos.targ, "Xsnp", "snp2R")))] %>%
  left_join(wolbachia_st) %>%
  filter(Infection_Status == "n")-> flt.dat

########   

### CTmax

flt.dat %>%
  filter(Xsnp %in% c("A/A","G/G")) %>%
  filter(!is.na(HighThermalToleranceExtreme_VaryingWithTemperature_F)) %>%
  t.test(HighThermalToleranceExtreme_VaryingWithTemperature_F~as.factor(Xsnp), data = .)
flt.dat %>%
  filter(Xsnp %in% c("A/A","G/G")) %>%
  filter(!is.na(HighThermalToleranceExtreme_VaryingWithTemperature_M)) %>%
  t.test(HighThermalToleranceExtreme_VaryingWithTemperature_M~as.factor(Xsnp), data = .)

flt.dat %>%
  filter(snp2R %in% c("A/A","C/C")) %>%
  filter(!is.na(HighThermalToleranceExtreme_VaryingWithTemperature_M)) %>%
  t.test(HighThermalToleranceExtreme_VaryingWithTemperature_M~as.factor(snp2R), data = .)

####
  flt.dat %>%
    filter(snp2R %in% c("A/A","C/C") & Xsnp %in% c("A/A","G/G")) %>%
    select(snp2R, Xsnp, HighThermalToleranceExtreme_VaryingWithTemperature_M, HighThermalToleranceExtreme_VaryingWithTemperature_F ) %>%
    melt(id = c("snp2R", "Xsnp"), variable.name = "pheno") %>%
    melt(id = c("pheno", "value"), variable.name = "genot", value.name = "alleles") %>%
    mutate(allele_type = case_when(genot == "snp2R" & alleles == "C/C" ~ "tropical",
                                   genot == "snp2R" & alleles == "A/A" ~ "temperate",
                                   genot == "Xsnp" & alleles == "A/A" ~ "tropical",
                                   genot == "Xsnp" & alleles == "G/G" ~ "temperate",
                                   )) %>%
    ggplot(
      aes(
        x=allele_type,
        y=value,
        fill=pheno
      )
    ) + geom_boxplot() + facet_grid(pheno~genot, scales = "free") +
    theme(legend.position = "bottom")-> 
    boxplot_prop
  ggsave(boxplot_prop, file = "boxplot_prop.pdf",
         h = 6, w = 6)

####
  flt.dat %>%
    filter(snp2R %in% c("A/A","C/C") & Xsnp %in% c("A/A","G/G")) %>%
    select(snp2R, Xsnp, 
           HighThermalToleranceExtreme_VaryingWithTemperature_M, 
           HighThermalToleranceExtreme_VaryingWithTemperature_F,
           Proportion) %>%
    melt(id = c("snp2R", "Xsnp", "Proportion"), variable.name = "pheno") %>%
    melt(id = c("pheno", "value", "Proportion"), variable.name = "genot", value.name = "alleles") %>%
    mutate(allele_type = case_when(genot == "snp2R" & alleles == "C/C" ~ "tropical",
                                   genot == "snp2R" & alleles == "A/A" ~ "temperate",
                                   genot == "Xsnp" & alleles == "A/A" ~ "tropical",
                                   genot == "Xsnp" & alleles == "G/G" ~ "temperate",
    )) %>%
    ggplot(
      aes(
        x=value, y =Proportion,
        color=pheno
      )
    ) + geom_point() + geom_smooth(method = "lm") +
    facet_grid(genot~allele_type) +
    theme(legend.position = "bottom")-> 
    hatch_adult
  ggsave(hatch_adult, file = "hatch_adult.pdf",
         h = 6, w = 6)
  
####
  flt.dat %>%
    filter(snp2R %in% c("A/A","C/C") & Xsnp %in% c("A/A","G/G")) %>%
    select(snp2R, Xsnp, 
           HighThermalToleranceExtreme_VaryingWithTemperature_M, 
           HighThermalToleranceExtreme_VaryingWithTemperature_F,
           Proportion) %>%
    melt(id = c("snp2R", "Xsnp", "Proportion"), variable.name = "pheno") %>%
    mutate(joint_G = paste(snp2R, Xsnp)) %>%
    ggplot(
      aes(
        x=value, y =Proportion,
        color=pheno
      )
    ) + geom_point() + geom_smooth(method = "lm", se = F) +
    facet_grid(.~joint_G) +
    theme(legend.position = "bottom")-> 
    hatch_adult_joint
  ggsave(hatch_adult_joint, file = "hatch_adult_joint.pdf",
         h = 6, w = 6)
  

  
  
  
  
########                       
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
  filter(Infection_Status == "n") %>%
  ggplot(aes(
    x= snp_stat,
    y= value,
    fill =Infection_Status
  )) + geom_boxplot() + theme_bw() + 
  facet_grid(~variable) ->
  phenos.plot

ggsave(phenos.plot, file = "phenos.plot.pdf",
       h = 3.5, w = 5)

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
  dplyr::select(HighThermalToleranceExtreme_VaryingWithTemperature_F,
                line_id, snp_stat) %>%
  left_join(wolbachia_st) ->
  CTmax.data

CTmax.data %<>%
  mutate(mean.phen = mean(CTmax.data$`HighThermalToleranceExtreme_VaryingWithTemperature_F`, na.rm = T))
embryos %<>%
  mutate(emb.phen = mean(embryos$Proportion, na.rm = T))

  
CTmax.data %>%
  full_join(embryos) ->
  joint.data

joint.data %>%
  filter(Infection_Status == "n") %>%
  ggplot(aes(
    x=HighThermalToleranceExtreme_VaryingWithTemperature_F,
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
  

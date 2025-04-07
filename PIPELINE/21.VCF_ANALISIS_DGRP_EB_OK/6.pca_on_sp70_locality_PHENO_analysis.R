## PLOT PCA
#Load libraries
library(tidyverse)
library(vcfR)
library(data.table)
library(adegenet)
library(FactoMineR)
library(magrittr)
library(reshape2)
library(zoo)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(ggrepel)
library(lubridate)
library(foreach)
library(SeqArray)
library(doMC)
registerDoMC(8)
library(fastglm)
library(lmtest )
library(gmodels)
library(lme4)

# Load the VA vcf
VAvcf <- read.vcfR(
  "/gpfs2/scratch/jcnunez/fst_brent/slice_VCFs/top2R_1MB_VA.recode.vcf.gz")
VAgl <- vcfR2genlight(VAvcf) 

###
tab(VAgl, NA.method = "asis") %>%
  as.data.frame() %>%
  t() ->
  VA_data

VA_data %<>%
  as.data.frame() %>%
  mutate(SNP_id = rownames(.)) %>%
  separate(SNP_id, remove = FALSE,
           into = c("chr","pos","feature"))

VA_data %>%
  filter(pos > 20551633-5000) %>%
  filter(pos < 20551633+5000) ->
  VA_focal
  

colnames(VA_data) -> VA_samps
data.frame(sampleid=VA_samps) %>%
  separate(sampleid, 
           remove = F, into = c("loc", 
                                "samp", 
                                "time"),
           sep = "_") ->
  VA_meta

VA_meta %<>%
  mutate( dat_spec = case_when(loc == "CM" ~ time,
                               loc %in% c("OW1","OW2") ~ NA
  )) %>% 
  separate(dat_spec, into = c("MO", "DAY"), sep = 2) %>%
  mutate(Date_corr = as.Date(paste(DAY, MO, 2016), 
                             format = "%d %m %Y"))


#####

WORLDvcf <- read.vcfR(
  "/gpfs2/scratch/jcnunez/fst_brent/slice_VCFs/top2R_1MB_World.recode.vcf.gz")
WORLDgl <- vcfR2genlight(WORLDvcf) 

###
tab(WORLDgl, NA.method = "asis") %>%
  as.data.frame() %>%
  t() ->
  WORLD_data

WORLD_data %<>%
  as.data.frame() %>%
  mutate(SNP_id = rownames(.)) %>%
  separate(SNP_id, remove = FALSE,
           into = c("chr","pos","feature"))

WORLD_data %>%
  filter(pos > 20551633-5000) %>%
  filter(pos < 20551633+5000) ->
  WORLD_focal

ME<-fread("/gpfs2/scratch/jcnunez/fst_brent/slice_VCFs/Sample_guides/All_samps_ME.txt", header = F)
PA<-fread("/gpfs2/scratch/jcnunez/fst_brent/slice_VCFs/Sample_guides/All_samps_PA.txt", header = F)
NET<-fread("/gpfs2/scratch/jcnunez/fst_brent/slice_VCFs/Sample_guides/All_samps_Netherlands.txt", header = F)
ZIM<-fread("/gpfs2/scratch/jcnunez/fst_brent/slice_VCFs/Sample_guides/All_samps_Zimb.txt", header = F)

samps_Wo <- c(ME$V1, PA$V1, NET$V1, ZIM$V1)

WORLD_focal[,which(colnames(WORLD_focal) %in% samps_Wo)] %>%
  as.data.frame() %>%
  mutate(SNP_id = rownames(.)) %>%
  separate(SNP_id, remove = FALSE,
           into = c("chr","pos","feature"))->
  WORLD_focal

#### DGRP
DGRPvcf <- read.vcfR(
  "/gpfs2/scratch/jcnunez/fst_brent/slice_VCFs/top2R_1MB_DGRP.recode.vcf.gz")
DGRPgl <- vcfR2genlight(DGRPvcf) 

###
tab(DGRPgl, NA.method = "asis") %>%
  as.data.frame() %>%
  t() ->
  DGRP_data

DGRP_data %<>%
  as.data.frame() %>%
  mutate(SNP_id = rownames(.)) %>%
  separate(SNP_id, remove = FALSE,
           into = c("chr","pos","feature"))

DGRP_data %<>% mutate(pos = as.numeric(pos) + 4112495)
DGRP_data %>% filter(pos == 20551633)

####
VA_focal$pos = as.numeric(VA_focal$pos)
WORLD_focal$pos = as.numeric(WORLD_focal$pos)

inner_join(VA_focal, WORLD_focal, by = c("chr","pos","feature")) %>%
  inner_join(DGRP_data,  by = c("chr","pos","feature")) ->
  joint_data

joint_data %>%
  filter(SNP_id.y == "2R_20551633_SNP") %>%
  t() %>% data.frame(SP70 = ., sampleid = rownames(.)) %>%
  filter(SP70 %in% c(0,2)) ->
  homozyg_SP70

joint_data[,which(colnames(joint_data) %in% 
                    homozyg_SP70$sampleid)] ->
  dat_for_pca


dat_for_pca %>% t() ->
  dat_for_pca_t

save(dat_for_pca_t, file = "data.for.SP70_30k.pca_joint.Rdata")

dat_for_pca_t %>%
  PCA(graph = F) ->
  PCA_obj

###save


####
PCA_obj$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleid = rownames(.)) %>%
  left_join(homozyg_SP70) %>%
  mutate(pop = case_when( grepl("LN", sampleid) ~ "PA",
                          grepl("CM", sampleid) ~ "VA",
                          grepl("ME", sampleid) ~ "ME",
                          grepl("Nether", sampleid) ~ "NL",
                          grepl("ZI", sampleid) ~ "ZI",
                          grepl("line", sampleid) ~ "DGRP"
  )) ->
  pca_coords_annot

pca_coords_annot %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    shape = pop,
    color = SP70
  )) + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2.1)+
  theme_bw() +
  scale_color_manual(values = c("firebrick","steelblue"))->
  pca.plot
ggsave(pca.plot, file = "pca.plot.pdf",
       w=6, h =3)

#### Phenotyping 
phenos <- readRDS("wideform.fixed.phenotable.RDS")
names(phenos)

MCT = "HighThermalToleranceExtreme_VaryingWithTemperature_F"
FCT = "HighThermalToleranceExtreme_VaryingWithTemperature_M"

data.frame(
  sampleid = paste("line",phenos$ral_id, sep = "_"),
  CT_F = phenos[,"HighThermalToleranceExtreme_VaryingWithTemperature_F"],
  CT_M = phenos[,"HighThermalToleranceExtreme_VaryingWithTemperature_M"]
  ) -> phenos_ct

left_join(phenos_ct, homozyg_SP70) %>%
  .[complete.cases(.),] ->
  phenos_CT_matched

names(phenos_CT_matched) = c("sampleid", "CTF", "CTM", "SP70")
wolb = fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/wolbachia.status.txt")
names(wolb)[1] = "sampleid"
left_join(phenos_CT_matched, wolb) -> phenos_CT_matched
save(phenos_CT_matched, file = "phenos_CT_matched.Rdata")

##### Bring in the PCA data

phenos_CT_matched  %>%
  left_join(pca_coords_annot) ->
  phenos_CT_matched_plus_PCA

phenos_CT_matched_plus_PCA %>% melt(id = c("sampleid",
                                           "SP70", "Infection_Status",
                                           "Dim.1","Dim.2","Dim.3",
                                           "Dim.4","Dim.5","pop")) ->
  phenos_CT_matched_plus_PCA.melt


phenos_CT_matched_plus_PCA.melt %>%
  ggplot(aes(
    x=SP70,
    y=value,
    fill=Infection_Status 
  )) + geom_boxplot() + 
  facet_grid(variable~Infection_Status) +
  theme_bw()->
  pheno_vio_all
ggsave(pheno_vio_all, file = "pheno_vio_all.pdf",
       w= 3, h = 4)

anova(lm(value ~ SP70*Infection_Status*variable, 
    data = phenos_CT_matched_plus_PCA.melt))

phenos_CT_matched %>%
  filter(Infection_Status == "n") %>%
  t.test(CTM~SP70, data = .)
phenos_CT_matched %>%
  filter(Infection_Status == "y") %>%
  t.test(CTF~SP70, data = .)

#### plot PCAs
phenos_CT_matched_plus_PCA.melt %>%
  ggplot(aes(
    x=Dim.1,
    y=value,
    color= variable  ,
    shape = SP70,
    linetype = SP70
  )) + geom_point() + geom_smooth(method = "lm", se = F) +
  facet_grid(Infection_Status~SP70)->
  pc_ctmax_plot

ggsave(pc_ctmax_plot, file = "pc_ctmax_plot.pdf",
       w = 4, h =3)

phenos_CT_matched_plus_PCA.melt %>%
  ggplot(aes(
    x=Dim.2,
    y=value,
    color= variable  ,
    shape = SP70,
    linetype = SP70
  )) + geom_point() + geom_smooth(method = "lm", se = F) +
  facet_grid(Infection_Status~SP70)->
  pc2_ctmax_plot

ggsave(pc2_ctmax_plot, file = "pc2_ctmax_plot.pdf",
       w = 4, h =3)


#### correlations
phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.1+CTF, data = .)
phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.1+CTM, data = .)


phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.2+CTF, data = .)
phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.2+CTM, data = .)

phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 2) %>%
  cor.test(~Dim.2+CTF, data = .)
phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 2) %>%
  cor.test(~Dim.2+CTM, data = .)

save(phenos_CT_matched_plus_PCA, file = "phenos_CT_matched_plus_PCA.Rdata")
############ Eliza and Olin's Data
############ Eliza and Olin's Data
############ Eliza and Olin's Data
############ Eliza and Olin's Data
############ Eliza and Olin's Data
############ Eliza and Olin's Data
load("phenos_CT_matched_plus_PCA.Rdata")

E_O_data <- fread("Eliza_Olins_data.txt")

inner_join(E_O_data, phenos_CT_matched_plus_PCA) ->
  EO_dat_Lecheta_dat

mean(filter(EO_dat_Lecheta_dat, Infection_Status == "n")$Proportion)
sd(filter(EO_dat_Lecheta_dat, Infection_Status == "n")$Proportion)
min(filter(EO_dat_Lecheta_dat, Infection_Status == "n")$Proportion)
max(filter(EO_dat_Lecheta_dat, Infection_Status == "n")$Proportion)


model <- glm(cbind(Hatched, (Sample_size-Hatched)) ~ CTM*CTF*SP70_allele, 
             data = filter(EO_dat_Lecheta_dat, Infection_Status == "n") , family= binomial)

anova(model)

filter(EO_dat_Lecheta_dat, Infection_Status == "n", 
       SP70_allele == "C/C") %>%
  cor.test(~Proportion+CTF, data = .)

filter(EO_dat_Lecheta_dat, Infection_Status == "n", 
       SP70_allele == "C/C") %>%
  cor.test(~Proportion+CTM, data = .)

filter(EO_dat_Lecheta_dat, Infection_Status == "n", 
       SP70_allele == "A/A") %>%
  cor.test(~Proportion+CTF, data = .)

filter(EO_dat_Lecheta_dat, Infection_Status == "n", 
       SP70_allele == "A/A") %>%
  cor.test(~Proportion+CTM, data = .)

### PCA correlations
filter(EO_dat_Lecheta_dat, Infection_Status == "n", 
       SP70_allele == "C/C") %>%
  cor.test(~Proportion+Dim.1, data = .)
filter(EO_dat_Lecheta_dat, Infection_Status == "n", 
       SP70_allele == "C/C") %>%
  cor.test(~Proportion+Dim.2, data = .)

filter(EO_dat_Lecheta_dat, Infection_Status == "n", 
       SP70_allele == "A/A") %>%
  cor.test(~Proportion+Dim.1, data = .)
filter(EO_dat_Lecheta_dat, Infection_Status == "n", 
       SP70_allele == "A/A") %>%
  cor.test(~Proportion+Dim.2, data = .)


#####
EO_dat_Lecheta_dat.wBinom <-
  foreach(i=1:dim(EO_dat_Lecheta_dat)[1],
          .combine = "rbind",
          .errorhandling = "remove")%do%{
            
            tmp <- EO_dat_Lecheta_dat[i,]
            
            test <- prop.test(x=tmp$Hatched, n=tmp$Sample_size, 
                              p =0.5)
            
            data.frame(tmp,
                       p.val=test$p.value,
                       uci=test$conf.int[1],
                       lci=test$conf.int[2]
            )
          }

####

EO_dat_Lecheta_dat.wBinom %>%
  filter(Infection_Status == "n") %>%
  dplyr::select(CTF,CTM,Proportion,lci,uci,SP70_allele) %>%
  melt(id = c("Proportion","lci","uci","SP70_allele")) %>%
  ggplot(aes(
    x=value,
    y=Proportion,
    ymin=lci,
    ymax=uci,
    color=variable,
  )) + 
  geom_errorbar(width = 0.1) +
  geom_point(size = 2) + geom_smooth(method = "lm", se = F) +
  facet_grid(~SP70_allele) + theme_bw() ->
  adultvsegg.hatch
ggsave(adultvsegg.hatch, file = "adultvsegg.hatch.pdf",
       w= 8, h = 2.0)

####
  
EO_dat_Lecheta_dat.wBinom %>%
  filter(Infection_Status == "n") %>%
  dplyr::select(Dim.1,Dim.2,Proportion,lci,uci,SP70_allele) %>%
  melt(id = c("Proportion","lci","uci","SP70_allele")) %>%
  ggplot(aes(
    x=value,
    y=Proportion,
    ymin=lci,
    ymax=uci,
    color=SP70_allele,
  )) + 
  geom_errorbar(width = 0.5) +
  geom_point() + geom_smooth(method = "lm", se = F) +
    facet_grid(~variable)->
  pcs.hatch
ggsave(pcs.hatch, file = "pcs.hatch.pdf",
       w= 6, h = 3)



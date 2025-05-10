## PLOT PCA
#Load libraries

##### THIS CODE HAS A CHECKPOINT FOR GRAPHING!!! se below... 

library(scales)
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
library(forcats)


win = 1600
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
  filter(pos > 20551633-win) %>%
  filter(pos < 20551633+win) ->
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
  filter(pos > 20551633-win) %>%
  filter(pos < 20551633+win) ->
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

save(joint_data, file = "joint_data.Rdata")

joint_data %>%
  filter(SNP_id.y == "2R_20551633_SNP") %>%
  t() %>% data.frame(SP70 = ., sampleid = rownames(.)) %>%
  filter(SP70 %in% c(0,2)) ->
  homozyg_SP70
save(homozyg_SP70, file = "homozyg_SP70.Rdata")

joint_data[,which(colnames(joint_data) %in% 
                    homozyg_SP70$sampleid)] ->
  dat_for_pca

dat_for_pca %>% t() ->
  dat_for_pca_t

colnames(dat_for_pca_t) = joint_data$SNP_id.y

save(dat_for_pca_t, file = "data.for.SP70_1600bp.pca_joint.Rdata")

dat_for_pca_t %>%
  PCA(graph = F) ->
  PCA_obj

###save
save(PCA_obj, file = "PCA_obj.haplo.Rdata")
load("PCA_obj.haplo.Rdata")
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

save(pca_coords_annot, file = "data.for.plotting.pca_coords_annot.Rdata")

##### START HERE!!!
##### START HERE!!!##### START HERE!!!##### START HERE!!!##### START HERE!!!
##### START HERE!!!
##### START HERE!!!
##### START HERE!!!

load("homozyg_SP70.Rdata")
load("data.for.plotting.pca_coords_annot.Rdata")

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

### Phenotyping 
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
#save(phenos_CT_matched, file = "phenos_CT_matched.Rdata")

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
  filter(Infection_Status == "n") -> 
  phenos_CT_matched_plus_PCA.melt.no_wolb

phenos_CT_matched_plus_PCA.melt.no_wolb %>%
  ggplot(aes(
    x=SP70,
    y=value,
  )) + geom_boxplot() + 
  facet_grid(.~variable) +
  theme_bw()->
  pheno_vio_all
ggsave(pheno_vio_all, file = "pheno_vio_all.pdf",
       w= 5, h = 4)

t.test(
  value ~ SP70,
  data = filter(phenos_CT_matched_plus_PCA.melt.no_wolb, 
                variable == "CTM")
)
t.test(
  value ~ SP70,
  data = filter(phenos_CT_matched_plus_PCA.melt.no_wolb, 
                variable == "CTF")
)

#### correlations
phenos_CT_matched_plus_PCA %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.1+CTF, data = .) -> F10
phenos_CT_matched_plus_PCA %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.2+CTF, data = .) -> F20

phenos_CT_matched_plus_PCA %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.1+CTM, data = .) -> M10
phenos_CT_matched_plus_PCA %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.2+CTM, data = .) -> M20

phenos_CT_matched_plus_PCA %>%
  filter(SP70 == 2) %>%
  cor.test(~Dim.1+CTF, data = .) -> F11
phenos_CT_matched_plus_PCA %>%
  filter(SP70 == 2) %>%
  cor.test(~Dim.2+CTF, data = .) -> F21

phenos_CT_matched_plus_PCA %>%
  filter(SP70 == 2) %>%
  cor.test(~Dim.1+CTM, data = .) -> M11
phenos_CT_matched_plus_PCA %>%
  filter(SP70 == 2) %>%
  cor.test(~Dim.2+CTM, data = .) -> M21
#
#data.frame(
#  type = rep(c(rep("F", 3), rep("M", 3)),2),
#  geno = c(rep(0, 6), rep(2, 6)),
#  pc = c(rep(1:3, 2), rep(1:3, 2)),
#  cor = c(F10$estimate,F20$estimate,F30$estimate,
#          M10$estimate,M20$estimate,M30$estimate,
#          F11$estimate,F21$estimate,F31$estimate,
#          M11$estimate,M21$estimate,M31$estimate
#          ),
#  P   = c(F10$p.value,F20$p.value,F30$p.value,
#          M10$p.value,M20$p.value,M30$p.value,
#          F11$p.value,F21$p.value,F31$p.value,
#          M11$p.value,M21$p.value,M31$p.value
#  )
#)
#
#####
#### plot PCAs
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 2
cols = gg_color_hue(n)


phenos_CT_matched_plus_PCA.melt %>%
  filter(Infection_Status == "n") %>%
  filter(variable == "CTF" & SP70 == 0) %>%
  ggplot(aes(
    x=Dim.1,
    y=value,
  )) + geom_point(color = cols[1], size = 2.5 ) + 
  geom_smooth(method = "lm", se = F, color = "black") +
  theme_bw() ->
  ctf_0_pc1
ggsave(ctf_0_pc1, file = "ctf_0_pc1.pdf",
       w = 4.5, h =4)

phenos_CT_matched_plus_PCA.melt %>%
  filter(Infection_Status == "n") %>%
  filter(variable == "CTM" & SP70 == 2) %>%
  ggplot(aes(
    x=Dim.1,
    y=value,
  )) + geom_point(color = cols[2], size = 2.5 ) + 
  geom_smooth(method = "lm", se = F, color = "black") +
  theme_bw() ->
  ctf_0_pc1
ggsave(ctf_0_pc1, file = "ctm_2_pc1.pdf",
       w = 4.5, h =4)

phenos_CT_matched_plus_PCA.melt %>%
  filter(Infection_Status == "n") %>%
  filter(variable == "CTM" & SP70 == 2) %>%
  ggplot(aes(
    x=Dim.1,
    y=value,
  )) + geom_point(color = cols[2], size = 2.5 ) + 
  geom_smooth(method = "lm", se = F, color = "black") +
  theme_bw() ->
  ctf_0_pc1
ggsave(ctf_0_pc1, file = "ctm_2_pc1.pdf",
       w = 4.5, h =4)

E_O_data <- fread("Eliza_Olins_data.txt")
names(E_O_data)[1] = "sampleid"
inner_join(E_O_data, phenos_CT_matched_plus_PCA) ->
  EO_dat_Lecheta_dat
EO_dat_Lecheta_dat %>%
  filter(Infection_Status == "n" & SP70 == 2) %>%
  ggplot(aes(
    x=Dim.1,
    y=Proportion,
    ymin=
    ymax=
    color=SP70
  )) + geom_point(size = 2.5, color = "brown" ) + 
  geom_smooth(method = "lm", se = F) +
  theme_bw() ->
  ehs_aa
ggsave(ehs_aa, file = "ehs_aa.pdf",
       w = 4.5, h =4)

phenos_CT_matched_plus_PCA.melt %>%
  filter(Infection_Status == "n") %>%
  filter(variable == "CTF" & SP70 == 0) %>%
  ggplot(aes(
    x=Dim.2,
    y=value,
  )) + geom_point(color = cols[1], size = 2.5 ) + 
  geom_smooth(method = "lm", se = F, color = "black") +
  theme_bw() ->
  ctf_0_pc2
ggsave(ctf_0_pc2, file = "ctf_0_pc2.pdf",
       w = 4.5, h =4)

phenos_CT_matched_plus_PCA.melt %>%
  filter(Infection_Status == "n") %>%
  filter(variable == "CTM" & SP70 == 0) %>%
  ggplot(aes(
    x=Dim.2,
    y=value,
  )) + geom_point(color = cols[2], size = 2.5 ) + 
  geom_smooth(method = "lm", se = F, color = "black") +
  theme_bw() ->
  ctm_0_pc2
ggsave(ctm_0_pc2, file = "ctm_0_pc2.pdf",
       w = 4.5, h =4)


#####
phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.1+CTF, data = .)
phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.2+CTF, data = .)
phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 2) %>%
  cor.test(~Dim.1+CTF, data = .)
phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 2) %>%
  cor.test(~Dim.2+CTF, data = .)


phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.1+CTM, data = .)
phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.2+CTM, data = .)
phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 2) %>%
  cor.test(~Dim.1+CTM, data = .)
phenos_CT_matched_plus_PCA %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 2) %>%
  cor.test(~Dim.2+CTM, data = .)

####
E_O_data <- fread("Eliza_Olins_data.txt")
names(E_O_data)[1] = "sampleid"
inner_join(E_O_data, phenos_CT_matched_plus_PCA) ->
  EO_dat_Lecheta_dat
EO_dat_Lecheta_dat %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.1+Proportion, data = .)
EO_dat_Lecheta_dat %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 0) %>%
  cor.test(~Dim.2+Proportion, data = .)

EO_dat_Lecheta_dat %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 2) %>%
  cor.test(~Dim.1+Proportion, data = .)
EO_dat_Lecheta_dat %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 2) %>%
  cor.test(~Dim.2+Proportion, data = .)


EO_dat_Lecheta_dat %>%
  filter(Infection_Status == "n") %>%
  filter(SP70 == 2) %>%
  ggplot(aes(
    x=Dim.1,
    y=Proportion
  )) + geom_point(color = cols[1]  ) + 
  geom_smooth(method = "lm", se = F, color = "black") +
  facet_grid(~SP70) + theme_bw() ->
  pc1_hatch_plot

ggsave(pc1_hatch_plot, file = "pc1_hatch_plot.pdf",
       w = 2, h =2)


save(phenos_CT_matched_plus_PCA, file = "phenos_CT_matched_plus_PCA.Rdata")

### Phenotypic space
phenos_CT_matched  %>%
  full_join(pca_coords_annot) ->
  phenos_CT_matched_plus_PCA_FULL

full_join(phenos_CT_matched_plus_PCA_FULL,  VA_meta)->
  phenos_CT_matched_plus_PCA_FULL

ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = phenos_CT_matched_plus_PCA_FULL,
             aes(
               x=Dim.1,
               y=Dim.2,
               shape = SP70
             ), alpha = 0.5, size = 2.5, fill = "grey") +
  geom_jitter(data = filter(phenos_CT_matched_plus_PCA_FULL, pop == "DGRP" & 
                             !is.na(CTM) & SP70 == 0
                             ),
             aes(
               x=Dim.1,
               y=Dim.2,
               fill=CTM,
               shape = SP70
             ), size = 4) +
  scale_fill_gradient2(midpoint = 39.9, high = "red", low = "blue") +
  scale_shape_manual(values = 21:22) +
  theme_bw() ->
  plottack1
ggsave(plottack1, file = "plottack1.pdf",
       h = 3, w = 7)

phenos_CT_matched_plus_PCA_FULL %>%
  group_by(pop) %>%
  summarise(
    m = ci(Dim.1)[1],
    uci = ci(Dim.1)[2],
    lci = ci(Dim.1)[3],
  ) %>%
  filter(!is.na(pop)) %>%
  ggplot(aes(
    x=pop,
    y=m,
    ymin=uci,
    ymax=lci
  )) + geom_errorbar(w = 0.01) +
    geom_point(size = 3) + theme_bw()->
  pop_est

ggsave(pop_est, file = "pop_est.pdf")


##

############ Eliza and Olin's Data
############ Eliza and Olin's Data
############ Eliza and Olin's Data
############ Eliza and Olin's Data
############ Eliza and Olin's Data
############ Eliza and Olin's Data

load("phenos_CT_matched_plus_PCA.Rdata")

E_O_data <- fread("Eliza_Olins_data.txt")
names(E_O_data)[1] = "sampleid"

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


load("homozyg_SP70.Rdata")
load("data.for.plotting.pca_coords_annot.Rdata")

full_join(pca_coords_annot, EO_dat_Lecheta_dat.wBinom) ->
  pca_coord_w_lecheta


ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = pca_coord_w_lecheta,
             aes(
               x=Dim.1,
               y=Dim.2,
               shape = SP70
             ), alpha = 0.5, size = 2.5, fill = "grey") +
  geom_point(data = filter(pca_coord_w_lecheta, pop == "DGRP" & 
                              !is.na(Proportion) 
  ),
  aes(
    x=Dim.1,
    y=Dim.2,
    fill=Proportion,
    shape = SP70
  ), size = 4) +
  scale_fill_gradient2(midpoint = 0.3, high = "red", low = "blue") +
  scale_shape_manual(values = 21:22) +
  theme_bw() ->
  plottack_hatch
ggsave(plottack_hatch, file = "plottack_hatch.pdf",
       h = 3, w = 7)

####

EO_dat_Lecheta_dat.wBinom %>%
  dplyr::select(CTF,CTM,Proportion,lci,uci,SP70_allele, Infection_Status) %>%
  melt(id = c("Proportion","lci","uci","SP70_allele","Infection_Status")) %>%
  ggplot(aes(
    x=value,
    y=Proportion,
    ymin=lci,
    ymax=uci,
    color=variable,
    #shape = Infection_Status
  )) + 
  geom_errorbar(width = 0.1) +
  geom_point(size = 2) + geom_smooth(method = "lm", se = F) + 
  facet_grid(.~SP70_allele) + 
  theme_bw() ->
  adultvsegg.hatch
ggsave(adultvsegg.hatch, file = "adultvsegg.hatch.pdf",
       w= 6, h = 2.5)

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


##### Where SNPs are
load("PCA_obj.haplo.Rdata")
dimdesc(PCA_obj, axes = 1:2, proba = 1.0) ->
  corr_study

corr_study$Dim.1 %>%
  as.data.frame() %>%
  mutate(p.adj = p.adjust(quanti.p.value, "bonferroni")) %>%
  mutate(SNPid = rownames(.))%>%
  mutate(PC = 1)->
  pc1_dim

corr_study$Dim.2 %>%
  as.data.frame() %>%
  mutate(p.adj = p.adjust(quanti.p.value, "bonferroni")) %>%
  mutate(SNPid = rownames(.)) %>%
  mutate(PC = 2)->
  pc2_dim

#####
pc1_dim %>%
  filter(abs(quanti.correlation) > 0.6) %>%
  mutate(SNPid= rownames(.)) %>%
  .$SNPid -> ids1
#####
pc2_dim %>%
  filter(abs(quanti.correlation) > 0.6) %>%
  mutate(SNPid= rownames(.)) %>%
  .$SNPid -> ids2

#####
load("joint_data.Rdata")
load("homozyg_SP70.Rdata")
load("phenos_CT_matched_plus_PCA.Rdata")

joint_data %>%
  filter(SNP_id.y == "2R_20550762_SNP") %>%
  t() %>% data.frame(SNP2 = ., sampleid = rownames(.)) %>%
  filter(SNP2 %in% c(0,2)) ->
  homozyg_target2

joint_data %>%
  filter(SNP_id.y %in% c(ids1, ids2)) %>%
  t() %>% as.data.frame() %>%
  mutate(sampleid = rownames(.)) %>%
  filter(grepl("line", sampleid)) ->
  other_snps_genos

colnames(other_snps_genos) = c(ids1, ids2, "sampleid")
other_snps_genos[complete.cases(other_snps_genos),] -> other_snps_genos
  

apply(other_snps_genos[,-which(names(other_snps_genos) == "sampleid")], 
      1, function(x) paste(x, collapse = ";") ) -> joint_geno_Calls
data.frame(joint_geno_Calls,sampleid=other_snps_genos$sampleid) ->
  homozyg_target_other

left_join(homozyg_SP70, homozyg_target_other) %>%
  filter(!is.na(joint_geno_Calls)) %>%
  mutate(joint_geno = paste(SP70,joint_geno_Calls,
                            sep = "_") ) %>%
  left_join(phenos_CT_matched_plus_PCA) %>%
  filter(!is.na(SP70)) %>%
  filter(Infection_Status == "n") ->
  data_other_SNPs

data_other_SNPs$joint_geno %>% 
  table %>% .[. >= 3] %>% names -> common_haps

data_other_SNPs %>%
  filter(joint_geno %in% common_haps) %>%
  ggplot(aes(
    x=fct_reorder(joint_geno, CTM),
    y=CTM
  )) + geom_boxplot() +
  coord_flip()  ->
  box_joint

ggsave(box_joint, file = "box_joint.pdf")

####


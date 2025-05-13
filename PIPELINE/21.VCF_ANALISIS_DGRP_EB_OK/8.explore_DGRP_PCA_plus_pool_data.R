#### Where our samples are

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

load("data.for.SP70_1800bp.pca_joint.Rdata")
load("homozyg_SP70.Rdata")
load("PCA_obj.haplo.Rdata")
wolb = fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/wolbachia.status.txt")

dat_for_pca_t %>% 
  .[which(grepl("line", rownames(.))),]  ->
  dat_for_pca_t2

dat_for_pca_t2 %>%
  PCA(graph = F #, 
      #ind.sup= which(rownames(dat_for_pca_t2) %in% rownames(dat))
  ) ->
  PCA_obj2

PCA_obj2$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleid = rownames(.)) ->
  main.pca

#####
### Explore Correlations
E_O_data <- fread("Eliza_Olins_data.txt")
names(E_O_data)[1] = "sampleid"
E_O_data <-
  foreach(i=1:dim(E_O_data)[1],
          .combine = "rbind",
          .errorhandling = "remove")%do%{
            
            tmp <- E_O_data[i,]
            
            test <- prop.test(x=tmp$Hatched, n=tmp$Sample_size, 
                              p =0.5)
            
            data.frame(tmp,
                       p.val=test$p.value,
                       uci=test$conf.int[1],
                       lci=test$conf.int[2]
            )
          }


phenos <- readRDS("wideform.fixed.phenotable.RDS")
#names(phenos)

MCT = "HighThermalToleranceExtreme_VaryingWithTemperature_F"
FCT = "HighThermalToleranceExtreme_VaryingWithTemperature_M"

data.frame(
  sampleid = paste("line",phenos$ral_id, sep = "_"),
  CT_F = phenos[,"HighThermalToleranceExtreme_VaryingWithTemperature_F"],
  CT_M = phenos[,"HighThermalToleranceExtreme_VaryingWithTemperature_M"]
) -> phenos_ct
names(phenos_ct) = c("sampleid", "CTF", "CTM")

names(wolb)[1] = "sampleid"
left_join(phenos_ct, homozyg_SP70) %>%
  left_join(wolb) %>%
  full_join(E_O_data[,c("sampleid","Proportion",
                        "uci", "lci")]) %>%
  filter()->
  ALL_PHENO_DATA_DGRP_SP70

full_join(main.pca, ALL_PHENO_DATA_DGRP_SP70) %>%
  filter(!is.na(SP70))->
  PCA_plus_ALL_PHENO_DATA_DGRP_SP70

save(PCA_plus_ALL_PHENO_DATA_DGRP_SP70, file = "PCA_plus_ALL_PHENO_DATA_DGRPonly_SP70.Rdata")
load("PCA_plus_ALL_PHENO_DATA_DGRPonly_SP70.Rdata")

PCA_plus_ALL_PHENO_DATA_DGRP_SP70 %>%
  filter(Infection_Status == "n") ->
  PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb
  
####
cor.test(~CTF+Dim.1,
         data = 
           filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  SP70 == 0))
cor.test(~CTF+Dim.1,
         data = 
           filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  SP70 == 2))
cor.test(~CTM+Dim.1,
         data = 
           filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  SP70 == 0))
cor.test(~CTM+Dim.1,
         data = 
           filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  SP70 == 2))



cor.test(~CTF+Dim.2,
         data = 
           filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  SP70 == 0))
cor.test(~CTF+Dim.2,
         data = 
           filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  SP70 == 2))
cor.test(~CTM+Dim.2,
         data = 
           filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  SP70 == 0))
cor.test(~CTM+Dim.2,
         data = 
           filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  SP70 == 2))



cor.test(~Proportion+Dim.1,
         data = 
           filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  SP70 == 0))
cor.test(~Proportion+Dim.2,
         data = 
           filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  SP70 == 0))

cor.test(~Proportion+Dim.1,
         data = 
           filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  SP70 == 2))
cor.test(~Proportion+Dim.2,
         data = 
           filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  SP70 == 2))
#######
#######
#######
#######
####### == plots

ggplot() +
  geom_point(
    data = PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
    aes(
      x=Dim.1,
      y=Dim.2,
      shape = SP70
    ), color = "grey"
  ) +
  geom_point(
    data = filter(PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb,
                  !is.na(CTF)),
    aes(
      x=Dim.1,
      y=Dim.2,
      fill=CTM,
      shape = SP70
    ), size = 3
  ) + theme_bw() +
  scale_shape_manual(values = 21:24) +
  scale_fill_gradient2(midpoint = 39.7, high = "red", low = "blue") ->
  DGRP_dim1_2_ctm
ggsave(DGRP_dim1_2_ctm, file = "DGRP_dim1_2_ctm.pdf",
       w = 5, h = 2.5)

#####
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 2
cols = gg_color_hue(n)

PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb %>%
  select(Dim.1, Dim.2, SP70, CTF, CTM, Proportion) %>%
  melt(id = c("Dim.1", "Dim.2", "SP70")) %>%
  filter(!is.na(value)) %>%
  as.data.frame() %>%
  ggplot(aes(
    x=Dim.1,
    y=value,
    shape=SP70,
    linetype = SP70,
    color=variable
  )) + geom_point() + theme_bw() + scale_shape_manual(values = c(1,0)) +
  scale_color_manual(values = c(cols,"brown")) +
  geom_smooth(method = "lm", se = F, color = "black") +
  facet_wrap(~variable, scales = "free_y")->
  pc1_vars
ggsave(pc1_vars, file = "pc1_vars.pdf",
       w= 6, h = 2)

PCA_plus_ALL_PHENO_DATA_DGRP_SP70_NoWolb %>%
  select(Dim.1, Dim.2, SP70, CTF, CTM, Proportion) %>%
  melt(id = c("Dim.1", "Dim.2", "SP70")) %>%
  as.data.frame() %>%
  filter(!is.na(value)) %>%
  ggplot(aes(
    x=Dim.2,
    y=value,
    shape=SP70,
    linetype = SP70,
    color=variable
  )) + geom_point() + theme_bw() + scale_shape_manual(values = c(1,0)) +
  scale_color_manual(values = c(cols,"brown")) +
  geom_smooth(method = "lm", se = F, color = "black") +
  facet_wrap(~variable, scales = "free_y")->
  pc2_vars
ggsave(pc2_vars, file = "pc2_vars.pdf",
       w= 6, h = 2)



###### supplementary data


#####
PCA_obj$var %>%
  as.data.frame %>% rownames(.) ->
  snps_of_interest

####

samps <- fread("/netfiles02/lockwood_lab/IntrogressionProject/Population_files/samp.guide.files.introg.txt")
samps %>% 
  filter(Pop %in% c("SK","VT8",
                    "SKF","VT8F",
                    "SKF2","VT8F2",
                    "SKF3","VT8F3"#,
                    #"ME_temp","PAN_trop"
  )) %>%
  mutate(type = case_when(
    #Pop %in% c("ME_temp", "PAN_trop") ~ "wild",
    Pop %in% c("SK") ~ "SK-parent",
    Pop %in% c("VT8") ~ "VT8-parent",
    Pop %in% c("VT8F", "VT8F2","VT8F3") ~ "F16-VT-mother",
    Pop %in% c("SKF", "SKF2","SKF3") ~ "F16-SK-mother"
  ))  -> samps.flt

#####
genofile.path <- "/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds"
genofile <- seqOpen(genofile.path)

seqResetFilter(genofile)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile)) %>%
  mutate(SNP_id = paste(chr, pos, "SNP", sep = "_"))

snps.dt <- snps.dt[nAlleles==2]
seqSetFilter(genofile, variant.id=snps.dt$variant.id, sample.id = samps.flt$Sample_id)
snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

seqSetFilter(genofile,
             snps.dt[SNP_id%in%c(
               snps_of_interest
             )][missing<.05]$variant.id)

### get allele frequency data
ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")

dat <- ad$data/dp
dim(dat)  

colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), "SNP", sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")

#####
dim(dat)
dim(dat_for_pca_t)
dat_for_pca_t_cor = dat_for_pca_t/2

SNPs_taregt <- colnames(dat_for_pca_t)[which(colnames(dat_for_pca_t) %in% colnames(dat))]

dat_for_pca_t_cor[,SNPs_taregt] %>% 
  .[which(grepl("line", rownames(.))),] %>%
  rbind(dat[,SNPs_taregt]) ->
  dat_for_pca_t2

dat_for_pca_t2 %>%
  PCA(graph = F , 
      ind.sup= which(rownames(dat_for_pca_t2) %in% rownames(dat))
  ) ->
  PCA_obj2

PCA_obj2$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleid = rownames(.)) ->
  main.pca

PCA_obj2$ind.sup$coord %>%
  as.data.frame() %>%
  mutate(sampleid = rownames(.)) %>%
  separate(remove = F,
           sampleid,
           into = c("parent", "rep", " tech"),
           sep = "_") %>%
  separate(parent, into = c("pop","etc"), sep = 2)->
  supp.pca

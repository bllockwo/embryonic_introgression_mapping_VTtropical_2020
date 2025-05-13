
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


E_O_data <- fread("Eliza_Olins_data.txt")
names(E_O_data)[1] = "sampleid"

##### Where SNPs are
##### Where SNPs are##### Where SNPs are
##### Where SNPs are
##### Where SNPs are
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
  filter(abs(quanti.correlation) > 0.59) %>%
  mutate(SNPid= rownames(.)) %>%
  .$SNPid -> ids1
#####
pc2_dim %>%
  filter(abs(quanti.correlation) > 0.59) %>%
  mutate(SNPid= rownames(.)) %>%
  .$SNPid -> ids2

top_snps = c("2R_20550549_SNP",
             "2R_20551182_SNP",
             "2R_20551519_SNP",
             "2R_20550323_SNP",
             "2R_20552286_SNP",
             "2R_20551486_SNP",
             "2R_20551129_SNP",
             "2R_20550070_SNP",
             "2R_20550046_SNP",
             "2R_20550403_SNP")

c(ids1, ids2) %in% top_snps

pc1_dim %>%
  filter(SNPid %in% c("2R_20550549_SNP",
                      "2R_20551182_SNP",
                      "2R_20551519_SNP",
                      "2R_20550323_SNP",
                      "2R_20552286_SNP",
                      "2R_20551486_SNP",
                      "2R_20551129_SNP",
                      "2R_20550070_SNP",
                      "2R_20550046_SNP",
                      "2R_20550403_SNP"))

pc2_dim %>%
  filter(SNPid %in% c("2R_20550549_SNP",
                      "2R_20551182_SNP",
                      "2R_20551519_SNP",
                      "2R_20550323_SNP",
                      "2R_20552286_SNP",
                      "2R_20551486_SNP",
                      "2R_20551129_SNP",
                      "2R_20550070_SNP",
                      "2R_20550046_SNP",
                      "2R_20550403_SNP"))



#####
load("joint_data.Rdata")
load("homozyg_SP70.Rdata")
load("phenos_CT_matched_plus_PCA.Rdata")

#joint_data %>%
#  filter(SNP_id.y == "2R_20550762_SNP") %>%
#  t() %>% data.frame(SNP2 = ., sampleid = rownames(.)) %>%
#  filter(SNP2 %in% c(0,2)) ->
#  homozyg_target2

joint_data %>%
  filter(SNP_id.y %in% top_snps) ->
  top_snps_data

top_snps_data %>%
  t() %>% as.data.frame() %>%
  mutate(sampleid = rownames(.)) %>%
  filter(grepl("line", sampleid)) ->
  other_snps_genos

colnames(other_snps_genos) = c(top_snps, "sampleid")
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
  table %>% table

data_other_SNPs$joint_geno %>% 
  table %>% .[. >= 3] %>% names -> common_haps

data_other_SNPs %>%
  filter(joint_geno %in% common_haps) %>%
  select(joint_geno, CTF,CTM) %>%
  melt(id = c("joint_geno")) ->
  dat.for.plot

t.test(value~joint_geno,
       data = filter(dat.for.plot, variable == "CTF" &
                       joint_geno %in% 
                       c("0_0;0;0;0;0;0;0;0;0;0",
                         "2_0;0;0;0;0;0;0;0;0;0")))
t.test(value~joint_geno,
       data = filter(dat.for.plot, variable == "CTM" &
                       joint_geno %in% 
                       c("0_0;0;0;0;0;0;0;0;0;0",
                         "2_0;0;0;0;0;0;0;0;0;0")))

## rescue pheno
t.test(value~joint_geno,
       data = filter(dat.for.plot, variable == "CTF" &
                       joint_geno %in% 
                       c("2_2;2;2;2;0;2;0;2;2;2",
                         "2_0;0;0;0;0;0;0;0;0;0")))


dat.for.plot %>%
  ggplot(aes(
    x=joint_geno,
    y=value,
    fill=variable,
  )) + geom_boxplot() + theme_classic() +
  facet_grid(~variable)  ->
  box_joint

ggsave(box_joint, file = "box_joint.pdf",
       w= 3, h = 3)


data_other_SNPs %>%
  separate(joint_geno_Calls, into = paste("snp", 1:10, sep = "_"), sep = ";" ) %>%
  mutate(SNP_distance = (as.numeric(snp_1) +
                           as.numeric(snp_2) +
                           as.numeric(snp_3) +
                           as.numeric(snp_4) +
                           as.numeric(snp_5) +
                           as.numeric(snp_6) +
                           as.numeric(snp_7) +
                           as.numeric(snp_8) +
                           as.numeric(snp_9) +
                           as.numeric(snp_10))) %>%
  mutate(SNP_distance_rel = SNP_distance/2) -> 
  distance_analysis

anova(lm(CTF~SNP_distance_rel*SP70, data = distance_analysis))
anova(lm(CTM~SNP_distance_rel*SP70, data = distance_analysis))

distance_analysis %>%
  ggplot(aes(
    x=SNP_distance_rel,
    y=CTF,
    color=SP70
  )) + geom_point(size = 2) + 
  scale_color_brewer(palette = "Dark2") +
  geom_smooth(method = "lm", se =F) + theme_bw() ->
  snp_dist_plot

ggsave(snp_dist_plot, file = "snp_dist_plot.pdf",
       w= 3.5, h = 3)

#### Embryos....
#### Embryos....
#### Embryos....
#### Embryos....

left_join(homozyg_SP70, homozyg_target_other) %>%
  filter(!is.na(joint_geno_Calls)) %>%
  mutate(joint_geno = paste(SP70,joint_geno_Calls,
                            sep = "_")) %>%
  left_join(E_O_data) %>%
  filter(!is.na(Proportion)) %>%
  filter(Wolbachia_Status == "No") ->
  data_other_SNPs_Embrios

data_other_SNPs_Embrios %>%
  ggplot(aes(
    x=joint_geno,
    y=Proportion,
  )) + geom_boxplot() + 
  theme_classic()  ->
  box_joint_Embrios

ggsave(box_joint_Embrios, file = "box_joint_Embrios.pdf",
       w= 3, h = 3)



data_other_SNPs_Embrios %>%
  separate(joint_geno_Calls, into = paste("snp", 1:10, sep = "_"), sep = ";" ) %>%
  mutate(SNP_distance = (as.numeric(snp_1) +
                           as.numeric(snp_2) +
                           as.numeric(snp_3) +
                           as.numeric(snp_4) +
                           as.numeric(snp_5) +
                           as.numeric(snp_6) +
                           as.numeric(snp_7) +
                           as.numeric(snp_8) +
                           as.numeric(snp_9) +
                           as.numeric(snp_10))) %>%
  mutate(SNP_distance_rel = SNP_distance/2) %>%
  ggplot(aes(
    x=SNP_distance_rel,
    y=Proportion,
    color=SP70
  )) + geom_point() + geom_smooth(method = "lm") ->
  snp_dist_plot_Emb

ggsave(snp_dist_plot_Emb, file = "snp_dist_plot_Emb.pdf",
       w= 6, h = 6)

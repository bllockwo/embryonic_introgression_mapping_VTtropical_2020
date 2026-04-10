library(tidyverse)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(adegenet)
library(reshape2)
library(FactoMineR)
require(gtools)
require(foreach)

genofile <- seqOpen("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/dgrp2.gds.gz")
seqResetFilter(genofile)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      alleles=seqGetData(genofile, "allele"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile)) %>%
                      mutate(snp_id = paste(chr,pos, sep = "_"))

snps.dt %>% 
filter(snp_id %in%
c(
"2R_16439019",
"2R_16439138",
"X_15496974",
"X_15501637",
"X_15741847")) -> snps.sel

seqSetFilter(genofile,  
             variant.id=snps.sel$variant.id)                      

allelic.data <- data.table(
					  chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile),
                      alleles=seqGetData(genofile, "allele")
                      ) %>%
                      mutate(snp_id = paste(chr,pos, sep = "_"))


snpgdsGetGeno(genofile, snp.id=snps.sel$variant.id,
              verbose=TRUE, with.id=TRUE) ->
  GENO.matrix


GENO.matrix$genotype %>%
  as.data.frame() ->
  GENO.matrix.dt

colnames(GENO.matrix.dt) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), sep="_")
rownames(GENO.matrix.dt) <- seqGetData(genofile, "sample.id")
### load data
egg_data <- fread("/gpfs2/scratch/jcnunez/Embryo_fix/OK_EB_data.Aug16.2026.txt")
all_phenos <- readRDS("/gpfs2/scratch/jcnunez/fst_brent/GWAS_EB_OK/wideform.fixed.phenotable.RDS")

egg_data.wBinom <-
  foreach(i=1:dim(egg_data)[1],
          .combine = "rbind",
          .errorhandling = "remove")%do%{
            
            tmp <- egg_data[i,]
            
            test <- prop.test(x=tmp$Hatched, n=tmp$Sample_size, 
                              p =0.5)
            
            data.frame(tmp,
                       p.val=test$p.value,
                       uci=test$conf.int[1],
                       lci=test$conf.int[2]
            )
          }


###
wolb <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/wolbachia.status.txt")
inv <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/Inversion.status.txt")
names(wolb)[1] = "Line"
names(inv)[1] = "Line"


###
all_phenos %>% 
mutate(Line = paste("line_", ral_id, sep = "")) -> CTdats

####X--15,607,604 in dm6 is
#### -------> 15501637

GENO.matrix.dt %>%
mutate(Line = rownames(.)) %>%
full_join(egg_data.wBinom) %>% 
  full_join(CTdats) %>% 
  full_join(wolb) %>% 
  full_join(inv) ->
GENO.matrix.pheno


GENO.matrix.pheno %>%
  filter(!is.na(Proportion)) %>%
  mutate(gen_type = case_when(
    `2R_16439138` == 0 & `X_15501637` == 2 ~ "tropical",
    `2R_16439138` == 2 & `X_15501637` == 0 ~ "temperate",
    TRUE ~ "mixed"
  )) %>%
  ggplot(
    aes(
      x=fct_reorder(Line, Proportion),
      ymin=lci,
      ymax=uci,
      y=Proportion,
      fill = gen_type
    )
  ) + geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_errorbar(size = 0.1) + geom_point(shape = 21, size = 2) + theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values = c("grey","blue","red"))->
  egg_bimon_plot
ggsave(egg_bimon_plot, file = "egg_bimon_plot.pdf", h =2.5, w =1.5)

#####
####
GENO.matrix.pheno %>%
filter(!is.na(Proportion)) %>%
  filter(!is.na(`X_15501637`)) %>%
filter(!is.na(`2R_16439138`)) ->
clean.dat

clean.dat %>%
  group_by(`2R_16439138`,`X_15501637`) %>%
  summarize(N=n())

clean.dat %>%
  group_by(`2R_16439138`) %>%
  summarize(N=n())

clean.dat %>%
  group_by(`X_15501637`) %>%
  summarize(N=n())

### Plot
clean.dat %>%  
  ggplot(aes(
    x=paste( `2R_16439138`),
    y=Proportion,
  )) + geom_boxplot(width = 0.5) +
  ylim(0, 0.8) +
  ggtitle("2R Embryonic Thermal") + theme_bw() ->
  EmbRes2R

clean.dat %>%  
  ggplot(aes(
    x=paste(`X_15501637`),
    y=Proportion,
  )) + geom_boxplot(width = 0.5) +
  ylim(0, 0.8) +
  ggtitle("X Embryonic Thermal") + theme_bw()->
  EmbResX

### Plot both
clean.dat %>%  
ggplot(aes(
x=paste( `2R_16439138`,`X_15501637`),
y=Proportion,
)) + geom_boxplot(width = 0.5) +
  ylim(0, 0.8) +
ggtitle("Embryonic Thermal") + theme_bw()->
EmbRes

ggsave(EmbRes, file = "EmbRes.Box.evidence.pdf", w = 4, h = 3)
ggsave(EmbRes2R, file = "EmbRes2R.Box.evidence.pdf", w = 2, h = 3)
ggsave(EmbResX, file = "EmbResX.Box.evidence.pdf", w = 2, h = 3)


###
anovaFun <- function(m1, m2) {
  #  m1 <- t0; m2<- t1.year
  ll1 <- as.numeric(logLik(m1))
  ll2 <- as.numeric(logLik(m2))
  
  parameter <- abs(attr(logLik(m1), "df") -  attr(logLik(m2), "df"))
  
  chisq <- -2*(ll1-ll2)
  
  1-pchisq(chisq, parameter)
  
}

names(CTdats)[-which(names(CTdats) %in% c("Line","ral_id"))] -> phenos_allnams

pheno_scan=
foreach(i=phenos_allnams, .combine = "rbind", .errorhandling = "remove")%do%{
  
  clean.dat %>%
    select(i, `2R_16439138`, `X_15501637`, `In(2L)t`, `In(2R)NS`, Infection_Status ) ->
    tmp
  names(tmp)[1] = "var"
  lm(var~`2R_16439138`*`X_15501637`+`In(2L)t`+`In(2R)NS`+
       Infection_Status, data = tmp) -> res_inter
  lm(var~`In(2L)t`+`In(2R)NS`+
       Infection_Status, data = tmp) -> res_null
  
  data.table(
    i=i,
    p=anovaFun(res_null, res_inter)
  )
}


pheno_scan %>% filter(p < 0.05)

clean.dat %>%  
  ggplot(aes(
    x=paste( `2R_16439138`,`X_15501637`),
    y=`DayBoutNumber_standard_Female`,
  )) + geom_boxplot(width = 0.5) +
  ggtitle("Other") + theme_bw()->
  Other

ggsave(Other, file = "Other.pdf", w = 4, h = 3)

lm(Proportion~`2R_16439138`*`X_15501637`+`In(2L)t`+`In(2R)NS`+
     Infection_Status, data = clean.dat) -> res_inter

anova(res_inter)
##


###
### Adults

GENO.matrix.pheno %>%
  filter(!is.na(HighThermalToleranceExtreme_VaryingWithTemperature_F)) %>%
  filter(!is.na(X_15501637)) %>%
  filter(!is.na(`2R_16439138`)) %>%
  select(HighThermalToleranceExtreme_VaryingWithTemperature_F,
         HighThermalToleranceExtreme_VaryingWithTemperature_M,
         X_15501637,`2R_16439138`) %>%
  melt(id =c("X_15501637","2R_16439138")) %>%
  mutate(sex = case_when(variable == "HighThermalToleranceExtreme_VaryingWithTemperature_F" ~ "F",
                         variable == "HighThermalToleranceExtreme_VaryingWithTemperature_M" ~ "M",
                         ))->
  clean.dat.adult
  
clean.dat.adult %>%  
  ggplot(aes(
    x=paste( `2R_16439138`),
    y=value,
    fill = sex
  )) + geom_boxplot(width = 0.5) +
  ggtitle("2R adult Thermal") + theme_bw()  + theme(legend.position = "none")  ->
  AdRes2R

clean.dat.adult %>%  
  ggplot(aes(
    x=paste( `X_15501637`),
    y=value,
    fill = sex
  )) + geom_boxplot(width = 0.5) +
  ggtitle("X adult Thermal") + theme_bw() + theme(legend.position = "none") ->
  AdResX


### Plot both
clean.dat.adult %>%  
  ggplot(aes(
    x=paste( `2R_16439138`,`X_15501637`),
    y=value,
    fill = sex
  )) + geom_boxplot(width = 0.5) +
  ggtitle("adult Thermal") + theme_bw()  + theme(legend.position = "none") ->
  AdRes

ggsave(AdRes, file = "AdRes.Box.evidence.pdf", w = 4, h = 3)
ggsave(AdResX, file = "AdResX.Box.evidence.pdf", w = 2, h = 3)
ggsave(AdRes2R, file = "AdRes2R.Box.evidence.pdf", w = 2, h = 3)


### cor
GENO.matrix.pheno %>%
  filter(!is.na(X_15501637)) %>%
  filter(!is.na(`2R_16439138`)) %>%
ggplot(
  aes(
    x=HighThermalToleranceExtreme_VaryingWithTemperature_F,
    y=Proportion,
    color = as.character(`X_15501637`)
  )
) + geom_point() + geom_smooth(method = "lm", color = "black") + facet_grid(~`2R_16439138`) ->
  cor_data

ggsave(cor_data, file = "cor_data.pdf", w = 6, h= 2)

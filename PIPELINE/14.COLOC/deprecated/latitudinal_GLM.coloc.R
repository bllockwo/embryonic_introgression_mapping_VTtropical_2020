library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)
library(foreach)
library(data.table)
library(gmodels)
library(GenomicRanges)
library(rtracklayer)
library(SeqArray)
library(permute)
library(forcats)
nPerm=100

#### FUNCTIONS
getData <- function(variant, samps2use=seasonal.sets$sampleId) {
  # variant=snp.dt[chr=="2L"][pos==14617051]$variant.id;
  # samps2use=seasonal.sets$sampleId
  # variant=595225
  
  ### filter to target
  seqResetFilter(genofile)
  seqSetFilter(genofile, variant.id=variant, sample.id=samps2use)
  
  ### get frequencies
  message("Allele Freqs")
  
  ad <- seqGetData(genofile, "annotation/format/AD")$data
  dp <- seqGetData(genofile, "annotation/format/DP")
  
  af <- data.table(ad=expand.grid(ad)[,1],
                   dp=expand.grid(dp)[,1],
                   sampleId=rep(seqGetData(genofile, "sample.id"),
                                dim(ad)[2]),
                   variant.id=rep(seqGetData(genofile, "variant.id"),
                                  each=dim(ad)[1]),
                   chr=seqGetData(genofile, "chromosome"),
                   pos=seqGetData(genofile, "position"))
  
  ### tack them together
  af[,af:=ad/dp]
  
  ### calculate effective read-depth
  afis <- merge(af, samps[,c("sampleId", "nFlies"), with=F], by="sampleId")
  
  afis[chr=="X", nEff:=round((dp*nFlies)/(dp+nFlies-1))]
  afis[chr!="X", nEff:=round((dp*2*nFlies )/(dp+2*nFlies- 1))]
  afis[,af_nEff:=round(af*nEff)/nEff]
  
  ### return
  afis
}

### ------> BEGIN ANALYSIS

# needs CURL to install
#install.packages("./rtracklayer_1.64.0.tar.gz", repos = NULL, type="source")

##### data FET
real_dat <-"/netfiles02/lockwood_lab/IntrogressionProject/FET_output/Fisher_test_results_All_crosses.txt"
real_FET <- fread(real_dat)
names(real_FET)[1] = "chr"

real_FET %>%
  filter(sk_all_f.test_pval < 5e-08) ->
  boinfe_cutoff
boinfe_cutoff %>%
  filter( (chr == "2R" & POS > 19115753 & POS < 20615753) | (chr == "X" & POS > 15123410 & POS < 16123410) ) %>%
  dplyr::select(CHR=chr,   POS_DM6=POS) ->
  windows_interest

windows_interest_for_DEST <- windows_interest
names(windows_interest_for_DEST) = c("chr","pos")

windows_interest %>%
  filter(CHR == "2R") ->
  interest_2R

windows_interest %>%
  filter(CHR == "X") ->
  interest_X

#####

### open GDS file
genofile <- seqOpen("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2021.DEST.v1.release/dest.all.PoolSNP.001.50.10Nov2020.ann.gds")

### get target populations
samps <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2021.DEST.v1.release/samps_10Nov2020.csv")

samps <- rbind(samps[set=="DrosRTEC"],
               samps[set=="DrosEU"],
               samps[set=="dgn"]
)


### get subsample of data to work on
seqResetFilter(genofile)
seqSetFilter(genofile, sample.id=samps$sampleId)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))

snps.dt %>%
filter( (chr == "2R" & pos > 19115753 & pos < 20615753) | 
          (chr == "X" & pos > 15123410 & pos < 16123410)) ->
  snps.dt.flt
  
left_join(windows_interest_for_DEST, snps.dt.flt) %>%
  filter(!is.na(variant.id)) ->
  snps.dt.flt2

## choose number of alleles
snps.dt.flt2 <- snps.dt.flt2[nAlleles==2]
seqSetFilter(genofile, sample.id=samps$sampleId, 
             variant.id=snps.dt.flt2$variant.id)

snps.dt.flt2[,af:=seqGetData(genofile, "annotation/info/AF")$data]

### select sites
seqSetFilter(genofile, sample.id=samps$sampleId,
             snps.dt.flt2[chr%in%c("2R", "X")][missing<.05][af>.2]$variant.id)

### get allele frequency data
ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")

dat <- ad$data/dp
dim(dat)  

## Add metadata
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")

####
samps %>%
  filter(continent == "NorthAmerica") %>%
  .$sampleId -> N_Ame_samps

samps %>%
  filter(continent == "Europe") %>%
  filter(set %in% c("DrosRTEC","DrosEU") ) %>%
  .$sampleId

##########

message("Permutations")
o <- foreach(i=1:dim(windows_interest)[1], 
             .combine="rbind", .errorhandling="remove")%do%{
  #i <- 1; tmp.ids <- 786996
  message(paste(i, dim(windows_interest)[1], sep=" / "))
  
  ### get allele frequency data
  af <- getData(variant=snps.dt[chr==windows_interest$CHR[i]][pos==windows_interest$POS_DM6[i]]$variant.id,
                samps2use =N_Ame_samps)
  af <- merge(af, samps, by="sampleId")
  af[,year_pop:=as.factor(interaction(locality, year))]
  af <- af[!is.na(af_nEff)]
  af <- af[af_nEff>0 & af_nEff<1]
  ## af <- getData(variant=2178993)
  
  o <- foreach(j=1:nPerm, 
               .combine="rbind", 
               .errorhandling="remove")%do%{
    # j<-0
    
                 if(j != 0){
                   tmp_perm <- af
                   p_lrt=-999
                   seas.AIC = -999
                   null.AIC = -999
                   
                   ###
                   tmp_perm %<>%
                     mutate(
                       lat_perm = lat[shuffle(lat)]
                     )
                   ###
                   
                   null.perm <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ 1,    data = tmp_perm, family= binomial)
                   effect.perm <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ lat_perm, data = tmp_perm, family= binomial)
                   t3.sing <- t4.sing <- NA
                   
                   p_lrt=anova(effect.perm, null.perm, test="Chisq")[2,5]
                   effect.AIC = extractAIC(effect.perm)[2]
                   null.AIC = extractAIC(null.perm)[2]
                   
                   data.frame(
                     chr=windows_interest$CHR[i],
                     pos=windows_interest$POS_DM6[i],
                     variant.id=snps.dt[chr==windows_interest$CHR[i]][pos==windows_interest$POS_DM6[i]]$variant.id, 
                     perm=j,
                     b_lat=summary(effect.perm)$coef[2,1], se_temp=summary(effect.perm)$coef[2,2],
                     nObs=dim(af)[1],
                     nFixed=sum(af$af_nEff==0) + sum(af$af_nEff==1),
                     af=mean(af$af_nEff), neff=mean(af$nEff),
                     p_lrt=p_lrt,
                     seas.AIC = effect.AIC,
                     null.AIC = null.AIC
                   ) ### close deliverable
                   
                 } ### close if j !=0
                 
  } ## perms
  
} ## grand loop

save(o, file = "LAT_GLMs_Brent.PERMs.Rdata")

### REAL
message("real")
o_real <- foreach(i=1:dim(windows_interest)[1], 
                  .combine="rbind", 
                  .errorhandling="remove")%do%{
                    #i <- 1; tmp.ids <- 786996
                    message(paste(i, dim(windows_interest)[1], sep=" / "))
                    
                    ### get allele frequency data
                    af <- getData(variant=snps.dt[chr==windows_interest$CHR[i]][pos==windows_interest$POS_DM6[i]]$variant.id,
                                  samps2use =N_Ame_samps)
                    af <- merge(af, samps, by="sampleId")
                    af[,year_pop:=as.factor(interaction(locality, year))]
                    af <- af[!is.na(af_nEff)]
                    af <- af[af_nEff>0 & af_nEff<1]
                    ## af <- getData(variant=2178993)
                    j =0
                    
                                     tmp <- af
                                     p_lrt=-999
                                     seas.AIC = -999
                                     null.AIC = -999
                                     
                                     null.real <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ 1,    data = tmp, family= binomial)
                                     effect.real <- glm(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ lat, data = tmp, family= binomial)
                                     t3.sing <- t4.sing <- NA
                                     
                                     p_lrt=anova(effect.real, null.real, test="Chisq")[2,5]
                                     effect.AIC = extractAIC(effect.real)[2]
                                     null.AIC = extractAIC(null.real)[2]
                                     
                                     data.frame(
                                       chr=windows_interest$CHR[i],
                                       pos=windows_interest$POS_DM6[i],
                                       variant.id=snps.dt[chr==windows_interest$CHR[i]][pos==windows_interest$POS_DM6[i]]$variant.id, 
                                       perm=j,
                                       b_lat=summary(effect.real)$coef[2,1], se_temp=summary(effect.real)$coef[2,2],
                                       nObs=dim(af)[1],
                                       nFixed=sum(af$af_nEff==0) + sum(af$af_nEff==1),
                                       af=mean(af$af_nEff), neff=mean(af$nEff),
                                       p_lrt=p_lrt,
                                       seas.AIC = effect.AIC,
                                       null.AIC = null.AIC
                                     ) ### close deliverable

                                 } ## grand loop

save(o_real, file = "LAT_GLMs_Brent.REAL.Rdata")


####### CHECKPOINT!!
####### CHECKPOINT!!
####### CHECKPOINT!!
####### CHECKPOINT!!
####### CHECKPOINT!!

load("LAT_GLMs_Brent.PERMs.Rdata")
load("LAT_GLMs_Brent.REAL.Rdata")

#### sort
o_real %>%
  filter(chr == "2R") %>%
  arrange(log10(p_lrt)) %>%
  mutate(Rank = 1:dim(.)[1]) ->
  chr2R_o

o_real %>%
  filter(chr == "X") %>%
  arrange(log10(p_lrt)) %>%
  mutate(Rank = 1:dim(.)[1]) ->
  chrX_o

####
rbind(chr2R_o, chrX_o) -> o_real_mod

###
o %>%
  group_by(chr, pos) %>%
  summarise(uci=quantile(p_lrt, 0.05)) ->
  o_sum

###
left_join(o_real_mod, o_sum) ->
  o_real_mod_perm

#####
o_real_mod_perm %<>%
  mutate(lat_signif = p_lrt < uci)

#### ANNOTATIONS <_____ WAIT... THERE IS A CHECKPOINT
getANNOTS <- function(variant) {
  # variant=snps.dt[chr=="2L"][pos==14617051]$variant.id;
  # samps2use=seasonal.sets$sampleId
  # variant=595225
  
  ### filter to target
  seqResetFilter(genofile)
  seqSetFilter(genofile, variant.id=variant)
  
  ### get frequencies
  message("Annotations")
  
  annotation_data=seqGetData(genofile,"annotation/info/ANN")
  
  str_split(annotation_data$data, "\\|") -> annot_list
  
  annot_df = 
  foreach(i=1:length(annot_list), .combine = "rbind",
          .errorhandling = "remove")%do%{
            tmp <- annot_list[[i]]
            
            data.frame(
              SNP=tmp[1],
              effect=tmp[2],
              impact=tmp[3],
              gene_name=tmp[4],
              fbgn=tmp[5],
              type=tmp[8],
              DNA_mut=tmp[10],
              AA_mut=tmp[11]
            ) %>%
              mutate(Ann_rank = case_when(impact == "HIGH" ~ 1,
                                          impact == "MODERATE" ~ 2,
                                          impact == "LOW" ~ 3,
                                          impact == "MODIFIER" ~ 4,
                                          ))
          }
  
  data.table(chr=seqGetData(genofile, "chromosome"),
             pos=seqGetData(genofile, "position"),
             variant.id=seqGetData(genofile, "variant.id"),
             annot_df
             )
  
}

#### Explore significant SNPS
heat_snps_annot =
  foreach(i=1:dim(o_real_mod_perm)[1],
          .combine = "rbind",
          .errorhandling = "remove")%do%{
            getANNOTS(variant=snps.dt[chr==o_real_mod_perm$chr[i]][pos==o_real_mod_perm$pos[i]]$variant.id)  
          }


heat_snps_annot %>%
  group_by( chr, pos) %>%
  slice_max(Ann_rank, with_ties = F)  %>%
  mutate(effect_simple = case_when(
    effect %in% c("downstream_gene_variant","upstream_gene_variant") ~ "intergenic_region",
    TRUE ~ effect
  )) ->
  heat_snps_annot_simple

save(heat_snps_annot_simple, file = "heat_snps_annot_simple.Rdata")

#### CHECKPOINT PT2
load("heat_snps_annot_simple.Rdata")

####
left_join(o_real_mod_perm, heat_snps_annot_simple) ->
  o_real_mod_perm_annot

#### Annotations ^^^
o_real_mod_perm_annot %>%
  filter(lat_signif == T) %>%
  group_by(chr) %>%
  arrange(p_lrt) %>%
  slice_head(n=10) %>%
  as.data.frame()



#### PLOT
ggplot(data=filter(o_real_mod_perm_annot, effect_simple != "non_coding_transcript_exon_variant")) +
  geom_line(
    #data=o,
    aes(
      x=Rank,
      y=-log10(uci)
    )
  ) +
  geom_point(
    #data=o_real,
    aes(
      x=Rank,
      y=-log10(p_lrt),
      color = p_lrt < uci,
      shape = effect_simple, 
      alpha = ifelse(effect_simple == "intergenic_region", 0.5, 1)
    )
  ) +
  facet_grid(~chr, scales = "free_x") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) ->
  plor_lat_glm

ggsave(plor_lat_glm, file = "plor_lat_glm.pdf", w = 8, h = 3)

######
#X 15607604
#getData()
SNP_X_15607604=
getData(variant=snps.dt[chr=="X"][pos==15607604]$variant.id,
        samps2use =samps$sampleId)
left_join(SNP_X_15607604, samps) ->
  SNP_X_15607604

SNP_X_15607604 %>%
  filter( (continent %in% c("NorthAmerica","North_America") ) | (continent == "Europe" & lat > 25) ) %>%
  mutate(continent = case_when(continent == "North_America" ~ "NorthAmerica",
                               TRUE ~ continent
                               )) %>%
  group_by(locality,continent) %>%
  summarise(mLat = mean(lat),
            maf = mean(af)
            ) %>%
  ggplot(aes(
    x=mLat,
    y=maf
  )) + geom_point() + geom_smooth(method = "lm") + 
  ggtitle("sog", subtitle = "X:15607604") + 
  facet_wrap(~continent, scales = "free_x") ->
  SNP_X_15607604.plot
ggsave(SNP_X_15607604.plot, file ="SNP_X_15607604.plot.pdf")
#########
#########
####
SNP_2R_19489934=
  getData(variant=snps.dt[chr=="2R"][pos==19489934]$variant.id,
          samps2use =samps$sampleId)
left_join(SNP_2R_19489934, samps) ->
  SNP_2R_19489934

SNP_2R_19489934 %>%
  filter( (continent %in% c("NorthAmerica","North_America") ) | (continent == "Europe" & lat > 25) ) %>%
  mutate(continent = case_when(continent == "North_America" ~ "NorthAmerica",
                               TRUE ~ continent
  )) %>%
  group_by(locality,continent) %>%
  summarise(mLat = mean(lat),
            maf = mean(af)
  ) %>%
  ggplot(aes(
    x=mLat,
    y=maf
  )) + geom_point() + geom_smooth(method = "lm") + 
  ggtitle("TBCB", subtitle = "2R:19489934") + 
  facet_wrap(~continent, scales = "free_x") ->
  SNP_2R_19489934.plot
ggsave(SNP_2R_19489934.plot, file ="SNP_2R_19489934.plot.pdf")

#########
####
SNP_2R_20274917=
  getData(variant=snps.dt[chr=="2R"][pos==20274917]$variant.id,
          samps2use =samps$sampleId)
left_join(SNP_2R_20274917, samps) ->
  SNP_2R_20274917

SNP_2R_20274917 %>%
  filter( (continent %in% c("NorthAmerica","North_America") ) | (continent == "Europe" & lat > 25) ) %>%
  mutate(continent = case_when(continent == "North_America" ~ "NorthAmerica",
                               TRUE ~ continent
  )) %>%
  group_by(locality,continent) %>%
  summarise(mLat = mean(lat),
            maf = mean(af)
  ) %>%
  ggplot(aes(
    x=mLat,
    y=maf
  )) + geom_point() + geom_smooth(method = "lm") + 
  ggtitle("TBCB", subtitle = "2R:19489934") + 
  facet_wrap(~continent, scales = "free_x") ->
  SNP_2R_20003914.plot
ggsave(SNP_2R_20003914.plot, file ="SNP_2R_20003914.plot.pdf")

rbind(
  SNP_2R_19489934,
  SNP_2R_20274917,
  SNP_X_15607604
) ->
  SNPs_Examples

save(SNPs_Examples, file = "SNPs_Examples.X.2R.Rdata")


#/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.Nunez_et_al_Supergene_paper/Phenotyping/GWAS_ldpruned_grm/gwas.tar.gz
#cp /netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.Nunez_et_al_Supergene_paper/Phenotyping/GWAS_ldpruned_grm/gwas.tar.gz ./
# tar -xvzf gwas.tar.gz

### GWAS COLOC
#### Lift over -- chain
LO_c <- "/netfiles/nunezlab/D_melanogaster_resources/liftOver_files/dm6ToDm3.over.chain"
chainObject <- import.chain(LO_c)

gwas_fol <- "/gpfs2/scratch/jcnunez/fst_brent/Coloc_mods/GWAS/lpgwas"
gwas_fils <- system(paste("ls", gwas_fol), intern = T)


colocGWAS =
  foreach(i=gwas_fils, .combine = "rbind",
          .errorhandling = "remove")%do%{
            message(i)
            
            tmp <- fread(paste(gwas_fol, i, sep = "/") )
            tmp %>%
              mutate(p.adj = p.adjust(PVAL, "fdr")) ->
              tmp.adj
            
            grObject <- GRanges(seqnames=c("chr2R"), 
                                ranges=IRanges(start=19115753, 
                                               end=20615753))
            
            LO2r <- as.data.frame(liftOver(grObject, chainObject))
            delta2r = 19115753-LO2r$start
            
            tmp.adj %>% 
              filter(CHR == "2R") %>%
              filter(POS > LO2r$start & POS < LO2r$end) %>%
              mutate(POS_DM6 = POS+delta2r)->
              slice2R
            
            left_join(interest_2R ,slice2R) %>%
              filter(PVAL < 0.05) %>%
              dim() %>% .[1] -> NSPS_2R
            
            ####
            grObjectX <- GRanges(seqnames=c("chrX"), 
                                ranges=IRanges(start=15123410, 
                                               end=16123410))
            LOX <- as.data.frame(liftOver(grObjectX, chainObject))
            deltaX = 15123410-LOX$start
            
            tmp.adj %>% 
              filter(CHR == "X") %>%
              filter(POS > LOX$start & POS < LOX$end) %>%
              mutate(POS_DM6 = POS+deltaX)->
              sliceX
            
            left_join(interest_X ,sliceX) %>%
              filter(PVAL < 0.05) %>%
              dim() %>% .[1] -> NSPS_X
            
            Pheno = i
            
            message(paste(NSPS_2R, NSPS_X, sep = "|"))
            data.table(
              Pheno=Pheno,
              NSPS_2R=NSPS_2R,
              NSPS_X=NSPS_X
            )
            
}


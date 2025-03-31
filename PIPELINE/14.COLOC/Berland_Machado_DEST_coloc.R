library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)
library(forcats)
library(data.table)
library(gmodels)

library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(rnaturalearth)
library(lmtest)

library(VennDiagram)

load("/netfiles02/lockwood_lab/IntrogressionProject/FET_output/Brent_checkpoint_coloc.Rdata")
### load "genes_in_area"          "putative_targets.annot"
### generated in ... PIPELINE/9.2.FET.summaries/Brents_revisiting_analysis.R 
##########
##########
Machado2021_dm6_srt <- get(load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2021.Machado/Machado2021_dm6_srt.Rdata"))
bergland2014_dm6_srt <- get(load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2014.Bergland/bergland2014_dm6_srt.Rdata"))
names(bergland2014_dm6_srt)[3] = "POS"
names(Machado2021_dm6_srt)[3] = "POS"

load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2024.Nunez_et_al_Genetics/DatFor.Haplotypes.trajectory.time.weather.Rdata")
weather.ave %>%
  filter(mod == 2) ->
  weather.2

names(weather.2)[1] = "sampleId_orig"

weather.ave %>%
  filter(mod == 10) ->
  weather.10

names(weather.10)[1] = "sampleId_orig"

#
putative_targets.annot %>%
  left_join(bergland2014_dm6_srt) %>%
  filter(!is.na(clinal.q)) %>%
  ggplot(aes(
    x=POS,
    y=-log10(clinal.q),
    color = effect
  )) + geom_point() +
  geom_hline(yintercept = -log10(0.05))->
  cline_bergland

ggsave(cline_bergland, file = "cline_bergland.pdf")

####
putative_targets.annot %>%
  left_join(bergland2014_dm6_srt) %>%
  filter(!is.na(clinal.q)) %>%
  filter(clinal.q < 0.05) %>% dim

putative_targets.annot %>%
  left_join(Machado2021_dm6_srt) %>%
  filter(!is.na(seas.p)) %>% 
  filter(seas.p < 0.004) %>% dim

putative_targets.annot %>%
  left_join(Machado2021_dm6_srt) %>%
  left_join(bergland2014_dm6_srt) %>%
  select(POS, clinal.q, seas.p,simple_eff, GeneName) ->
  coloc_snps

coloc_snps %>%
  melt(id = c("POS","simple_eff", "GeneName")) %>%
  ggplot(aes(
    x=POS,
    y=-log10(value),
    #color = simple_eff
  )) + geom_point() +
  facet_grid(~variable) +
  xlim(19400000+1000,  20615753+10) + 
  geom_hline(yintercept = -log10(0.004), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05)) +
  theme_bw() ->
  seas_machado_cline_berg

ggsave(seas_machado_cline_berg, 
       file = "seas_machado_cline_berg.pdf",
       w= 6, h = 3)


coloc_snps %>%
  ggplot(aes(
    x=-log10(clinal.q),
    y=-log10(seas.p),
    #color = simple_eff
  )) + geom_point() +
  geom_hline(yintercept = -log10(0.004), 
             linetype = "dashed") +
  geom_vline(xintercept = -log10(0.05))+
  theme_bw() ->
  two_study_ps

ggsave(two_study_ps, 
       file = "two_study_ps.pdf",
       w= 3, h = 3)


coloc_snps %>%
  filter(
    clinal.q < 0.05,
    seas.p < 0.004
  )

# POS    clinal.q      seas.p         simple_eff GeneName
#<num>       <num>       <num>             <char>   <char>
#  1: 20551633 0.009609646 0.001097584 synonymous_variant   BORCS7
####
Machado2021_dm6_srt %>%
  filter(POS == 20551633 & chr == "2R")

#Get_DEST DATA
dest_loc <-"/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/dest.all.PoolSNP.001.50.24Aug2024.ann.gds"
genofile <- seqOpen(dest_loc)
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.csv")

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))

snps.dt %>%
  filter(chr == "2R") %>%
  filter(pos == 20551633) ->
  focal

seqSetFilter(genofile,
             focal$variant.id)

data.table(ann=seqGetData(genofile, "annotation/info/ANN")) ->
  focal_ann

data.table(ann=seqGetData(genofile, "allele")) ->
  alleles


ad <- seqGetData(genofile, "annotation/format/AD")
ad <- ad$data
dp <- seqGetData(genofile, "annotation/format/DP")
print("Create dat object")

dat = ad/dp
####
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , 
                       sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")


dat %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>%
  left_join(samps) ->
  focal_snp_data

######## now bring in the cross
###### ====> CROSS
genofile_cross <- seqOpen("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds")
seqResetFilter(genofile_cross)

snps.CROSS <- data.table(chr=seqGetData(genofile_cross, "chromosome"),
                         pos=seqGetData(genofile_cross, "position"),
                         variant.id=seqGetData(genofile_cross, "variant.id"),
                         alleles=seqGetData(genofile_cross, "allele"),
                         nAlleles=seqNumAllele(genofile_cross))

snps.CROSS %>%
  filter(chr == "2R") %>%
  filter(pos == 20551633) ->
  focal.CROSS

seqSetFilter(genofile_cross, variant.id=focal.CROSS$variant.id, 
             sample.id = c("SKF2_2_1",
                           "SKF3_3_1",
                           "SKF_1_1",
                           "SK_0_1",
                           "VT8F2_2_1",
                           "VT8F3_3_1",
                           "VT8F_1_1",
                           "VT8_0_1"))

data.table(ann=seqGetData(genofile_cross, "allele")) ->
  alleles.cross

message("Allele Freqs")
ad.cross <- seqGetData(genofile_cross, "annotation/format/AD")$data
dp.cross <- seqGetData(genofile_cross, "annotation/format/DP")

dat.cross <- ad.cross/dp.cross
colnames(dat.cross) <- paste(seqGetData(genofile_cross, "chromosome"), seqGetData(genofile_cross, "position"), sep="_")
rownames(dat.cross) <- seqGetData(genofile_cross, "sample.id")

dat.cross %>%
  as.data.frame() %>%
  mutate(sampleId = row.names(.)) %>%
  mutate(x_g = case_when(
    sampleId %in% c("SK_0_1") ~ "SK",
    sampleId %in% c("VT8_0_1") ~ "VT",
    TRUE ~ "cross"
  )) %>%
  group_by(x_g) %>%
  summarise(me.AF.c = mean(`2R_20551633`),
            #med.AF.c = median(`2R_20551633`),
            )-> dat.cross.ag

####
focal_snp_data %>%
  filter(Recommendation == "Pass") %>%
  left_join(weather.10) %>%
  filter(continent == "North_America") %>%
  #filter(sr_season %in% c("spring","fall")) %>%
  #filter(long > -100) %>%
  filter(!is.na(temp.var)) ->
  NAM.weather

#cor.test(~`2R_20551633`+humidity.ave, data = NAM.weather)
#cor.test(~`2R_20551633`+humidity.ave, data = filter(NAM.weather, sr_season == "fall"))
#cor.test(~`2R_20551633`+humidity.ave, data = filter(NAM.weather, sr_season == "spring"))

cor.test(~`2R_20551633`+precip.ave, data = NAM.weather)
cor.test(~`2R_20551633`+precip.ave, data = filter(NAM.weather, sr_season == "fall"))
cor.test(~`2R_20551633`+precip.ave, data = filter(NAM.weather, sr_season == "spring"))

#cor.test(~`2R_20551633`+temp.var, data = NAM.weather)
#cor.test(~`2R_20551633`+temp.var, data = filter(NAM.weather, sr_season == "fall"))
#cor.test(~`2R_20551633`+temp.var, data = filter(NAM.weather, sr_season == "spring"))

#### to REPRODUCE paper ... use only N. America > -90 long samples (i.e., easter seaboard)
foreach(var = vars,
        .combine = "rbind")%do%{
          
          dat_A <- NAM.weather
          AF.a <- dat_F$`2R_20551633`
          wea.a <- dat_F[,which(names(dat.weather) == var)]
          
          dat_F <- NAM.weather %>% filter(sr_season == "fall")
          AF.f <- dat_F$`2R_20551633`
          wea.f <- dat_F[,which(names(dat.weather) == var)]
   
          dat_S <- NAM.weather %>% filter(sr_season == "spring")
          AF.s <- dat_S$`2R_20551633`
          wea.s <- dat_S[,which(names(dat.weather) == var)]
          
          cor.test(AF.f, wea.f) -> tmp.f
          cor.test(AF.s, wea.s) -> tmp.s
          cor.test(AF.a, wea.a) -> tmp.a
          
          data.frame(
            var=var,
            p.fall=tmp.f$p.value,
            s.fall=tmp.s$p.value,
            a.fall=tmp.a$p.value
          )
        }

focal_snp_data %>%
  filter(Recommendation == "Pass") %>%
  left_join(weather.2) %>%
  filter(continent != "Asia") %>%
  #filter(long > -100) %>%
  filter(!is.na(temp.var)) ->
  NAM.weather

NAM.weather %>%
  filter(sr_season %in% c("fall","spring")) ->
  seas.NAM.weather

ggplot() +
  geom_point(
    data = seas.NAM.weather,
      aes(
    x= humidity.ave,
    y= `2R_20551633`,
    color = sr_season,
    shape = as.factor(cluster2.0_k4)
  )) + geom_smooth(
    data = seas.NAM.weather,
    aes(
      x= humidity.ave,
      y= `2R_20551633`,
      color = sr_season,
      linetype =as.factor(cluster2.0_k4)
    ), method = "lm") +
  geom_hline(data = dat.cross.ag, 
             aes(yintercept = me.AF.c, color = x_g)) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
  facet_wrap(~cluster2.0_k4) + theme_bw() ->
  NAM.humidity.ave


ggplot() +
  geom_point(
    data = seas.NAM.weather,
    aes(
      x= temp.var,
      y= `2R_20551633`,
      color = sr_season,
      shape = as.factor(cluster2.0_k4)
    )) + geom_smooth(
      data = seas.NAM.weather,
      aes(
        x= temp.var,
        y= `2R_20551633`,
        color = sr_season,
        linetype =as.factor(cluster2.0_k4)
      ), method = "lm") +
  geom_hline(data = dat.cross.ag, 
             aes(yintercept = me.AF.c, color = x_g)) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
  facet_wrap(~cluster2.0_k4) + theme_bw() ->
  NAM.temp.var

ggplot() +
  geom_point(
    data = seas.NAM.weather,
    aes(
      x= precip.ave,
      y= `2R_20551633`,
      color = sr_season,
      shape = as.factor(cluster2.0_k4)
    )) + geom_smooth(
      data = seas.NAM.weather,
      aes(
        x= precip.ave,
        y= `2R_20551633`,
        color = sr_season,
        linetype =as.factor(cluster2.0_k4)
      ), method = "lm") +
  geom_hline(data = dat.cross.ag, 
             aes(yintercept = me.AF.c, color = x_g)) +
  geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
  facet_wrap(~cluster2.0_k4) + theme_bw() ->
  NAM.precip.ave

ggsave(NAM.humidity.ave, file = "NAM.humidity.ave.pdf",
       w = 9, h = 2.3)
ggsave(NAM.temp.var, file = "NAM.temp.var.pdf",
       w = 9, h = 2.3)
ggsave(NAM.precip.ave, file = "NAM.precip.ave.pdf",
       w = 9, h = 2.3)

###

cor.test(~lat+`2R_20551633`, data =
filter(seas.NAM.weather, continent == "North_America"))

cor.test(~lat+`2R_20551633`, data =
           filter(seas.NAM.weather, continent == "Europe"))
cor.test(~long+`2R_20551633`, data =
           filter(seas.NAM.weather, continent == "Europe"))

#####

#filter(!set %in% c("dgn","dest_plus","DrosEU_3","DrosEU_3_sa") ) %>%

focal_snp_data %>%
  filter(Recommendation == "Pass") %>%
  filter(!is.na(cluster2.0_k4)) %>%
  group_by(city) %>%
  summarise(m.AF = mean(`2R_20551633`, na.rm = T),
            lat=mean(lat),
            long=mean(long),
            cluster2.0_k4 = mean(cluster2.0_k4)
            ) %>%
  filter(!is.na(m.AF) & cluster2.0_k4 %in% 1:4 )-> snp_dat_ag

world <- ne_countries(scale = 110, returnclass = "sf") 

  ggplot() +
    geom_sf(data = world, 
            color = "black", 
            fill = "lightgray", linewidth = 0.1) +
   xlim(-128,160) + ylim(-40, 68) +
    geom_point(
    data = snp_dat_ag,
    size = 1.5,
         aes(
    y=lat,
    x=long,
    fill=m.AF,
    shape=as.factor(cluster2.0_k4)
  )) +
    #geom_vline(xintercept = -100) +
  scale_shape_manual(values = 21:24) +  
  ggtitle("Frequency of C in 2R:20551633") +
  scale_fill_gradient2(midpoint = 0.5)->
  lat_cline

ggsave(lat_cline, file = "lat_cline.pdf")

###Some plots of latitudinallity


####
####
# BEST MODEL
####
####
best_mod <- get(load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2024.Nunez_et_al_Genetics/Revision_Best_Models/bestAIC.v2.Rdata"))
best_mod %>% filter(chr == "2R") %>% 
  filter(cluster == "5.Cville") %>%
  filter(perm_strategy == "new") %>%
  mutate(mod_n = paste(mod,var, sep = "_")) %>%
  ggplot(aes(
    x=fct_reorder(mod_n, rr),
    y=rr,
    color=sig,
    shape = inv
  ), size = 3.5) + geom_hline(yintercept = 1) + 
  geom_point() + theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) ->
  rr_nunez

ggsave(rr_nunez, file = "rr_nunez.pdf",
       w=4, h = 3)

####
####
#### --> seasonal
####
####


focal_snp_data %>%
  filter(Recommendation == "Pass") %>%
  filter(set == "cville")  %>%
  left_join(weather.10) ->
  dat.weather.cville

dat.weather.cville[,c("2R_20551633","precip.ave")] = 
  scale(dat.weather.cville[,c("2R_20551633","precip.ave")])

grangertest(dat.weather.cville$precip.ave,
              dat.weather.cville$`2R_20551633`, 
            order = 1)

cor.test(dat.weather.cville$precip.ave,
         dat.weather.cville$`2R_20551633`)

dat.weather.cville %>%
  group_by(sampleId) %>%
  mutate(Date = strsplit(sampleId, "_")[[1]][5]) ->
  dat.weather.cville

ggplot() +
  geom_line(data =dat.weather.cville,
    aes(
      x= as.Date(Date),
      y= `2R_20551633`,
      group = as.factor(year))
  ) +
    geom_line(data = dat.weather.cville,
      aes(
        x= as.Date(Date),
        y= precip.ave,
        group = as.factor(year)),
      linetype = "dashed"
      ) + facet_wrap(~year, scales = "free_x") ->
  seas_cville

ggsave(seas_cville, file = "seas_cville.pdf",
       w=9, h = 2.5)


###### In this analysis I will ask how many SNPs overall show colocalization.
###### Done on Mar 27, 2025
all_FET <-"/netfiles02/lockwood_lab/IntrogressionProject/FET_output/Fisher_test_results_All_crosses.txt"
all_FET_d <- fread(all_FET)
names(all_FET_d)[1] = "chr"
all_FET_d %<>%
  mutate(SNP_id_dm6 = paste(chr,POS, sep = "_")) %>%
  mutate(pval.adj = p.adjust(sk_all_f.test_pval, 'bonferroni'))
Machado2021_dm6_srt <- get(load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2021.Machado/Machado2021_dm6_srt.Rdata"))
bergland2014_dm6_srt <- get(load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2014.Bergland/bergland2014_dm6_srt.Rdata"))
names(bergland2014_dm6_srt)[3] = "POS"
names(Machado2021_dm6_srt)[3] = "POS"


bergland2014_dm6_srt %>%
  filter(clinal.q < 0.05) %>%
  .$SNP_id_dm6 -> clinal_hits

Machado2021_dm6_srt %>%
  filter(seas.p < 0.004) %>%
  .$SNP_id_dm6 -> seas_hits
sum(clinal_hits %in% seas_hits)

all_FET_d %>%
  filter(pval.adj < 0.01) %>%
  .$SNP_id_dm6 -> intro_hits

x <- list(cline=clinal_hits , 
          seas=seas_hits , 
          intro=intro_hits)


v0 <- venn.diagram( x, filename=NULL, 
                    fill = c("red", "blue", "green"),
                    alpha = 0.50,
                    col = "transparent")

pdf("venn.pdf")
grid.draw(v0)
dev.off()

overlaps <- calculate.overlap(x)
overlaps$a5


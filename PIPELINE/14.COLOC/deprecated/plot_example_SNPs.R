###
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

library(rnaturalearth)

real_dat <-"/netfiles02/lockwood_lab/IntrogressionProject/FET_output/Fisher_test_results_All_crosses.txt"
real_FET <- fread(real_dat)
names(real_FET)[1:2] = c("chr","pos")
real_FET %>%
  filter( (chr == "2R" & pos > 19115753 & pos < 20615753) | (chr == "X" & pos > 15123410 & pos < 16123410) ) %>%
  filter(sk_all_f.test_pval < 5e-08) ->
  boinfe_cutoff

boinfe_cutoff$pos = as.numeric(boinfe_cutoff$pos)

dim(boinfe_cutoff)


#####
load("SNPs_Examples.X.2R.Rdata")
#SNPs_Examples
SNPs_Examples %>%
  filter( (continent %in% c("NorthAmerica","North_America") ) | (continent == "Europe" & lat > 25) ) %>%
  mutate(continent = case_when(continent == "North_America" ~ "NorthAmerica",
                               TRUE ~ continent
  )) %>%
  group_by(chr, pos, continent, locality) %>%
  summarise(
    m.Lat = mean(lat, na.rm = T),
    m.Long = mean(long, na.rm = T),
    m.af = mean(af, na.rm = T)
  ) %>% filter(!is.na(m.af)) ->
  SNPs_Examples.sum


world <- ne_countries(scale = 110, returnclass = "sf") 

ggplot() +
  geom_sf(data = world, 
          color = "black", fill = "lightgray", linewidth = 0.1) +
  theme_classic() +
  geom_point(data = SNPs_Examples.sum,
             aes(
               x=m.Long,
               y=m.Lat,
               fill=m.af
             ), shape = 21, size = 2.5
             ) +
  xlim(-128,36) + ylim(20, 68) +
  facet_grid(paste(chr,pos)~.) +
  scale_fill_gradient2(midpoint = 0.5) +
  ggtitle("Candidate SNPs in DEST 1.0") ->
  DEST_cands_plot

ggsave(DEST_cands_plot, file = "DEST_cands_plot.pdf", w= 6, h = 8)

###### ====> Allele analysis
cand_snps <- SNPs_Examples$variant.id %>% unique

###### ====> DEST
genofile_DEST <- seqOpen("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2021.DEST.v1.release/dest.all.PoolSNP.001.50.10Nov2020.ann.gds")
seqSetFilter(genofile_DEST, variant.id=cand_snps)

snps.dt <- data.table(chr=seqGetData(genofile_DEST, "chromosome"),
                      pos=seqGetData(genofile_DEST, "position"),
                      variant.id=seqGetData(genofile_DEST, "variant.id"),
                      alleles=seqGetData(genofile_DEST, "allele"),
                      nAlleles=seqNumAllele(genofile_DEST))


###### ====> CROSS
genofile_cross <- seqOpen("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds")
seqResetFilter(genofile_cross)

snps.CROSS <- data.table(chr=seqGetData(genofile_cross, "chromosome"),
                      pos=seqGetData(genofile_cross, "position"),
                      variant.id=seqGetData(genofile_cross, "variant.id"),
                      alleles=seqGetData(genofile_cross, "allele"),
                      nAlleles=seqNumAllele(genofile_cross))

left_join(boinfe_cutoff, snps.CROSS) ->
  snps.CROSS.boinfe_cutoff

seqSetFilter(genofile_cross, variant.id=snps.CROSS.boinfe_cutoff$variant.id, 
             sample.id = c("SKF2_2_1",
                           "SKF3_3_1",
                           "SKF_1_1",
                           "SK_0_1",
                           "VT8F2_2_1",
                           "VT8F3_3_1",
                           "VT8F_1_1",
                           "VT8_0_1"
                           ))


message("Allele Freqs")
ad <- seqGetData(genofile_cross, "annotation/format/AD")$data
dp <- seqGetData(genofile_cross, "annotation/format/DP")

dat <- ad/dp
colnames(dat) <- paste(seqGetData(genofile_cross, "chromosome"), seqGetData(genofile_cross, "position") , paste("snp", seqGetData(genofile_cross, "variant.id"), sep = ""), sep="_")
rownames(dat) <- seqGetData(genofile_cross, "sample.id")

dat %>% 
  t() %>%
  as.data.frame() %>%
  mutate(snp_id = row.names(.)) %>%
  separate(snp_id, into = c("chr","pos","feature"), sep = "_") %>%
  dplyr::select(-feature) %>%
  melt(id = c("chr","pos")) ->
  F16_ests

##
SNPs_Examples %>%
filter( (continent %in% c("NorthAmerica","North_America") ) | (continent == "Europe" & lat > 25) ) %>%
  mutate(continent = case_when(continent == "North_America" ~ "NorthAmerica",
                               TRUE ~ continent
  )) ->
  SNPs_Examples_slice

#%>%
#  group_by(locality,continent,chr,pos) %>%
#  summarise(mLat = mean(lat),
#            maf = mean(af, na.rm = T)
#  ) -> SNP_summaries

ggplot() +
  geom_point(
    data=filter(SNPs_Examples_slice, season %in% c("spring","fall")),
    aes(
    x=lat,
    y=af,
    shape=season, color = season
  )) +
  geom_hline(data=F16_ests,
             aes(yintercept = value,
                 color = variable)) +
  facet_grid(paste(chr,pos)~continent, scales = "free_x") ->
  pattern_plots

ggsave(pattern_plots, file = "pattern_plots.pdf")
  


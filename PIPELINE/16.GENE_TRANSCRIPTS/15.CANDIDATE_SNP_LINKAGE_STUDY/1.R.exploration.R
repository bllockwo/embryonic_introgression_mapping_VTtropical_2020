library(vcfR)
library(adegenet)
library(tidyverse)
library(FactoMineR)
library(data.table)
library(foreach)
library(GenomicRanges)
library(rtracklayer)
library(SeqArray)
library(gdsfmt)
library(forcats)
library(SNPRelate)
chain_3to6="/netfiles/nunezlab/D_melanogaster_resources/liftOver_files/dm3ToDm6.over.chain"
chainObject <- import.chain(chain_3to6)
ivs <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/Inversion.status.txt")
names(ivs)[1] = "line"
ivs_2r = ivs %>% dplyr::select(line, `In(2R)NS`)
samps <- fread("/netfiles02/lockwood_lab/IntrogressionProject/Population_files/samp.guide.files.introg.txt")
###


DGRP.vcf <- read.vcfR(
  "./Embryo.snp.dgrp.maf.recode.vcf")

DGRP.gl <- vcfR2genlight(DGRP.vcf) 

tab(DGRP.gl, NA.method = "asis") %>%
  as.data.frame() %>%
  mutate(row.names = row.names(.))  ->
  DGRP.gl.dat
  
  snp_info_dgrp =
  data.frame(
  line =rownames(DGRP.gl.dat),
  snp=DGRP.gl.dat$`2R_16439138_SNP`
  )

write.table(filter(snp_info_dgrp, snp == 2)$line, 
			file = "snp2.txt", 
			append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
write.table(filter(snp_info_dgrp, snp == 0)$line, 
			file = "snp0.txt", 
			append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
#####
##### Investigating LD patterns
####
ld.dgrp <- fread("plink.ld") %>%
  mutate(delta=BP_A-BP_B)

inversion_markers <- fread("/netfiles/nunezlab/D_melanogaster_resources/Inversions/inversion_makers_kapun.dm5.txt")
names(inversion_markers)[c(2,3)] = c("CHR_B","BP_B")
ld.dgrp %>%
  left_join(inversion_markers) -> ld.dgrp

ld.dgrp %>%
  filter(!is.na(inversion))

ggplot() +
  geom_point(
    data=filter(ld.dgrp, R2 > 0.01),
    aes(
      x=abs(delta+1),
      y=R2,
      color = inversion
    ), size = 2.5) +
  geom_smooth(
    data=filter(ld.dgrp),
    aes(
      x=abs(delta+1),
      y=R2,
    ), span = 1/1000) + 
  geom_point(
    data=filter(ld.dgrp,!is.na(inversion)),
    aes(
      x=abs(delta+1),
      y=R2,
      color = inversion
    ), size = 2.5) +
  scale_x_continuous(trans='log10') +
  theme_classic() ->
  LD.DGRP.plot

ggsave(LD.DGRP.plot, file = "LD.DGRP.pdf",
       h= 3.5, w = 4.0)

# Extracting near LD

ld.dgrp %>%
  filter(R2 > .4) %>%
  dplyr::select(chr=CHR_B, dm3.pos=BP_B, delta =delta, R2) ->
  high.ld.partners

#### FST
fst <- fread("Embryo.FST.snp.dgrp.windowed.weir.fst")
fst %>%
  filter(BIN_START > 16439138-1e6 & BIN_END < 16439138+1e6) %>%
  arrange(-WEIGHTED_FST) %>% head(6)

fst %>%
  filter(BIN_START > 16439138-1e6 & BIN_END < 16439138+1e6) %>%
  ggplot(aes(
    x=(BIN_START+BIN_END)/2,
    y=WEIGHTED_FST,
  )) + geom_line(linewidth = 1.5, color = "red") +
  theme_classic() +
  geom_vline(xintercept = 16439138) ->
  fst.DGRP.plot

ggsave(fst.DGRP.plot, file = "fst.DGRP.plot.pdf",
       w= 3.5, h = 3.5)


#### MAKE A PCA
ex_col = which(colnames(DGRP.gl.dat) == "row.names")

  DGRP.gl.dat[-ex_col] %>%
  PCA(graph = F) ->
  pca.info
 
save(pca.info, file = "DGRP.pca.Rdata")
load("DGRP.pca.Rdata")

    pca.info$ind$coord %>%
    as.data.frame %>%
    mutate(line = rownames(.)) %>%
    left_join(snp_info_dgrp) %>%
    left_join(ivs_2r) ->
    DGRP_pca_info
    
  DGRP_pca_info %>%
  ggplot(aes(
  x=Dim.1,
  y=Dim.2,
  shape=`In(2R)NS`,
  color = as.character(snp)
  )) + geom_point() ->
  dgrp_pca_plot
  
  ggsave(dgrp_pca_plot, file = "dgrp_pca_plot.pdf")
  
#####
pi.dgrp0 <- fread("Embryo0.snp.dgrp.windowed.pi") %>% mutate(snp=0)
pi.dgrp2 <- fread("Embryo2.snp.dgrp.windowed.pi") %>% mutate(snp=2)

rbind(pi.dgrp0,pi.dgrp2) %>%
  group_by(snp) %>%
  summarise(m.pi = median(PI))
  
  rbind(pi.dgrp0,pi.dgrp2) %>%
  ggplot(aes(
  x=(BIN_START+BIN_END)/2,
  y=PI,
  color=as.character(snp)
  )) + geom_line() +
  xlim(16439138-5e5,16439138+5e5) +
  geom_vline(xintercept = 16439138) ->
  pi.DGRP.plot
  
  ggsave(pi.DGRP.plot, file = "pi.DGRP.plot.pdf")


####
load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2024.Nunez_et_al_Genetics/DatFor.Haplotypes.trajectory.time.weather.Rdata")
weather.ave %>%
    filter(mod == 2) ->
    weather.2
names(weather.2)[1] = "sampleId_orig"

dest1_loc <-"/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds"
genofile <- seqOpen(dest1_loc)
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/main/populationInfo/dest_v2.samps_8Jun2023.csv")

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))
                      
snps.dt %>%
  filter(chr == "2R") %>%
  filter(pos %in% 20551633) ->
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
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,                      sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")

dat %>%
  as.data.frame() %>%
  mutate(sampleId =rownames(.)) %>%
  left_join(samps)->
  dat_met

dat_met %>%
  group_by(continent) %>%
  summarise(m.AF = median(`2R_20551633`, na.rm = T))

  dat_met %>% filter(set == "cville") %>% 
    left_join(weather.2) %>%
  ggplot(aes(x=temp.var,
             y=`2R_20551633`)
         )+
    geom_point() + geom_smooth(method = "lm") ->
    point_cville
  ggsave(point_cville, file = "point_cville.pdf")
  


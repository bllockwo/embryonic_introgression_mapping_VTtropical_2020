#### PCA Analysis
#### 

#Load libraries
library(tidyverse)
library(vcfR)
library(data.table)
library(adegenet)
library(FactoMineR)
library(magrittr)
library(reshape2)
library(zoo)
library(rtracklayer)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(ggrepel)


# Load the DGRP vcf
DGRP <- read.vcfR(
  "/netfiles/nunezlab/Drosophila_resources/DGRP2/dgrp2.vcf.gz")
DGRPgl <- vcfR2genlight(DGRP) 

###
tab(DGRPgl, NA.method = "asis") %>%
  as.data.frame() %>%
  t() ->
  all_dgrp_data

##### <---- do lifover
rownames(all_dgrp_data) -> SNPs_dgrp.dmel3

data.frame(SNP_id = SNPs_dgrp.dmel3) %>% 
  separate(SNP_id, into = c("chr", "pos", "type")) ->
  GDRP.pos.df

GDRP.pos.df %>%
  filter(chr == "2R") %>%
  filter(type == "SNP") %>%
  mutate(start = pos, end = pos, chr = paste("chr", chr, sep = "")) -> GDRP.pos.df.flt

GDRP.pos.Granges = 
makeGRangesFromDataFrame(GDRP.pos.df.flt,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)

dm3to6 = "/netfiles/nunezlab/Drosophila_resources/liftovers/dm3ToDm6.over.chain"
chain <- import.chain(dm3to6)
liftOver(GDRP.pos.Granges, chain) -> DGRP2.lo.pos
unlist(DGRP2.lo.pos) %>% as.data.frame() -> DGRP2.lo.pos.df

which(lengths(DGRP2.lo.pos) == 0) -> missing_snps

cbind(
as.data.frame(GDRP.pos.Granges)[-missing_snps,],
DGRP2.lo.pos.df) -> lift_over_file

lift_over_file[,c(1,2,7)] -> lift_over_file.sp
names(lift_over_file.sp) = c("chr", "dm3.pos", "dm6.pos")

lift_over_file.sp %<>%
  mutate(dm3.id = paste(gsub("chr", "", chr),dm3.pos,"SNP", sep = "_"),
         dm6.id = paste(gsub("chr", "", chr),dm6.pos,"SNP", sep = "_")
         )
####

all_dgrp_data[lift_over_file.sp$dm3.id,] -> all_dgrp_data.dm3to6

rownames(all_dgrp_data.dm3to6) = lift_over_file.sp$dm6.id

#Separate data based on inversion set
#Metadata
invst <- fread("/netfiles/nunezlab/Drosophila_resources/DGRP2/Inversion.status.txt")

standard  <- invst$DGRP_Line[which(invst$`In(2R)NS` == "ST")]
inversion <- invst$DGRP_Line[which(invst$`In(2R)NS` == "INV")]
heterokaryotype <- invst$DGRP_Line[which(invst$`In(2R)NS` == "INV/ST")]

#####
#Make standard set
all_dgrp_data.dm3to6 %>%
  .[,which(colnames(.) %in% standard)] ->
  standard_lines

all_dgrp_data.dm3to6 %>%
  .[,which(colnames(.) %in% inversion)] ->
  inverted_lines

all_dgrp_data.dm3to6 %>%
  .[,which(colnames(.) %in% heterokaryotype)] ->
  heterokaryotype_lines

######## Define functions to stud

filter_drgp_sets = function(x, miss_thresh){
  
  count_NA = function(x) {
    x  %>%   
      .[is.na(.)] %>%
      length() ->
      numNA
    return(numNA/length(x))
  }
  
  apply(x,
        2, 
        FUN=count_NA
  ) -> x_NA_lines
  
  x_keep_samps <- names(x_NA_lines[x_NA_lines <= miss_thresh])
  
  x %>%
    .[,which(colnames(.) %in%  x_keep_samps)] -> 
    x_lines_samps_filt
  
  #count NA per SNPs
  apply(x_lines_samps_filt,
        1, 
        FUN=count_NA
  ) -> x_NA_SNPs
  
  x_keep_snps <- names(x_NA_SNPs[x_NA_SNPs <= miss_thresh])
  
  x_lines_samps_filt %>%
    .[which(rownames(.) %in%  x_keep_snps),] -> 
    x_lines_samps_SNPs_filt
  
  x_lines_samps_SNPs_filt_NAimp = na.aggregate(x_lines_samps_SNPs_filt,
                                               FUN= median)
  
  return(x_lines_samps_SNPs_filt_NAimp)
  
}

standard_lines_flt = filter_drgp_sets(standard_lines, 1)

inverted_lines_flt = filter_drgp_sets(inverted_lines, 1)

heterokaryotype_lines_flt = filter_drgp_sets(heterokaryotype_lines, 1)

standard_lines_flt %>% dim
inverted_lines_flt %>% dim
heterokaryotype_lines_flt %>% dim

###
cbind(standard_lines_flt,inverted_lines_flt) -> DGRP.2R.INVSTD
t(DGRP.2R.INVSTD) -> DGRP.2R.INVSTD.mtrx


####

### Load the genotype locality
genofile.path <- "/netfiles/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds"

### load the genofile into memory
genofile <- seqOpen(genofile.path)

###
#### Load the filtering SNP object -- which JCBN created
filtering.dt <- get(load("/netfiles/lockwood_lab/IntrogressionProject/SNPcalling_output/SNP.filtering.guide.Rdata"))

filtering.dt %>%
  filter(is.na(libs)) %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) ->
  chosen.snps

###

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))


snps.dt %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) %>%
  right_join(chosen.snps, by = c("chr", "pos", "SNP_id") ) ->
  snps.dt

###
###
###

snps.dt <- snps.dt[nAlleles==2]

seqSetFilter(genofile, variant.id=snps.dt$variant.id)

snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

seqSetFilter(genofile,
             snps.dt[chr%in%c("2L", "2R", "3L", "3R", "X")][missing<.05]$variant.id)


#### extract matrices of cov and calls
### get allele frequency data
ad <- seqGetData(genofile, "annotation/format/AD")$data
dp <- seqGetData(genofile, "annotation/format/DP")

### divide the matrices to generate allele frequency calls
dat <- ad/dp

### check the dimensions of the allele frequency matrix
dim(dat)  

colnames(dat) <- paste(seqGetData(genofile, "chromosome"), 
                       seqGetData(genofile, "position"), 
                       "SNP",  sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")

dat[,which(colnames(dat) %in% lift_over_file.sp$dm6.id)] -> dat.2R

####
rbind(DGRP.2R.INVSTD.mtrx[,colnames(dat.2R)],
      dat.2R
      ) -> merged.dat

####
merged.dat %>% 
  PCA(graph = F, scale.unit = F) ->
  PCA.DGRP.EMBR.obj

invst.samp = invst
names(invst.samp)[1] = "sampleId"

PCA.DGRP.EMBR.obj$ind$coord %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>%
  left_join(invst.samp) ->
  pca.merged.samps
  
setDT(pca.merged.samps)

  ggplot() +
  geom_point(
    data = pca.merged.samps[grep("line", pca.merged.samps$sampleId )],
    aes(
    x=Dim.1,
    y=Dim.2,
    color = `In(2R)NS`
  )) +
  geom_point(
    data = pca.merged.samps[grep("0_1", pca.merged.samps$sampleId )],
    aes(
      x=Dim.1,
      y=Dim.2,
      color = `In(2R)NS`
    )) +
  geom_text_repel(data = pca.merged.samps[grep("0_1", pca.merged.samps$sampleId )],
                  aes(x=Dim.1,
                      y=Dim.2,
                      label = sampleId),
                  nudge_x = 0.08) ->
  DGRP.IN.plot

ggsave(DGRP.IN.plot, file = "DGRP.IN.plot.pdf")

###

ggplot() +
  geom_point(
    data = pca.merged.samps[grep("line", pca.merged.samps$sampleId )],
    aes(
      x=Dim.2,
      y=Dim.3,
      color = `In(2R)NS`
    )) +
  geom_point(
    data = pca.merged.samps[grep("0_1", pca.merged.samps$sampleId )],
    aes(
      x=Dim.2,
      y=Dim.3,
      color = `In(2R)NS`
    )) +
  geom_text_repel(data = pca.merged.samps[grep("0_1", pca.merged.samps$sampleId )],
                  aes(x=Dim.2,
                      y=Dim.3,
                      label = sampleId),
                  nudge_x = 0.08) ->
  DGRP.IN.plot3

ggsave(DGRP.IN.plot3, file = "DGRP.IN.plot3.pdf")

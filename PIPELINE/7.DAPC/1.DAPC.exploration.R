
library(SeqArray)
library(tidyverse)
library(magrittr)
library(reshape2)
library(vroom)
library(data.table)
library(matrixStats)
library(LEA)
library(gdsfmt)
library(SNPRelate)
library(rnaturalearth)
library(rnaturalearthdata)
library(FactoMineR)
library(factoextra)
library(patchwork)
library(adegenet)
library(ggrepel)
library(foreach)
library(doParallel)
library(doMC)
registerDoMC(4)


genofile.path <- "/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds"

genofile <- seqOpen(genofile.path)


snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))


snps.dt <- snps.dt[nAlleles==2]

seqSetFilter(genofile, variant.id=snps.dt$variant.id)

snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]

snps.dt %>%
  group_by(chr) %>%
  summarise(N = n())

snps.dt %>%
  filter(missing<.05) %>%
  filter(af > 0.05) %>%
  #group_by(chr) %>%
  summarise(N = n())

snps.dt %>%
  filter(missing<.05 & af > 0.5) %>%
  filter(chr%in%c("2L", "2R", "3L", "3R")) %>%
  summarise(N = n())

snps.dt %>%
  filter(missing<.05 & af > 0.5) %>%
  filter(chr%in%c("X")) %>%
  summarise(N = n())


seqSetFilter(genofile,
             snps.dt[chr%in%c("2L", "2R", "3L", "3R")][missing<.05]$variant.id)

### get allele frequency data
ad <- seqGetData(genofile, "annotation/format/AD")
dp <- seqGetData(genofile, "annotation/format/DP")

dat <- ad$data/dp
dim(dat)  

colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , paste("snp", seqGetData(genofile, "variant.id"), sep = ""), sep="_")

rownames(dat) <- seqGetData(genofile, "sample.id")
####
dat %>%
  t() %>% 
  as.data.frame -> dat_AF_samps_target

## Some characterizations of AFs and subsequent filtering
MeanAF=c()
MinAF=c()
MaxAF=c()


apply(dat_AF_samps_target,
      1, FUN=mean, na.rm=TRUE ) -> MeanAF
apply(dat_AF_samps_target,
      1, FUN=min, na.rm=TRUE ) -> MinAF
apply(dat_AF_samps_target,
      1, FUN=max, na.rm=TRUE ) -> MaxAF


cbind(dat_AF_samps_target, MeanAF, MinAF, MaxAF) -> dat_AF_samps_target


### ---> only poly
dat_AF_samps_target %>%
  .[which(.$MeanAF > 0.00 & .$MeanAF < 1.00),] %>%
  .[which(.$MinAF > 0.001 & MaxAF < 0.999),] ->  ### This samples only polymorphic sites
  dat_AF_samps_target_poly

dat_AF_samps_target_poly %>% dim

dat_AF_samps_target_poly %>%
  dplyr::select(!c(MeanAF,MinAF,MeanAF,MinAF,MaxAF)) %>% 
  t %>%
  as.data.frame() ->
  dat.poly

#### --> only fixed
dat_AF_samps_target %>%
  filter(CH_0_1 %in% c(0,1)) %>%
  filter(SK_0_1 %in% c(0,1)) %>%
  filter(VT10_0_1 %in% c(0,1)) %>%
  filter(VT10_0_1 %in% c(0,1)) ->
  dat_AF_samps_target_fix

dat_AF_samps_target_fix %>%
  dplyr::select(!c(MeanAF,MinAF,MeanAF,MinAF,MaxAF)) %>% 
  t %>%
  as.data.frame() ->
  dat.fix

####
data.frame(sample.id = rownames(dat),
           type = c("F16x", "F16x", "F16x", "Par", "Wild", "Wild",
                    "F16x", "F16x", "F16x", "Par",
                    "F16x", "F16x", "F16x", "Par",
                    "F16x", "F16x", "F16x", "Par"
                    ),
           habi = c("exp", "exp", "exp", "Trop", "Temp", "Trop",
                    "exp", "exp", "exp", "Trop",
                    "exp", "exp", "exp", "Temp",
                    "exp", "exp", "exp", "Temp"
           ),
           orig = c("CH", "CH", "CH", "CH", "ME", "PAN",
                    "SK", "SK", "SK", "SK",
                    "VT", "VT", "VT", "VT",
                    "VT", "VT", "VT", "VT"
           )
           ) -> metadat


### all
dapc(dat, grp = metadat$habi, n.pca = 6, n.da = 3) -> habi.dapc
alpha.opt <- optim.a.score(habi.dapc) 

#pdf("catter.pdf")
#scatter.dapc(habi.dapc)
#dev.off()

habi.dapc$ind.coord %>%
  as.data.frame() %>%
  mutate(sample.id = rownames(.)) %>%
  left_join(metadat) ->
  dapc.data.df

dapc.data.df %>%
  ggplot(aes(
    x= LD1,
    y= LD2,
    fill=habi,
    shape = orig,
    label = orig
  )) +
  geom_jitter(size = 3, aes()) +
  theme_bw() + 
  geom_text_repel(data = subset(dapc.data.df, habi != "exp"),
                  size = 5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  scale_shape_manual(values = 21:26) ->
  dapc.samps

ggsave(dapc.samps, file = "dapc.samps.pdf", w = 5, h =4)

### poly
### ### poly
### poly
### poly

dapc(dat.poly, grp = metadat$habi, n.pca = 6, n.da = 3) -> habi.dapc.poly

#pdf("catter.pdf")
#scatter.dapc(habi.dapc)
#dev.off()

habi.dapc.poly$ind.coord %>%
  as.data.frame() %>%
  mutate(sample.id = rownames(.)) %>%
  left_join(metadat) ->
  dapc.data.df.poly

dapc.data.df.poly %>%
  ggplot(aes(
    x= LD1,
    y= LD2,
    fill=habi,
    shape = orig,
    label = orig
  )) +
  geom_jitter(size = 3, aes()) +
  theme_bw() + 
  geom_text_repel(data = subset(dapc.data.df.poly, habi != "exp"),
                  size = 5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  scale_shape_manual(values = 21:26) ->
  dapc.samps.poly

ggsave(dapc.samps.poly, file = "dapc.samps.poly.pdf", w = 5, h =4)


### fixed
### ### fixed
### fixed
### fixed

dapc(dat.fix, grp = metadat$habi, n.pca = 6, n.da = 3) -> habi.dapc.fix

#pdf("catter.pdf")
#scatter.dapc(habi.dapc)
#dev.off()

habi.dapc.fix$ind.coord %>%
  as.data.frame() %>%
  mutate(sample.id = rownames(.)) %>%
  left_join(metadat) ->
  dapc.data.df.fix

dapc.data.df.fix %>%
  ggplot(aes(
    x= LD1,
    y= LD2,
    fill=habi,
    shape = orig,
    label = orig
  )) +
  geom_jitter(size = 3, aes()) +
  theme_bw() + 
  geom_text_repel(data = subset(dapc.data.df.fix, habi != "exp"),
                  size = 5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines")) +
  scale_shape_manual(values = 21:26) ->
  dapc.samps.fix

ggsave(dapc.samps.fix, file = "dapc.samps.fix.pdf", w = 5, h =4)

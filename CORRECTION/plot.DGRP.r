### PlOT DGRP windows

## reproducible code
library(plyr); library(dplyr)
library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(liftOver)

#setwd(".../GitHub/embryonic_introgression_mapping_VTtropical_2020/CORRECTION/")

dat <- fread("DGRP.ProperlyRelabeled.data.txt")
## Generate a joint genotype label for Sp70 and sog; called "geno"
dat %<>%
  mutate(geno = paste("SP70=",SP70_allele,"-","sog=",SOG_allele, 
                      sep =""))


##### Add DGRP genotypes
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(adegenet)
library(FactoMineR)
library(gtools)
library(foreach)
library(ape)

DGRPf <- "/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/dgrp2.gds.gz"
DGRP <- seqOpen(DGRPf, allow.duplicate=T)

####
seqResetFilter(DGRP)
### RESRET
### RESRET
### 
snps.dt <- data.table(chr=seqGetData(DGRP, "chromosome"),
                      pos=seqGetData(DGRP, "position"),
                      variant.id=seqGetData(DGRP, "variant.id"),
                      alleles=seqGetData(DGRP, "allele"),
                      nAlleles=seqNumAllele(DGRP),
                      missing=seqMissing(DGRP))

#chr      beg      end    size nsnps Lindley score peak pos.
#1  2R 19359398 21441591 2082193  5567                20475826
#Lindley score peak value -log10(p-val) peak pos -log10(p-val) peak value
#1                 1949.076               19489892         14.5510783018925

beg=19359398-4112495
end=21441591-4112495

snps.dt %>%
  filter(chr == "2R") %>%
  filter(nAlleles == 2) %>%
  filter(pos > beg & pos < end) -> flt_2R

grObject <- GRanges(seqnames=c("chr2R"), 
                    ranges=IRanges(start=min(flt_2R$pos), 
                                   end=max(flt_2R$pos)))
chainObject <- import.chain("/netfiles/nunezlab/D_melanogaster_resources/liftOver_files/dm3ToDm6.over.chain")
# run liftOver
flt_2R.LO <- as.data.frame(
  liftOver(grObject, chainObject)) %>%
  mutate(    O.min=min(flt_2R$pos),
             O.max=max(flt_2R$pos))
flt_2R.LO

#chr      beg      end   size nsnps Lindley score peak pos.
#2   X 15432794 16309492 876698  1126                16009527
#Lindley score peak value -log10(p-val) peak pos -log10(p-val) peak value
#2                 601.2125               15843297         13.0762428629402

begx=15432794-105967
endx=16309492-105967

snps.dt %>%
  filter(chr == "X") %>%
  filter(nAlleles == 2) %>%
  filter(pos >= begx & pos <= endx) -> flt_X

grObject <- GRanges(seqnames=c("chrX"), 
                    ranges=IRanges(start=min(flt_X$pos), 
                                   end=max(flt_X$pos)))
chainObject <- import.chain("/netfiles/nunezlab/D_melanogaster_resources/liftOver_files/dm3ToDm6.over.chain")
# run liftOver
flt_X.LO <- as.data.frame(
  liftOver(grObject, chainObject)) %>%
  mutate(    O.min=min(flt_X$pos),
             O.max=max(flt_X$pos))
flt_X.LO

#flt2L
rbind(flt_2R, flt_X) -> filt.DGRP

filt.DGRP %<>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) %>%
  filter(missing < 0.05) %>%
  separate(alleles, remove = F, into = c("REF","ATL"), sep = ",")

filt.DGRP %>%
  filter(REF %in% c("A","C","T","G")) %>%
  filter(ATL %in% c("A","C","T","G")) ->
  filt.DGRP.flt

### Geno Step
## outliers -> "line_69", "line_304"

samps = str_replace(dat$line, "DGRP", "line")

snpgdsGetGeno(DGRP, 
              sample.id = samps[-which(samps %in% 
                                         c("line_69", 
                                           "line_304",
                                           "line_852"))],
              verbose=TRUE, 
              with.id=TRUE,
              snp.id =filt.DGRP.flt$variant.id) ->
  GENO.matrix

MATRIX <- GENO.matrix$genotype
MATRIX = as.data.frame(MATRIX)
rownames(MATRIX)  <- GENO.matrix$sample.id
names(MATRIX)  <- filt.DGRP.flt$SNP_id

#####
### PCA
MATRIX%>% 
  PCA(graph = F #, ind.sup= 206:207
  ) ->
  pca.mtrix

pca.mtrix$ind$coord %>%
  as.data.frame() %>%
  mutate(line=rownames(.)) %>%
  mutate(line = str_replace(line, "line", "DGRP")) %>%
  mutate(line = str_replace(line, "line", "DGRP")) %>%
  full_join(dat) ->
  pca.ind.dat

pca.ind.dat %<>%
  filter(!is.na(geno))

cor.test(pca.ind.dat$Dim.1, pca.ind.dat$Proportion ) ->a1
cor.test(pca.ind.dat$Dim.2, pca.ind.dat$Proportion ) ->a2
cor.test(pca.ind.dat$Dim.3, pca.ind.dat$Proportion ) ->a3
cor.test(pca.ind.dat$Dim.4, pca.ind.dat$Proportion ) ->a4

data.frame(
  #L1 = i,
  #L2 = j,
  D1p = a1$p.value,
  D1r = a1$estimate,
  D2p = a2$p.value,
  D2r = a2$estimate)

#}}

#save(cor_emb, file = "cor_emb.Rdata")
#load("cor_emb.Rdata")
pca.ind.dat %>%
  ggplot(aes(
    x=Dim.1, y = Dim.2,
    color=Proportion, shape = paste(geno)
  )) + geom_point(size = 3) + 
  scale_color_gradient2(midpoint = 0.3, high = "red", low = "blue") ->
  plot_pca
ggsave(plot_pca, file = "plot_pca.pdf")

###
pca.ind.dat %>%
  filter(!is.na(Dim.2)) %>%
  ggplot(aes(
    x=Dim.2,
    y=Proportion,
  )) + geom_point() + geom_smooth(method = "lm")  ->
  PCA.Heat.corr
ggsave(PCA.Heat.corr, 
       w=7,h=4,
       file = "Dim2.Prop.pdf")



pca.ind.dat %>%
  filter(!is.na(Dim.2)) %>%
  ggplot(aes(
    x=paste(geno),
    y=Proportion,
  )) + geom_boxplot() +
  facet_wrap(~Dim.2>0) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) ->
  PCA.box.emb

ggsave(PCA.box.emb, 
       w=7,h=4,
       file = "PCA.box.emb.pdf")


###
library(ggtree)
library(phytools)

dgrp_tree <- nj(dist(MATRIX))
pdf("tree.pdf")
plot(dgrp_tree)
dev.off()

p <- ggtree(dgrp_tree, branch.length="none") + geom_tiplab()

#p <- ggtree(dgrp_tree, branch.length="none") %<+% info


ggsave(p, file = "tree.pdf")



p <- ggtree(dgrp_tree, ignore.negative.edge=TRUE) %<+% 
  pca.ind.dat + xlim(-.1, 4)

p2 <- p + geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint(aes(shape = paste(SP70_allele, SOG_allele), color = Proportion, 
                    size = Proportion)) + 
  theme(legend.position = "right") + 
  scale_size_continuous(range = c(3, 10))

ggsave(p, file = "tree.pdf")



### ADD POOLS
POOLS <- "/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds"
genofile <- seqOpen(POOLS, allow.duplicate=T)
seqResetFilter(genofile)
seqGetData(genofile, "sample.id")

snps.dt.pool <- data.table(chr=seqGetData(genofile, "chromosome"),
                           pos=seqGetData(genofile, "position"),
                           variant.id=seqGetData(genofile, "variant.id"),
                           alleles=seqGetData(genofile, "allele"),
                           nAlleles=seqNumAllele(genofile),
                           missing=seqMissing(genofile))

snps.dt.pool %<>% mutate(SNP_id = paste(chr,pos, sep = "_"))

beg.P=19359398
end.P=21441591

begx.P=15432794
endx.P=16309492

snps.dt.pool %>%
  filter(chr == "2R") %>%
  filter(nAlleles == 2) %>%
  filter(pos >= 19359398 & pos < 21441591) -> flt_2R.pool

snps.dt.pool %>%
  filter(chr == "X") %>%
  filter(nAlleles == 2) %>%
  filter(pos >= 15432794 & pos <= 16309492) -> flt_X.pool

### Modify
flt_2R.pool %<>%
  mutate(cor.pos = pos-4112495)
flt_X.pool %<>%
  mutate(cor.pos = pos-105967)
rbind(flt_2R.pool, flt_X.pool) -> filt.exp.pools
filt.exp.pools %<>%
  mutate(liftover_id = paste(chr, cor.pos, sep = "_"))

select_join = filt.exp.pools$variant.id[which(filt.exp.pools$liftover_id 
                                              %in% 
                                                filt.pools$SNP_id)]
snpgdsGetGeno(genofile, 
              verbose=TRUE, with.id=TRUE,
              snp.id =select_join,
              sample.id= c("SK_0_1","VT8_0_1") ) ->
  POOLS.matrix

POOL.MATRIX <- POOLS.matrix$genotype
POOL.MATRIX = as.data.frame(POOL.MATRIX)
rownames(POOL.MATRIX)  <- c("SK","VT")
names(POOL.MATRIX)  <- filt.exp.pools$liftover_id[which(filt.exp.pools$liftover_id 
                                                        %in% 
                                                          filt.pools$SNP_id)]

rbind.fill(MATRIX, POOL.MATRIX) ->MP
rownames(MP) = c(rownames(MATRIX) , rownames(POOL.MATRIX))

#####
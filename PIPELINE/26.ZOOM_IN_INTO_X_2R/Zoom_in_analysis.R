#### Explore X-associated SNP
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
library(fastglm)
#library(rnaturalearth)

#### Annotations
#####
##create get data function
getAnnot <- function(chr.i="2L", pos.i=14617051) {
  # chr.i="2R"; pos.i=20551633;
  ### filter to target
  snp.tmp <- snp.dt %>%
    filter(chr == chr.i) %>%
    filter(pos == pos.i)
  seqSetFilter(genofile, variant.id=snp.tmp$variant.id)
  ### get annotations
  #message("Annotations")
  tmp <- seqGetData(genofile, "annotation/info/ANN")
  data.frame(
    chr=chr.i,
    pos=pos.i,
    annot = tmp$data) %>%
    separate(annot,
             into = c(
               "allele",
               "type",
               "effect",
               "gene_name",
               "fbgn",
               "mRNA",
               "fbtr",
               "feature",
               "shift",
               "SNP_notation",
               "AA_notation",
               "Alt_pos1",
               "Alt_pos2",
               "Alt_pos3",
               "dist_to_feature",
               "ect1"
             ), sep = "\\|") -> o
  
  return(o)
}
genofile <- seqOpen("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds")
seqResetFilter(genofile)

snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     variant.id=seqGetData(genofile, "variant.id"),
                     alleles=seqGetData(genofile, "allele"),
                     nAlleles=seqNumAllele(genofile))

####
#### Get Top SNPs
all_FET <-"/netfiles02/lockwood_lab/IntrogressionProject/FET_output/Fisher_test_results_All_crosses.txt"
all_FET_d <- fread(all_FET)
names(all_FET_d)[1] = "chr"
all_FET_d %<>%
  mutate(SNP_id_dm6 = paste(chr,POS, sep = "_")) %>%
  mutate(pval.adj_all = p.adjust(sk_all_f.test_pval, 'bonferroni'),
         pval.adj_skf1 = p.adjust(skf1_f.test_pval, 'bonferroni'),
         pval.adj_skf2 = p.adjust(skf2_f.test_pval, 'bonferroni'), 
         pval.adj_skf3 = p.adjust(skf3_f.test_pval, 'bonferroni'),
         pval.adj_vt8f1 = p.adjust(vt8f1_f.test_pval, 'bonferroni'), 
         pval.adj_vt8f2 = p.adjust(vt8f2_f.test_pval, 'bonferroni'), 
         pval.adj_vt8f3 = p.adjust(vt8f3_f.test_pval, 'bonferroni'))

all_FET_d %>%
  filter(pval.adj_all  < 0.01) ->
  All_FET_mutations

all_FET_d %>%
  filter(pval.adj_all  < 0.01) %>%
  filter(chr == "2R") %>%
  filter(POS >= 19115753 & POS <= 20615753) ->
  top_muts_2R

all_FET_d %>%
  filter(pval.adj_all  < 0.01) %>%
  filter(chr == "X") %>%
  filter(POS >= 15123410 & POS <= 16123410) ->
  top_muts_X

rbind(top_muts_2R, top_muts_X) -> top_muts_all

annots_top_fet = foreach(i=1:dim(top_muts_all)[1], 
                         .combine = "rbind")%do%{
                           
                           tmp = top_muts_all[i,]
                           
                           getAnnot(tmp$chr, tmp$POS)
                         }

annots_top_fet %>%
  filter(gene_name %in% c("sog","CG13430")) %>%
  group_by(chr) %>% slice_min(pos)
annots_top_fet %>%
  filter(gene_name %in% c("sog","CG13430")) %>%
  group_by(chr) %>% slice_max(pos)

annots_top_fet %>%
  filter(pos %in% c(15648874))

### for X = 15602941 -> 15649874
### for 2R = 20550633 -> 20556321

all_FET_d %>%
  filter(chr == "2R") %>%
  filter(POS > 20550633 & POS < 20556321) %>%
  mutate(V10 = "CG13430", V4 = POS)->
  FET2R

all_FET_d %>%
  filter(chr == "X") %>%
  filter(POS > 15601941 & POS < 15649874)%>%
  mutate(V10 = "CG9224", V4 = POS) ->
  FETX


####
### Gene Structures
structs <- fread("Target_genes.txt")
structs %>%
  filter(V3 == "exon") ->
  exons

exons$splice[grep("variant A", exons$V9)] = 1
exons$splice[grep("variant B", exons$V9)] = 2
exons$splice[grep("variant E", exons$V9)] = 3

#### 2R:20551633 and X:15602941
top_hits = data.frame(
  V10=c("CG13430","CG9224"),
  V4=c(20551633, 15602941)
)


ggplot() +
  geom_vline(data = top_hits,
             aes(xintercept = V4/1e6), color = "red") +
  geom_jitter(
    data = data.frame(rbind(FET2R, FETX)),
    aes(
      x=V4/1e6,
      y=-log10(pval.adj_all)
    )) +
geom_line(
    data = data.frame(rbind(FET2R, FETX)),
    aes(
      x=V4/1e6,
      y=-log10(pval.adj_all)
     )) +
  geom_segment(
    data=exons,
    linewidth = 1.5,
    aes(
      x=V4/1e6,
      xend=V5/1e6,
      y=-as.numeric(splice),
      yend=-as.numeric(splice)
    )) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = -log10(0.01), 
             linetype = "dashed") +
  facet_grid(~V10,scales = "free") +
  theme_classic()->
  gene_struct_plot

ggsave(gene_struct_plot, 
       w = 9, h = 2.5,
       file = "gene_struct_plot.pdf")

#####
annots_top_X = foreach(i=1:dim(filter(FETX, pval.adj_all < 0.01 ))[1], 
                         .combine = "rbind")%do%{
                           
                           tmp = filter(FETX, pval.adj_all < 0.01 )[i,]
                           
                           getAnnot(tmp$chr, tmp$POS)
                         }

unique(annots_top_X$pos)

annots_top_2r = foreach(i=1:dim(filter(FET2R, pval.adj_all < 0.01 ))[1], 
                       .combine = "rbind")%do%{
                         
                         tmp = filter(FET2R, pval.adj_all < 0.01 )[i,]
                         
                         getAnnot(tmp$chr, tmp$POS)
                       }


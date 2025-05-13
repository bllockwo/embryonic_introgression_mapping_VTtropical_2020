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

#####
genofile <- seqOpen("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds")
seqResetFilter(genofile)

snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                         pos=seqGetData(genofile, "position"),
                         variant.id=seqGetData(genofile, "variant.id"),
                         alleles=seqGetData(genofile, "allele"),
                         nAlleles=seqNumAllele(genofile))


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
               "shit",
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

annots_top_fet = foreach(i=1:dim(top_muts_all)[1], 
                         .combine = "rbind")%do%{
                           
                           tmp = top_muts_all[i,]
                           
                           getAnnot(tmp$chr, tmp$POS)
                         }

annots_top_fet %>%
  group_by(chr, pos, gene_name) %>%
  slice_head() ->
  sliced_genes_FET_top

write.table(sliced_genes_FET_top, 
            file = "sliced_genes_FET_top.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

annots_top_fet %>%
  ungroup() %>%
  group_by(chr, gene_name) %>%
  slice_head() %>%
  select(chr, gene_name, fbgn)->
  sliced_genes_FET_top_genes_only

write.table(sliced_genes_FET_top_genes_only, 
            file = "sliced_genes_FET_top_genes_only.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")



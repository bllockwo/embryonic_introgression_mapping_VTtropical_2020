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
##create get data function
#getAnnot("X", 15743371)
getAnnot("X", 15602941)
getAnnot("X", 15607604)
getAnnot("X", 15847814)
#getAnnot("2R", 20551633)
#chr      pos allele                    type   effect  gene_name
#1   X 15602941      G downstream_gene_variant MODIFIER        sog
#2   X 15602941      G downstream_gene_variant MODIFIER        sog
#3   X 15602941      G downstream_gene_variant MODIFIER        sog
#4   X 15602941      G       intergenic_region MODIFIER CG9220-sog
#1   X 15607604      T synonymous_variant    LOW       sog FBgn0003463
#2   X 15607604      T synonymous_variant    LOW       sog FBgn0003463
#3   X 15607604      T synonymous_variant    LOW       sog FBgn0003463
#1   X 15847814      C     3_prime_UTR_variant MODIFIER    CG8578 FBgn0030699
#2   X 15847814      C   upstream_gene_variant MODIFIER     ArgRS FBgn0027093
#3   X 15847814      C downstream_gene_variant MODIFIER    CG8578 FBgn0030699
#4   X 15847814      C downstream_gene_variant MODIFIER       sun FBgn0014391
#5   X 15847814      C downstream_gene_variant MODIFIER      UBL3 FBgn0026076
#6   X 15847814      C downstream_gene_variant MODIFIER    CG8578 FBgn0030699
#7   X 15847814      C downstream_gene_variant MODIFIER       sun FBgn0014391
#8   X 15847814      C downstream_gene_variant MODIFIER      UBL3 FBgn0026076
#9   X 15847814      C downstream_gene_variant MODIFIER      UBL3 FBgn0026076

#### FET Annotations

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

#### haplotypic mutations
haplo_SNPs = c("2R_20550549_SNP",
               "2R_20551182_SNP",
               "2R_20551519_SNP",
               "2R_20550323_SNP",
               "2R_20552286_SNP",
               "2R_20551486_SNP",
               "2R_20551129_SNP",
               "2R_20550070_SNP",
               "2R_20550046_SNP",
               "2R_20550403_SNP")

haplo_snos_o = 
data.frame(SNP_id = haplo_SNPs) %>%
  separate(SNP_id, remove = F, 
           into = c("chr", "POS", "feature"),
           sep = "_")

annots_haplo = foreach(i=1:dim(haplo_snos_o)[1], 
                         .combine = "rbind",
                       .errorhandling = "remove")%do%{
                           
                           tmp = haplo_snos_o[i,]
                           
                           getAnnot(tmp$chr, tmp$POS)
                         }

annots_haplo %>% 
  filter(gene_name == "CG13430") %>%
  group_by(pos) %>%
  slice_head() %>%
  mutate(snp_id = paste(paste("X",chr, sep = ""), pos, "SNP", sep = "_")) ->
  annotaions_snps


#####
load("/gpfs2/scratch/jcnunez/fst_brent/annotations/all_haplo_states.Rdata")

haplo_SNPs = c("X2R_20550549_SNP",
               "X2R_20551182_SNP",
               "X2R_20551519_SNP",
               "X2R_20550323_SNP",
               "X2R_20552286_SNP",
               "X2R_20551486_SNP",
               "X2R_20551129_SNP",
               "X2R_20550070_SNP",
               "X2R_20550046_SNP",
               "X2R_20550403_SNP")

analysis_o = 
foreach(i = haplo_SNPs,
        .combine = "rbind",
        .errorhandling = "remove")%do%{
          
          snpInt <- which(names(all_haplo_states) == i)
          sp70 <- which(names(all_haplo_states) == "SP70")
          pheno <- which(names(all_haplo_states) == "CTF")
          
          all_haplo_states[,c(sp70, snpInt, pheno)] %>%
            filter(!is.na(SP70)) -> tmp
          
          apply(tmp[,1:2], 
                1, function(x) paste(x, collapse = ";") ) -> joint_geno_tmps
          
          joint_geno_tmps %>% table -> tabsnps
          ####
          data.frame(tmp, joint_geno = joint_geno_tmps) %>%
            group_by(joint_geno) %>%
            summarise(mean = mean(CTF, na.rm = T)) ->
            tmp_ctf
          
          ####
          data.frame(
            snp_id = i,
            n00 = tabsnps[1],
            n02 = tabsnps[2],
            n20 = tabsnps[3],
            n22 = tabsnps[4],
            ntot = sum(tabsnps),
            ctmax0 = tmp_ctf$mean[1],
            ctmax0_mut = tmp_ctf$mean[2],
            ctmax2 = tmp_ctf$mean[3],
            ctmax2_mut = tmp_ctf$mean[4],
            delta_0 = tmp_ctf$mean[2] - tmp_ctf$mean[1] ,
            delta_2 = tmp_ctf$mean[4] - tmp_ctf$mean[3] 
          ) -> o
          
          return(o)
          
        }

left_join(analysis_o, annotaions_snps) -> 
  haplo_muts_annotated

save(haplo_muts_annotated, file = "haplo_muts_annotated.Rdata")


haplo_muts_annotated %>%
  ggplot(aes(
    x=snp_id,
    y=delta_0,
    fill=type
  )) + geom_bar(stat = "identity") + theme_bw() +
  ylim(-0.305,0.2)->
  performance_plot0
ggsave(performance_plot0, file = "performance_plot0.pdf",
       w=4, h = 3)

haplo_muts_annotated %>%
  ggplot(aes(
    x=snp_id,
    y=delta_2,
    fill=type
  )) + geom_bar(stat = "identity") + theme_bw() +
  ylim(-0.305,0.2)->
  performance_plot2
ggsave(performance_plot2, file = "performance_plot2.pdf",
       w=4, h = 3)

#####




library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)
library(foreach)
library(data.table)
library(gmodels)

library(SeqArray)
library(gdsfmt)
library(SNPRelate)


reco_windows_f="/gpfs2/scratch/jcnunez/fst_brent/simulations/reco_wins.master.txt"
reco_windows = fread(reco_windows_f)

real_dat <-"/netfiles02/lockwood_lab/IntrogressionProject/FET_output/Fisher_test_results_All_crosses.txt"
real_FET <- fread(real_dat)
names(real_FET)[1] = "chr"

### Explore whole dataset
real_FET %>%
  filter(sk_all_f.test_pval < 5e-08) ->
  boinfe_cutoff
dim(boinfe_cutoff)

boinfe_cutoff %>%
  filter(chr == "2R") %>%
  filter(sk_all_f.test_pval < 5e-08) %>%
  filter(POS > 18e6 & POS < 22e6)
#19375612-20983602

boinfe_cutoff %>%
  group_by(chr) %>%
  summarise(N=n())

##
boinfe_cutoff %>%
  filter(chr == "2R") %>%
  filter(sk_all_f.test_pval < 5e-08) %>%
  filter(POS > 18e6 & POS < 22e6) #-> loc2R_win
  
#### MAKE GOWINDA FILES
real_FET %>%
  dplyr::select(chr,   POS) ->
  all_SNP

write.table(all_SNP, 
            file = "/netfiles02/lockwood_lab/IntrogressionProject/GOWINDA_sets/GOWINDA_setsall_SNP.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

###
boinfe_cutoff %>%
  dplyr::select(chr,   POS) ->
  boinfe_cutoff_SNP

write.table(boinfe_cutoff_SNP, 
            file = "/netfiles02/lockwood_lab/IntrogressionProject/GOWINDA_sets/GOWINDA_GENOMEWIDE_outliers_SNP.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

###
#19375612-20983602
boinfe_cutoff_SNP %>%
  filter( (chr == "2R" & POS > 19115753 & POS < 20615753) | (chr == "X" & POS > 15123410 & POS < 16123410) ) %>%
  dplyr::select(chr,   POS) ->
  windows_interest

write.table(filter(windows_interest, chr == "2R"), 
            file = "/netfiles02/lockwood_lab/IntrogressionProject/GOWINDA_sets/GOWINDA_2R_window_SNP.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

write.table(filter(windows_interest, chr == "X"), 
            file = "/netfiles02/lockwood_lab/IntrogressionProject/GOWINDA_sets/GOWINDA_X_window_SNP.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

####
d=data.frame(xmin=c(19115753/1e6,15123410/1e6), 
             xmax=c(20615753/1e6,16123410/1e6), 
             ymin=c(0,0), 
             ymax=c(15,15), 
             chr=c('2R','X')
)


  ggplot() +
    geom_point(
    data=  FET_plot,
    aes(
    x=POS/1e6,
    y=-log10(sk_all_f.test_pval),
    color = -log10(sk_all_f.test_pval) > 7.5))  +
  geom_rect(data=d, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
            fill = "blue", alpha = 0.2) +
  ylim(0,15) +
  theme_minimal() +
    scale_color_manual(values = c("black", "purple")) +
  theme(legend.position = "none") +
  facet_grid(~chr, scales = "free_x") ->
  FET.pplot

ggsave(FET.pplot, file = "FET.pplot.pdf", w= 6, h = 3)

####### Concordant analyses
#### --> bring in the seasonal data from other datasets

####
reco_windows_x_2R =
  reco_windows %>%
  filter(V1 %in% c("2R","X"))

##
concordance = 
  foreach(i=1:dim(reco_windows_x_2R)[1], .combine = "rbind")%do%{
    
    reco_windows_x_2R[i,] -> guide.tmp
    
    chr.i=guide.tmp$V1
    sta=guide.tmp$V3
    end=guide.tmp$V4
    
    tmp =
      real_FET %>%
      filter(chr == chr.i ) %>%
      filter(POS >sta & POS < end )
    
    N = dim(tmp)[1]
    
    tmp %>%
      summarise(
    across( c("skf1_f.test_pval", "skf2_f.test_pval", "skf3_f.test_pval",
            "vt8f1_f.test_pval", "vt8f2_f.test_pval", "vt8f3_f.test_pval"), 
            ~ .x < 0.05 ) ) %>% as.data.frame() %>% rowSums()/6 -> vec.tmp
    
    data.frame(
      chr=chr.i,
      sta=sta,
      end=end,
      me_val = sum(vec.tmp==1),
      N = N
    )
  
  }


concordance %>%
  ggplot(aes(
    xmin=(sta+end)/2e6,
    xmax=(sta+end)/2e6,
    ymin=0,
    ymax=1,
    color=me_val/N
  )) + geom_rect() +
  theme_classic() +
  scale_color_gradient2(mid = "grey", high = "purple") +
  theme(legend.position = "bottom") +
  facet_grid(~chr, scales = "free_x") ->
  comcorc_plot

ggsave(comcorc_plot, file = "comcorc_plot.pdf",
       w= 6, h = 3.0)

concordance %>%
  mutate(per = me_val/N) %>%
  summarise( mean = mean(per)*100)

  
concordance %>%
mutate(per = me_val/N) %>%
  group_by(chr) %>%
  summarise( mean = mean(per)*100)

concordance %>%
  filter( (chr == "2R" & sta > 19115753 & end < 20615753) | 
            (chr == "X" & sta > 15123410 & end < 16123410) ) %>%
  mutate(per = me_val/N) %>%
  group_by(chr) %>%
  summarise( mean = mean(per)*100,
             N = sum(N), me_S = sum(me_val)
  )



####  
meta <- "/netfiles/nunezlab/D_melanogaster_resources/Datasets/2024.Nunez_et_al_Genetics/DEST_Cville_SNP_Metadata.Rdata"
metdat <- get(load(meta))
metdat$POS = as.double(metdat$POS)

####
real_FET %>% 
  filter(chr %in% c("2R", "X")) %>%
  filter(sk_all_f.test_pval <  5e-08) %>%
  filter(across( c("skf1_f.test_pval", "skf2_f.test_pval", "skf3_f.test_pval",
                   "vt8f1_f.test_pval", "vt8f2_f.test_pval", "vt8f3_f.test_pval"), 
                 ~ .x < 0.05 )) -> putative_targets

putative_targets %>% 
  left_join(metdat) ->
  putative_targets.annot

putative_targets.annot %<>%
  mutate(simple_eff = 
           case_when(effect %in% c("3_prime_UTR_variant","5_prime_UTR_variant") ~ "UTR",
                     effect %in% c("downstream_gene_variant","intergenic_region","upstream_gene_variant") ~ "intergenic",
                     effect %in% c("missense_variant&splice_region_variant") ~ "missense_variant",
                     effect %in% c("non_coding_transcript_exon_variant","splice_region_variant&non_coding_transcript_exon_variant") ~ "non-coding",
                     effect %in% c("splice_region_variant&intron_variant") ~ "intron_variant",
                     TRUE ~ effect
                     )) 

putative_targets.annot %>%
  group_by(chr, simple_eff) %>%
  summarise(N = n())

putative_targets.annot %>%
  filter(FeatureType == "transcript") %>%
  group_by(chr, GeneName) %>%
  summarise(N = n()) %>% 
  group_by(chr) %>% summarise(N = n()) 

putative_targets.annot %>%
  filter( (chr == "2R" & POS > 19115753 & POS < 20615753) | (chr == "X" & POS > 15123410 & POS < 16123410) ) %>%
  filter(!is.na(simple_eff)) %>%
  filter(simple_eff != "non-coding") %>%
  ggplot(aes(
    x=POS/1e6,
    y=as.factor(simple_eff),
    color = simple_eff,
    size= -log10(sk_all_f.test_pval)
  )) + geom_point() +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_grid(~chr, scales = "free_x") ->
  identiplot

ggsave(identiplot, file = "identiplot.pdf",
       w= 8, h = 3.0)

###
### GFF files
#gtf_dat <- fread("/netfiles/nunezlab/D_melanogaster_resources/Annotations/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gtf")
gtfCGS_dat <- fread("/netfiles/nunezlab/D_melanogaster_resources/Annotations/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic_CGStyle.gtf")

gtfCGS_dat %>%
  filter( (V1 == "2R" & V4 > 19115753 & V5 < 20615753) | (V1 == "X" & V4 > 15123410 & V5 < 16123410) ) %>% 
  filter(V3 == "CDS") ->
  genes_in_area

########## ****
##############################
#################### CHECKPOINT
########################################
save(putative_targets.annot, genes_in_area,
     file = "Brent_checkpoint_coloc.Rdata")




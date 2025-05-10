### libraries
library(data.table)
library(tidyverse)
library(foreach)

files = 
system("ls /gpfs2/scratch/jcnunez/fst_brent/GWAS_EB_OK/GWAS_out/batch.*.txt", 
       intern = T)

o = foreach(i=files,
            .combine = "rbind")%do%{
    tmp <- fread(i)
            }

####
files2 = 
  system("ls /gpfs2/scratch/jcnunez/fst_brent/GWAS_EB_OK/GWAS_out_CC/batch.*.txt", 
         intern = T)

o2 = foreach(i=files2,
            .combine = "rbind")%do%{
              tmp <- fread(i)
            }

o2 %>%
ggplot(aes(
  x=POS,
  y=-log10(PVAL))) + geom_point() +
  facet_grid(~CHR) ->
  CC_gwas
ggsave(CC_gwas, file = "CC_gwas.png", w=9, h=4)

####
files3 = 
  system("ls /gpfs2/scratch/jcnunez/fst_brent/GWAS_EB_OK/GWAS_out_AA/batch.*.txt", 
         intern = T)

o3 = foreach(i=files3,
             .combine = "rbind")%do%{
               tmp <- fread(i)
             }

o3 %>%
  ggplot(aes(
    x=POS,
    y=-log10(PVAL))) + geom_point() +
  facet_grid(~CHR) ->
  AA_gwas
ggsave(AA_gwas, file = "AA_gwas.png", w=9, h=4)


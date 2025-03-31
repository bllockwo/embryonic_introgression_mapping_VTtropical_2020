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

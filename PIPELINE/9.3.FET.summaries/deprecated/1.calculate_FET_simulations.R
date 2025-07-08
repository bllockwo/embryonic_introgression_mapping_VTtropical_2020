#### FET tests

library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)

countfiles = system("ls /gpfs2/scratch/jcnunez/fst_brent/simulations/collect_COUNTS_intro", intern = T)

args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1])
print(k)

i = countfiles[k]

inner=get(load(
  paste("/gpfs2/scratch/jcnunez/fst_brent/simulations/collect_COUNTS_intro/", i, sep = "")
))


inner %>%
  filter(AFsim > 0.05) %>%
  group_by(chr, real_pos) %>%
  mutate( p_fet = fisher.test(matrix(c(Counts_sim_pool, Cov_sim_pool,
                               VT8_0_1_ad, VT8_0_1_dp), 
                               nrow = 2))$p.value) ->
  outer.o

P_FET = outer.o[,c("chr", "real_pos", "p_fet")]
### Now collect the outputs

out = "/gpfs2/scratch/jcnunez/fst_brent/simulations/FET_SiMS_Ps"
name = paste("FET_P", k, "Rdata", sep = ".")

loc_sav= paste(out, name, sep = "/")

save(P_FET, file = loc_sav)




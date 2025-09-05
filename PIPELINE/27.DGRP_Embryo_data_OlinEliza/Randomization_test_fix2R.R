library(tidyverse)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(adegenet)
library(reshape2)
library(FactoMineR)
require(gtools)
require(foreach)

###
args = commandArgs(trailingOnly=TRUE)
#create job id from array
j = as.numeric(args[1])

###
anovaFun <- function(m1, m2) {
  #  m1 <- t0; m2<- t1.year
  ll1 <- as.numeric(logLik(m1))
  ll2 <- as.numeric(logLik(m2))
  
  parameter <- abs(attr(logLik(m1), "df") -  attr(logLik(m2), "df"))
  
  chisq <- -2*(ll1-ll2)
  
  1-pchisq(chisq, parameter)
  
}

###
genofile <- seqOpen("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/dgrp2.gds.gz")
seqResetFilter(genofile)

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile)) %>%
  mutate(snp_id = paste(chr,pos, sep = "_"))

### Fix 2R
#15483473-16017443
###
win.bp = 1e4
step.bp = 1e4+1

## 53
data.table(
  start=seq(from=15483473, to=16017443-win.bp, by=step.bp),
  end=seq(from=15483473, to=16017443-win.bp, by=step.bp) + win.bp) -> guide_table


SELECT_WINDOW = guide_table[j,]

snps.dt %>% 
  filter(chr == "X") %>%
  filter(pos > SELECT_WINDOW$start & pos < SELECT_WINDOW$end) ->
  picks1
snps.dt %>% 
  filter(snp_id %in%
           c("2R_16439138")) ->
  picks2

totpicks = c(picks1$variant.id,picks2$variant.id )

seqResetFilter(genofile)
seqSetFilter(genofile,  
             variant.id=totpicks)                      

snpgdsGetGeno(genofile, snp.id=totpicks,
              verbose=TRUE, with.id=TRUE) ->
  GENO.matrix

GENO.matrix$genotype %>%
  as.data.frame() ->
  GENO.matrix.dt

colnames(GENO.matrix.dt) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), sep="_")
rownames(GENO.matrix.dt) <- seqGetData(genofile, "sample.id")

### load data
egg_data <- fread("OK_EB_data.Aug16.2026.txt")
all_phenos <- readRDS("/gpfs2/scratch/jcnunez/fst_brent/GWAS_EB_OK/wideform.fixed.phenotable.RDS")
###
wolb <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/wolbachia.status.txt")
inv <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/Inversion.status.txt")
names(wolb)[1] = "Line"
names(inv)[1] = "Line"

GENO.matrix.dt %>%
  names() -> names_to_clean

names_to_clean %>% table() %>% .[.== 1] %>% names -> cleanna

GENO.matrix.dt[,which(names(GENO.matrix.dt) %in% cleanna)] %>% 
  mutate(Line = rownames(.)) %>%
  select(Line, `2R_16439138`) -> ANCHOR

GENO.matrix.dt[,which(names(GENO.matrix.dt) %in% cleanna)] %>% 
  mutate(Line = rownames(.)) %>%
  select(!`2R_16439138`) -> OTHERs

OTHERs %>%
  melt(id = "Line") %>%
  left_join(ANCHOR) %>% 
  full_join(egg_data) %>% 
  full_join(wolb) %>% 
  full_join(inv) ->
  GENO.matrix.pheno

unique(GENO.matrix.pheno$variable) -> targets

hold2r = 
foreach(i=targets,
        .combine = "rbind", .errorhandling = "remove")%do%{
          
          message(i)
          GENO.matrix.pheno %>%
            filter(variable == i) %>%
            filter(!is.na(`2R_16439138`)) %>%
            filter(!is.na(value)) %>%
            filter(!is.na(Proportion)) ->
            tmp
          
          lm(Proportion~`2R_16439138`+value+`In(2L)t`+`In(2R)NS`+
               Infection_Status, data = tmp) -> res_bas
          
          lm(Proportion~`2R_16439138`*value+`In(2L)t`+`In(2R)NS`+
               Infection_Status, data = tmp) -> res_inter
          
          data.frame(
            type = "real",
            i=i,
            AF=mean(tmp$value/2),
            p=anovaFun(res_bas, res_inter)
          ) ->real
          
          
          #res %>% as.data.frame() %>% mutate(var = rownames(.)) %>%
          #  mutate(i) %>% mutate(type = "real") -> real
          
          perms=
            foreach(k=1:50,
                  .combine = "rbind", .errorhandling = "remove")%do%{
                    message(k)

                    lm(sample(Proportion)~`2R_16439138`+value+`In(2L)t`+`In(2R)NS`+
                         Infection_Status, data = tmp) -> res_bas
                    
                    lm(sample(Proportion)~`2R_16439138`*value+`In(2L)t`+`In(2R)NS`+
                         Infection_Status, data = tmp) -> res_inter
                    
                    data.frame(
                      type = "perm",
                      i=i,
                      AF=mean(tmp$value/2),
                      p=anovaFun(res_bas, res_inter)
                    ) -> perm       
                  }
          
          data.frame(
            type = "empirical_p",
            i=i,
            AF=mean(tmp$value/2),
            p=sum(real$p > perms$p)/length(perms$p)
          ) ->emps
          
          rbind(real, perms, emps)
          
        }

hold2r %>%
  filter(AF > 0.05 & AF < 0.95) %>%
  group_by(i, type) %>%
  summarize(uci_1 = quantile(p, 0.01, na.rm = T),
            uci_5 = quantile(p, 0.05, na.rm = T)) %>%
  separate(i, into = c("chr","pos"), sep = "_") -> o

o %>%
  ggplot(
    aes(
      x=as.numeric(pos),
      y=-log10(uci_5),
      color=type
    )
  ) + geom_line() ->
  testplot
ggsave(testplot, file = "testplot.pdf")

root = "/gpfs2/scratch/jcnunez/fst_brent/Fall_2025_revisions/Exploration_DGRP/fix_2R"
addr = paste(root, paste("set", j, "Rdata" , sep = "."), sep = "/")

save(o, file =addr )




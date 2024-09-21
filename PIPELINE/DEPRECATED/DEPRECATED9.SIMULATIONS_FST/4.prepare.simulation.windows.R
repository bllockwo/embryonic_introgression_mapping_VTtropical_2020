##### Calculate FST in Introgression Lines

library(tidyverse)
library(magrittr)
library(poolfstat)

library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(data.table)
library(gtools)

library(foreach)
library(doMC)
library(doParallel)
registerDoMC(4)

#### REcombination windowing
filtering.dt <- get(load("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/SNP.filtering.guide.Rdata"))

filtering.dt %<>%
  mutate(simply_cm = round(cm_mb))

##filtering.dt %>% .$simply_cm %>% table
##filtering.dt %>%
##  ggplot(aes(
##    x=pos,
##    y=simply_cm,
##  )) + geom_line() + facet_wrap(~chr) -> CMbins
##ggsave(CMbins, file = "CMbins.pdf")


filtering.dt %>% group_by(chr) %>% arrange(pos) %>% 
  mutate(run=cumsum(c(0,diff(simply_cm))!=0)) %>% 
  subset(simply_cm!=2) %>% 
  group_by(chr,run) %>%
  summarise(level=max(simply_cm), start=min(pos), end=max(pos)) %>%
  select(-run) %>%
  mutate(Lenght = abs(start-end)) %>%
  mutate(rate_rho = level*1e-8) -> reco_wins
  #reco_wins$rate_rho %>% mean  

write.table(reco_wins, 
            file = "reco_wins.master.txt", 
            append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")

  #reco_wins %>% 
  #  group_by(chr)%>%
  #  summarize(n=n())
  
##reco_wins %>%
##  ggplot(aes(
##    x=(start+end)/2,
##    y=rate_rho,
##  )) + geom_line() + facet_wrap(~chr, nrow = 1) -> rhoins
##ggsave(rhoins, file = "rhobins.pdf", w = 12, h = 4)

####

##### objects needed for the analyses
samps = fread("/netfiles02/lockwood_lab/IntrogressionProject/Population_files/samp.guide.files.introg.txt")

###
genofile.path <- "/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds"

### load the genofile into memory
genofile <- seqOpen(genofile.path)

filtering.dt %>%
  filter(is.na(libs)) %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) ->
  chosen.snps
####

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile),
                      allele=seqGetData(genofile, "allele")) %>%
  separate(allele, into = c("ref_allele","alt_allele"), sep = ",")

#####
snps.dt %>%
  mutate(SNP_id = paste(chr, pos, sep = "_")) %>%
  right_join(chosen.snps, by = c("chr", "pos", "SNP_id") ) ->
  snps.dt

#### extrac
samps %<>%
  filter(Pop %in% c(
    "SK",
    "VT8"))

setDT(snps.dt)
snps.dt <- snps.dt[nAlleles==2]

######### window slicing
reco_windows=
foreach(k=1:dim(reco_wins)[1], 
        .combine = "rbind")%do%{

snps.dt %>%
  filter(chr == reco_wins$chr[k] & pos > reco_wins$start[k] & pos < reco_wins$end[k] ) ->
  snps.dt.win

#### ==
seqResetFilter(genofile)
seqSetFilter(genofile,
             sample.id=samps$Sample_id,
             variant.id=
               snps.dt.win$variant.id)
    
###
print("Create ad and dp objects")

ad <- seqGetData(genofile, "annotation/format/AD")
ad <- ad$data
dp <- seqGetData(genofile, "annotation/format/DP")

print("Create dat object")
dat <- ad/dp

#Add metadata ad
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") ,  sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")

dat %>% t %>%
  as.data.frame() %>%
  .[complete.cases(.),] %>%
  mutate(SNP_id = row.names(.)) %>% 
  mutate(F1_f = (SK_0_1 + VT8_0_1)/2) %>%
  mutate(VTparent_f = VT8_0_1) %>%
  separate(SNP_id, remove = F, into = c("chr","pos"), sep = "_") ->
  virtual_cross

virtual_cross$pos = as.numeric(virtual_cross$pos)

left_join(virtual_cross, snps.dt.win) %>% 
  mutate(adj_pos = pos-reco_wins[k,]$start ) -> virtual_cross

virtual_cross %>%
  dplyr::select(chr, pos, F1_f, VTparent_f, adj_pos) ->
  win.info
reco_wins[k,] -> win.meta

paste(reco_wins$chr[k], reco_wins$start[k], reco_wins$end[k],sep = "_") -> nameid

write.table(win.meta, 
            file = paste("win_meta/", nameid, ".meta.txt", sep = "")
            , append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")
write.table(win.info, 
            file = paste("win_data/", nameid, ".data.txt", sep = "")
            , append = FALSE, quote = F, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = F,
            col.names = F, qmethod = c("escape", "double"),
            fileEncoding = "")

}





####


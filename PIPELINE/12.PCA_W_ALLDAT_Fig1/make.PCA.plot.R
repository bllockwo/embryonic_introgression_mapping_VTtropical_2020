##### Build PCA with simulations

library(plyr)
library(tidyverse)
library(magrittr)
library(reshape2)
library(SeqArray)
library(FactoMineR)
library(factoextra)
library(data.table)
library(foreach)
library(ggrepel)

samps <- fread("/netfiles02/lockwood_lab/IntrogressionProject/Population_files/samp.guide.files.introg.txt")
samps %>% 
  filter(Pop %in% c("SK","VT8",
                    "SKF","VT8F",
                    "SKF2","VT8F2",
                    "SKF3","VT8F3"#,
                    #"ME_temp","PAN_trop"
  )) %>%
  mutate(type = case_when(
    #Pop %in% c("ME_temp", "PAN_trop") ~ "wild",
    Pop %in% c("SK") ~ "SK-parent",
    Pop %in% c("VT8") ~ "VT8-parent",
    Pop %in% c("VT8F", "VT8F2","VT8F3") ~ "F16-VT-mother",
    Pop %in% c("SKF", "SKF2","SKF3") ~ "F16-SK-mother"
  ))  -> samps.flt

#####
genofile.path <- "/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds"
genofile <- seqOpen(genofile.path)

sets =
list(all=c("2L", "2R", "3L", "3R", "X"),
chr2L="2L", 
chr2R="2R", 
chr3L="3L", 
chr3R="3R", 
chrX="X")

PCA_lists = list()
PCA_poly_lists = list()

PCA_all_OBJ = list()
PCA_poly_OBJ = list()

### Do PCAs
for(i in 1:length(sets)){
  
  chr_subset = sets[[i]]
  subset_name = names(sets)[i]
  message(paste(chr_subset, subset_name, sep = " "))
  
  seqResetFilter(genofile)
  snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                        pos=seqGetData(genofile, "position"),
                        variant.id=seqGetData(genofile, "variant.id"),
                        nAlleles=seqNumAllele(genofile),
                        missing=seqMissing(genofile))
  
  snps.dt <- snps.dt[nAlleles==2]
  seqSetFilter(genofile, variant.id=snps.dt$variant.id, sample.id = samps.flt$Sample_id)
  snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]
  
  seqSetFilter(genofile,
               snps.dt[chr%in%c(
                 chr_subset
               )][missing<.05]$variant.id)
  
  ### get allele frequency data
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")
  
  dat <- ad$data/dp
  dim(dat)  
  
  colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), sep="_")
  rownames(dat) <- seqGetData(genofile, "sample.id")
  
  ####
  dat %>%
    t() %>% 
    as.data.frame -> dat_t
  
  ## Some characterizations of AFs and subsequent filtering
  MeanAF=c()
  MinAF=c()
  MaxAF=c()
  
  apply(dat_t,
        1, FUN=mean, na.rm=TRUE ) -> MeanAF
  apply(dat_t,
        1, FUN=min, na.rm=TRUE ) -> MinAF
  apply(dat_t,
        1, FUN=max, na.rm=TRUE ) -> MaxAF
  
  cbind(dat_t, 
        MeanAF, MinAF, MaxAF) -> dat_t
  
  dat_t %>%
    .[which(MeanAF > 0.1),] ->
    #.[which(MeanAF != 1),] %>%
    #.[which(MeanAF != 0),] ->
    dat_t_flt
  
  #dat_t_flt %>% dim
  #dat_t_flt[grep("2L", rownames(dat_t_flt)),] %>% dim
  #dat_t_flt[grep("2R", rownames(dat_t_flt)),] %>% dim
  #dat_t_flt[grep("3L", rownames(dat_t_flt)),] %>% dim
  #dat_t_flt[grep("3R", rownames(dat_t_flt)),] %>% dim
  #dat_t_flt[grep("X", rownames(dat_t_flt)),] %>% dim
  
  dat_t_flt %>%
    dplyr::select(!c(MeanAF,MinAF,MeanAF,MinAF,MaxAF)) %>% 
    t %>%
    as.data.frame() ->
    dat_flt
  
  #### PCA
  random_snps = sample(colnames(dat_flt), 5000)
  
  dat_flt[,random_snps] %>%
    PCA(graph = F) ->
    pca.obj.all
  PCA_all_OBJ[[i]] = pca.obj.all
  
  ####
  pca.obj.all$ind$coord %>%
    as.data.frame() %>%
    mutate(Sample_id = rownames(.)) %>%
    left_join(samps.flt) %>%
    mutate(set_info = subset_name,
           var = "all"
           ) ->
    o
  PCA_lists[[i]] = o
  
  ####
  
  dat_t %>%
    .[which(.$MeanAF > 0.00 & .$MeanAF < 1.00),] %>%
    .[which(.$MinAF > 0.001 ),] %>%
    .[which(.$MaxAF < 0.999 ),] ->  ### This samples only polymorphic sites
    dat_t_poly
  
  dat_t_poly %>%
    dplyr::select(!c(MeanAF,MinAF,MeanAF,MinAF,MaxAF)) %>% 
    t %>%
    as.data.frame() ->
    dat_t_poly_flt
  
  #### PCA
  sampler = ifelse(dim(dat_t_poly_flt)[2] > 5000, 5000, dim(dat_t_poly_flt)[2])
  random_snps2 = sample(colnames(dat_t_poly_flt), sampler)
  
  dat_t_poly_flt[,random_snps2] %>%
    PCA(graph = F) ->
    pca.obj.poly
  
  PCA_poly_OBJ[[i]] = pca.obj.poly
  ####
  pca.obj.poly$ind$coord %>%
    as.data.frame() %>%
    mutate(Sample_id = rownames(.)) %>%
    left_join(samps.flt) %>%
    mutate(set_info = subset_name,
           var = "poly"
    ) ->
    o2
  PCA_poly_lists[[i]] = o2
}

pcas_all = do.call(rbind, PCA_lists)
pcas_poly = do.call(rbind, PCA_poly_lists)
####
names(PCA_all_OBJ) = sets
names(PCA_poly_OBJ) = sets
save(PCA_all_OBJ, file = "PCA_all_OBJ_eigeninfo.Rdata")
save(PCA_poly_OBJ, file = "PCA_poly_OBJ_eigeninfo.Rdata")

####
rbind(pcas_all, pcas_poly) ->
  pca_values_obj
save(pca_values_obj, file = "pca_values_obj.Rdata")

pca_values_obj %>%
  ggplot(
    aes(x= Dim.1, y=Dim.2,
        fill = type,
        shape = type)
  ) + geom_point(size = 3.0)  +
  facet_grid(var~set_info) +
  scale_shape_manual(values = 21:24) +
  theme_bw() + theme(legend.position="bottom")->
  pca_plot_intro

ggsave(pca_plot_intro, 
       file = "pca_plot_intro.pdf", 
       h = 4, w = 8)

pca_values_obj %>%
  ggplot(
    aes(x= Dim.3, y=Dim.4,
        fill = type,
        shape = type)
  ) + geom_point(size = 3.0)  +
  facet_grid(var~set_info) +
  scale_shape_manual(values = 21:24) +
  theme_bw() + theme(legend.position="bottom")->
  pca_plot_intro34

ggsave(pca_plot_intro34, 
       file = "pca_plot_intro34.pdf", 
       h = 4, w = 8)
######## 

eig_lists = list()
for(i in 1:length(PCA_all_OBJ)){
  PCA_all_OBJ[[i]]$eig[1:2,1:2] -> o
  eig_lists[[i]] = o
}

eig_poly_lists = list()
for(i in 1:length(PCA_all_OBJ)){
  PCA_poly_OBJ[[i]]$eig[1:2,1:2] -> o
  eig_poly_lists[[i]] = o
}


###### SLIDING WINDOW PCA ######
################################

chr_subset = c("2L", "2R", "3L", "3R", "X")

seqResetFilter(genofile)
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))

snps.dt <- snps.dt[nAlleles==2]
seqSetFilter(genofile, variant.id=snps.dt$variant.id, sample.id = samps.flt$Sample_id)
snps.dt[,af:=seqGetData(genofile, "annotation/info/AF")$data]
snps.dt[chr%in%c(
  chr_subset)] ->
    flt.snp.dt

win.bp <- 1e5
step.bp <- 5e4

setDT(flt.snp.dt)
setkey(flt.snp.dt, "chr")

## prepare windows
wins <- foreach(chr.i=c("2L","2R", "3L", "3R", "X"),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  tmp <- flt.snp.dt %>%
                    filter(chr == chr.i)
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                }

wins[,i:=1:dim(wins)[1]]
dim(wins)

#### LOOP OVER WINDOWs
euc_dist <- function(x1, x2){
  return(sqrt(sum((x1 - x2)^2)))
}

euc_analysis <- foreach(win.i=1:dim(wins)[1],
                   .errorhandling = "remove",
                   .combine = "rbind"
)%do%{
  message(paste(win.i, dim(wins)[1], sep=" / "))
  
  seqResetFilter(genofile)
  win.tmp <- flt.snp.dt %>%
    filter(chr == wins[win.i]$chr) %>%
    filter(pos %in% wins[win.i]$start:wins[win.i]$end)
  
  seqSetFilter(genofile,
               variant.id=win.tmp$variant.id, 
               sample.id = samps.flt$Sample_id)

  ### get allele frequency data
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")
  
  dat <- ad$data/dp
  #dim(dat)  
  
  colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position"), sep="_")
  rownames(dat) <- seqGetData(genofile, "sample.id")
  
  colMeans(dat, na.rm = T) -> filter.obj
  select_snps <- which(filter.obj > 0 & filter.obj < 1)
  
  ####
  dat[,select_snps] %>%
  PCA(graph = F) ->
    pca.obj.tmp
  ####
  pca.obj.tmp$ind$coord %>%
    as.data.frame() %>%
    mutate(Sample_id = rownames(.)) %>%
    left_join(samps.flt) ->
    pc.projs
  
  #Euclidean vs VT
  foreach(parent.i=pc.projs$Pop[c(4,8)],
          .errorhandling = "remove",
          .combine = "rbind"
  )%do%{
  foreach(samp.i=pc.projs$Pop[-c(4,8)],
          .errorhandling = "remove",
          .combine = "rbind"
  )%do%{
    pc.projs %>%
      filter(Pop == parent.i) -> parent
    parent$Dim.1 -> xp
    parent$Dim.2 -> yp
    pc.projs %>%
      filter(Pop == samp.i) -> F16
    F16$Dim.1 -> xo
    F16$Dim.2 -> yo
    euc_dist(c(xp, xo),c(yp, yo)) -> EUd
    
    o=
    data.frame(
      parent=parent.i,
      F16=samp.i,
      eud=EUd,
      chr = wins[win.i]$chr,
      win_start = wins[win.i]$start,
      win_end = wins[win.i]$end)
      return(o)
  }}

}

save(euc_analysis, file = "euc_analysis.Rdata")

euc_analysis %>%
  dcast(chr+win_start+win_end+F16~parent, 
        value.var = "eud") %>% 
  mutate(EU_Ratio = SK/VT8) %>%
  group_by(chr,win_start,win_end) %>%
  summarise(mean_EU_rat = mean(EU_Ratio)) ->
  dat_ucD

dat_ucD %>%
  ggplot(aes(
    x=(win_start+win_end)/2e6,
    y=log2(mean_EU_rat),
  )) +
  geom_hline(yintercept = 0) +
  geom_line(linewidth=0.4) +
  theme_bw() +
  facet_grid(.~chr, scales = "free_x") ->
  plot_euc

ggsave(plot_euc, file = "plot_euc.pdf",
       w=7, h = 2.5)
  
###
dat_ucD %>%
  filter(mean_EU_rat < 1) %>%
  group_by(chr) %>%
  summarise(N = n())
  

### investigate patterns in VCF files
#Load libraries
library(tidyverse)
library(vcfR)
library(data.table)
library(adegenet)
library(FactoMineR)
library(magrittr)
library(reshape2)
library(zoo)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(ggrepel)
library(lubridate)
library(foreach)
library(SeqArray)
library(doMC)
registerDoMC(8)
library(fastglm)
library(lmtest )
### Find out line genotypes

# Load the VA vcf
VAvcf <- read.vcfR(
  "/gpfs2/scratch/jcnunez/fst_brent/slice_VCFs/top2R_1MB_VA.recode.vcf.gz")
VAgl <- vcfR2genlight(VAvcf) 

###
tab(VAgl, NA.method = "asis") %>%
  as.data.frame() %>%
  t() ->
  VA_data

colnames(VA_data) -> VA_samps
data.frame(sampleid=VA_samps) %>%
  separate(sampleid, 
           remove = F, into = c("loc", 
                                "samp", 
                                "time"),
           sep = "_") ->
  VA_meta

VA_meta %<>%
  mutate( dat_spec = case_when(loc == "CM" ~ time,
                   loc %in% c("OW1","OW2") ~ NA
                   )) %>% 
  separate(dat_spec, into = c("MO", "DAY"), sep = 2) %>%
  mutate(Date_corr = as.Date(paste(DAY, MO, 2016), 
                             format = "%d %m %Y"))
  
VA_data %>%
  as.data.frame() %>%
  mutate(SNP_id = rownames(.)) %>%
  filter(SNP_id == "2R_20551633_SNP") %>%
  t() %>% as.data.frame() %>%
  mutate(sampleid = rownames(.)) %>%
  .[-dim(.)[1],] %>% 
  left_join(VA_meta)->
  SP70_data

SP70_data %>%
  group_by(MO) %>%
  filter(!is.na(Date_corr)) %>%
  summarise(N=n()*2,
            all.c = sum(as.numeric(`2R_20551633_SNP`))) %>%
  mutate(AF.t = prop.test(all.c, N)$estimate,
         AF.l = prop.test(all.c, N)$conf.int[1],
         AF.u = prop.test(all.c, N)$conf.int[2]
         ) -> MO.change

MO.change %>%
  filter(MO < 12) %>%
  ggplot(aes(
    x=as.numeric(MO),
    y=AF.t,
  )) + geom_point() + 
  geom_line() + theme_bw() +
  xlim(1,12) + ylim(0.3,0.7)->
  geno_AF_plot
ggsave(geno_AF_plot, file = "geno_AF_plot.pdf", 
       h = 3, w =6)



samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.csv")
samps %>%
  filter(city %in% "Charlottesville") %>% 
  filter(year %in% 2016:2018) %>% .$sampleId -> cville.SAMPStoUSE


### open GDS for common SNPs (PoolSNP)
genofile <- seqOpen("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/dest.all.PoolSNP.001.50.24Aug2024.ann.gds", allow.duplicate=T)

####
seqResetFilter(genofile)
snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     nAllele=seqGetData(genofile, "$num_allele"),
                     id=seqGetData(genofile, "variant.id"))                     
snp.dt <- snp.dt[nAllele==2]

seqSetFilter(genofile, snp.dt$id)
snp.dt[,global_af:=seqGetData(genofile, "annotation/info/AF")$data]

#####
getData <- function(chr="2L", 
                    start=14617051, end=14617051, 
                    samplesUse=samps$sampleId) {
  # chr="2L"; start=14617051; end=14617051; samplesUse=samps.seas$sampleI
  
  ### filter to target
  seqResetFilter(genofile)
  snp.tmp <- data.table(chr=chr, pos=start:end)
  setkey(snp.tmp, chr, pos)
  setkey(snp.dt, chr, pos)
  seqSetFilter(genofile, variant.id=snp.dt[J(snp.tmp), nomatch=0]$id, sample.id=samplesUse, verbose=T)
  
  ### get frequencies
  message("Allele Freqs")
  
  ad <- seqGetData(genofile, "annotation/format/AD")
  dp <- seqGetData(genofile, "annotation/format/DP")
  
  if(class(dp)[1]!="SeqVarDataList") {
    dp.list <- list()
    dp.list$data <- dp
    dp <- dp.list
  }
  
  af <- data.table(ad=expand.grid(ad$data)[,1],
                   dp=expand.grid(dp$data)[,1],
                   sampleId=rep(seqGetData(genofile, "sample.id"), dim(ad$data)[2]),
                   variant.id=rep(seqGetData(genofile, "variant.id"), each=dim(ad$data)[1]))
  
  ### tack them together
  message("merge")
  afi <- merge(af, snp.dt, by.x="variant.id", by.y="id")
  
  afi[,af:=ad/dp]
  
  ### calculate effective read-depth
  afis <- merge(afi, samps, by="sampleId")
  
  afis[chr=="X", nEff:=round((dp*nFlies - 1)/(dp+nFlies))]
  afis[chr!="X", nEff:=round((dp*2*nFlies - 1)/(dp+2*nFlies))]
  afis[,af_nEff:=round(af*nEff)/nEff]
  
  ### return
  afis[,c("sampleId", "af_nEff", "af", "nEff"), with=F]
}

Cville.data <- getData(chr="2R",
                       start=20551633,
                       end=20551633,
                       samplesUse=cville.SAMPStoUSE) %>%
  left_join(samps) %>%
  filter(year %in% 2016) 

Cville.data %>%
  group_by(MO=min_month) %>%
  summarise(AF.t = mean(af)) %>%
  mutate(source = "pool")-> pool_Ests

MO.change %>%
  dplyr::select(MO, AF.t) %>%
  mutate(source = "ind")-> ind_Ests

rbind(pool_Ests, ind_Ests) -> merged_ests

merged_ests %>%
  filter(MO > 6 & MO < 12) %>%
  dcast(MO~source, value.var = "AF.t") %>%
grangertest(pool~ind,data = .)

merged_ests %>%
  filter(MO > 6 & MO < 12) %>%
  ggplot(aes(
    x=as.numeric(MO),
    y=AF.t,
    linetype = source
  )) + geom_point() + 
  geom_line() + theme_bw() +
  xlim(1,12) ->
  geno_AF_plot2
ggsave(geno_AF_plot2, file = "geno_AF_plot2.pdf", 
       h = 3, w =6)

#### GET SAMPLE NAMES
#### GET SAMPLE NAMES
#### GET SAMPLE NAMES
#### GET SAMPLE NAMES
#### GET SAMPLE NAMES
#### GET SAMPLE NAMES
SP70_data %>%
  filter(`2R_20551633_SNP` %in% c(0)) %>%
  .$sampleid -> CC_samps_VA
write.table(CC_samps_VA,
            file = "CC_samps_VA.txt", append = FALSE, 
            quote = FALSE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
SP70_data %>%
  filter(`2R_20551633_SNP` %in% c(2)) %>%
  .$sampleid -> AA_samps_VA
write.table(AA_samps_VA,
            file = "AA_samps_VA.txt", append = FALSE, 
            quote = FALSE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

#####
# Load the VA vcf
WORLDvcf <- read.vcfR(
  "/gpfs2/scratch/jcnunez/fst_brent/slice_VCFs/top2R_1MB_World.recode.vcf.gz")
WORLDgl <- vcfR2genlight(WORLDvcf) 

###
tab(WORLDgl, NA.method = "asis") %>%
  as.data.frame() %>%
  t() ->
  WORLD_data

WORLD_data %>%
  as.data.frame() %>%
  mutate(SNP_id = rownames(.)) %>%
  filter(SNP_id == "2R_20551633_SNP") %>%
  t() %>% as.data.frame() %>%
  mutate(sampleid = rownames(.)) %>%
  .[-dim(.)[1],] ->
  SP70_dataWorld

SP70_dataWorld %>%
  filter(`2R_20551633_SNP` %in% c(0)) %>%
  .$sampleid -> CC_samps_World
write.table(CC_samps_World,
            file = "CC_samps_World.txt", append = FALSE, 
            quote = FALSE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")
SP70_dataWorld %>%
  filter(`2R_20551633_SNP` %in% c(2)) %>%
  .$sampleid -> AA_samps_World
write.table(AA_samps_World,
            file = "AA_samps_World.txt", append = FALSE, 
            quote = FALSE, sep = " ",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = FALSE, qmethod = c("escape", "double"),
            fileEncoding = "")

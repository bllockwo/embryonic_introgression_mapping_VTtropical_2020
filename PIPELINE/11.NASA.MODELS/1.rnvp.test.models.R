#### Cross-Enrichment 
#### 

rm(list = ls())

### libraries
library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(magrittr)
library(doMC)
registerDoMC(4)
library(SeqArray)
library(lubridate)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)
k=as.numeric(args[1]) ## 1-6
#k=6

####
#### ==> import FET -- -from Brent's data
intro.fet <- fread("/netfiles02/lockwood_lab/IntrogressionProject/QTL_mapping_Brent/Fisher_test_results_All_crosses.txt")
names(intro.fet)[1:2] = c("chr","pos")


#####
sets <- data.table(mod=c(1:11),
                   start=c(0,  0,  0,  7, 15, 30, 60, 15, 45,  0,  0),
                   end=	 c(7, 15, 30, 15, 30, 60, 90, 45, 75, 60, 90))
sets

### open GDS
genofile <- seqOpen("/netfiles/nunezlab/Drosophila_resources/Datasets/Nunez2023_Inversion.Seasonal/dest.all.PoolSNP.001.50.10Mar2021.ann.gds", allow.duplicate=T)


### samps
samps <- fread("/netfiles/nunezlab/Drosophila_resources/Datasets/Nunez2023_Inversion.Seasonal/DEST_10Mar2021_POP_metadata.csv")
samps[,month:=tstrsplit(collectionDate, "/")[[1]]]
samps[nchar(month)==1, month:=paste("0", month, sep="")]
samps[,day:=tstrsplit(collectionDate, "/")[[2]]]
samps[nchar(day)==1, day:=paste("0", day, sep="")]
samps <- samps[set!="dgn"]
samps[,Date:=date(paste(year, month, day, sep="-"))]

####
models=
  c(
    "temp.var;5;5.Cville",
    "temp.var;1;5.Cville",
    "temp.propmin;10;5.Cville",
    "temp.propMax;10;5.Cville",
    "temp.propMax;6;5.Cville",
    "temp.min;6;5.Cville",
    "temp.min;3;5.Cville"
  )

base <- "/netfiles/nunezlab/Drosophila_resources/Datasets/Nunez2023_Inversion.Seasonal"


######
######
### ----> loac Cville object ... this is a constant througt
#foreach(mod = models)%do%{
  
#### select model
  mod = models[k]
  
  tmp.fil <- paste(base, mod,  paste(mod,".glmRNP.Rdata",sep = ""), sep = "/" )
  print(tmp.fil)
  mod.obj <- get(load(tmp.fil)) 
  
  ###############
  ### windows ###
  ###############
  # generate a master index for window analysis
  ### define windows
  win.bp <- 1e5
  step.bp <- 5e4
  
  setkey(mod.obj, "chr")
  
  ## prepare windows
  wins <- foreach(chr.i=c("2L","2R", "3L", "3R"),
                  .combine="rbind", 
                  .errorhandling="remove")%dopar%{
                    
                    tmp <- mod.obj %>%
                      filter(chr == chr.i)
                    
                    data.table(chr=chr.i,
                               start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                               end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                  }
  
  wins[,i:=1:dim(wins)[1]]
  
  dim(wins)
  #####
  setkey(mod.obj, chr, pos)
  head(mod.obj)
  
  ### start the summarization process
  win.out <- foreach(win.i=1:dim(wins)[1], 
                     .errorhandling = "remove",
                     .combine = "rbind"
  )%do%{
    
    message(paste(win.i, dim(wins)[1], sep=" / "))
    
    
    win.tmp <- mod.obj[J(data.table(chr=wins[win.i]$chr, 
                                    pos=wins[win.i]$start:wins[win.i]$end, 
                                    key="chr,pos")), nomatch=0]
    
    
    #### Calculate Z score
    win.tmp[,Z:=qnorm(p_lrt, 0, 1)]
    #### Calculate Z rnp score
    win.tmp[,rnpZ:=qnorm(rnp, 0, 1)]
    
    
    seqSetFilter(genofile, 
                 variant.id=unique(win.tmp$variant.id),
                 sample.id=samps[locality=="VA_ch"][year>= 2016 ]$sampleId)
    
    #obtain AFs 
    af <- seqGetData(genofile, "annotation/format/FREQ")
    f.hat <- data.table(fhat=colMeans(af[[2]], na.rm=T), 
                        variant.id=seqGetData(genofile, "variant.id"))
    
    #merge AFs with object
    win.tmp <- merge(win.tmp, f.hat, by="variant.id")
    win.tmp[,het:=2*fhat*(1-fhat)]
    
    #thrs <- expand.grid(sapply(c(1:9), function(x) x*10^c(-5:-1)))[,1]
    
    pr.i <- c(
      0.05
    )
    
    #tmpo <- foreach(
    # pr.i=thrs, 
    # .errorhandling="remove", 
    # .combine="rbind")%do%{
    win.tmp %>% 
      filter(!is.na(rnp), pr.i == pr.i ) %>%
      group_by(perm, chr , variable, mod) %>%
      summarise(pos_mean = mean(pos),
                pos_mean = mean(pos),
                pos_min = min(pos),
                pos_max = max(pos),
                win=win.i,
                pr=pr.i,
                rnp.pr=c(mean(rnp<=pr.i)),
                rnp.binom.p=c(binom.test(sum(rnp<=pr.i), 
                                         length(rnp), pr.i)$p.value),
                wZa=sum(het*Z)/(sqrt(sum(het^2))),
                wZa.p=pnorm(sum(het*Z)/(sqrt(sum(het^2))), lower.tail=T),
                rnp.wZa=sum(het*rnpZ)/(sqrt(sum(het^2))),
                rnp.wZa.p=pnorm(sum(het*rnpZ)/(sqrt(sum(het^2))), lower.tail=T),
                min.p.lrt=min(p_lrt),
                min.rnp=min(rnp),
                nSNPs = n(),
                sum.rnp=sum(rnp<=pr.i),
      ) %>%
      mutate(
        model.pop = models[k],
        perm_type=ifelse(perm==0, "real","permuted"),
        invName=case_when(
          chr=="2L" & pos_min >	2225744	 & pos_max < 13154180	 ~ "2Lt",
          chr=="2R" & pos_min >	15391154 & pos_max < 	20276334 ~ 	"2RNS",
          chr=="3R" & pos_min >	11750567 & pos_max < 	26140370 ~ 	"3RK",
          chr=="3R" & pos_min >	21406917 & pos_max < 	29031297 ~ 	"3RMo",
          chr=="3R" & pos_min >	16432209 & pos_max < 	24744010 ~ 	"3RP",
          chr=="3L" & pos_min >	3173046	 & pos_max < 16308841	 ~ "3LP",
          TRUE ~ "noInv"
        )) -> win.out
    
    return(win.out)
    #}
    #tmpo
  }
  
###
  ### save
  out_folder <- "/gpfs1/home/j/c/jcnunez/scratch/Brent_Introgression/rnpv_analysis"
  
  message(paste(out_folder, "/Window_analysis_", models[k] , ".Rdata", sep=""))
  
  save(win.out, 
       file=paste(out_folder, "/Window_analysis_", models[k] , ".Rdata", sep=""))

  
###  
  

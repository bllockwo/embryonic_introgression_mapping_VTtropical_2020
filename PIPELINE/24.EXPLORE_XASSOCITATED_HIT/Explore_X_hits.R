#### Explore X-associated SNP

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
library(fastglm)
library(rnaturalearth)


####
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
  
  afis[chr=="X", nEff:=round((dp*nFlies)/(dp+nFlies - 1))]
  afis[chr!="X", nEff:=round((dp*2*nFlies)/(dp+2*nFlies - 1))]
  afis[,af_nEff:=round(af*nEff)/nEff]
  
  ### return
  afis[,c("sampleId", "af_nEff", "af", "nEff"), with=F]
}
#########
anovaFun <- function(m1, m2) {
  #  m1 <- t0; m2<- t1.year
  ll1 <- as.numeric(logLik(m1))
  ll2 <- as.numeric(logLik(m2))
  
  parameter <- abs(attr(logLik(m1), "df") -  attr(logLik(m2), "df"))
  
  chisq <- -2*(ll1-ll2)
  
  1-pchisq(chisq, parameter)
  
}
#####

####

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
  filter(chr == "X") %>%
  filter(POS >= 15123410 & POS <= 16123410) ->
  top_muts_X

top_muts_X %>%
  filter(pval.adj_skf1 < 0.01,
         pval.adj_skf2< 0.01,
         pval.adj_skf3 < 0.01,
         pval.adj_vt8f2 < 0.01,
         pval.adj_vt8f3 < 0.01,
         pval.adj_vt8f1 < 0.01)


top_muts_X %>%
  slice_min(pval.adj_all)
###
### open GDS
genofile <- seqOpen("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/dest.all.PoolSNP.001.50.8Jun2023.norep.AT_EScorrect.ann.gds", allow.duplicate=T)


#####
meta_git <- "https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.xa.csv"
samps <- fread(meta_git)
setDT(samps)

samps %>% 
  filter(city %in% 
           c("Charlottesville",
             "Odesa",
             "Yesiloz"
           )) %>%
  filter(set != "dgn") %>%
  filter(Recommendation == "Pass") %>%
  filter(!is.na(year)) %>%
  filter(!is.na(min_day)) %>%
  group_by(sampleId) %>%
  mutate(DATE_spec = as.Date(str_split(sampleId, "_")[[1]][5],
                             format = "%Y-%m-%d")) -> seasonal.samps

samps %>% 
filter(city %in% 
         c("Charlottesville"
         )) %>%
  filter(set != "dgn") %>%
  filter(Recommendation == "Pass") %>%
  filter(!is.na(year)) %>%
  filter(!is.na(min_day)) %>%
  group_by(sampleId) %>%
  mutate(DATE_spec = as.Date(str_split(sampleId, "_")[[1]][5],
                             format = "%Y-%m-%d")) -> Cville.seasonal.samps

samps %>% 
  filter(Recommendation == "Pass") %>%
  filter(!is.na(year)) %>%
  filter(!is.na(min_day)) %>%
  group_by(sampleId) %>%
  mutate(DATE_spec = as.Date(str_split(sampleId, "_")[[1]][5],
                             format = "%Y-%m-%d")) -> all.samps

####
seqResetFilter(genofile)
snp.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                     pos=seqGetData(genofile, "position"),
                     nAllele=seqGetData(genofile, "$num_allele"),
                     id=seqGetData(genofile, "variant.id"))                     
snp.dt <- snp.dt[nAllele==2]

seqSetFilter(genofile, snp.dt$id)
snp.dt[,global_af:=seqGetData(genofile, "annotation/info/AF")$data]

#### STOP! check for Robjects before re-generating data <---
#### STOP! check for Robjects before re-generating data <---
#### STOP! check for Robjects before re-generating data <---
#### STOP! check for Robjects before re-generating data <---

#### TEST for clinality
#### TEST for clinality
#### TEST for clinality
#### TEST for clinality
#### TEST for clinality
### samps

cline_test = 
foreach(snp=top_muts_X$POS,
        .combine = "rbind",
        .errorhandling = "remove")%do%{
          
          
          all.dataX <- getData(chr="X",
                               start=snp,
                               end=snp,
                               samplesUse=all.samps$sampleId) %>%
            left_join(all.samps) %>%
            mutate(year_pop = paste(year, city, sep = ":"))     
          
          all.dataX %>%
            group_by(city, cluster2.0_k4, continent) %>%
            summarise(
              m.AF = mean(af, na.rm = T),
              m.Lat = mean(lat),
              m.Long = mean(long)
            )-> cline_summaries
          
          cline_summaries %>% 
            filter(cluster2.0_k4 == 4) ->
            X_lat_NAm
          
          cor.test(~m.Lat+m.AF, data =  filter(X_lat_NAm, continent == "North_America")) -> o1
          cor.test(~m.Lat+m.AF, data =  filter(cline_summaries, continent == "Europe" & cluster2.0_k4 != 4)) -> o2
          cor.test(~m.Long+m.AF, data =  filter(cline_summaries, continent == "Europe" & cluster2.0_k4 != 4)) -> o3
          
          data.frame(POS=snp,
          corNA=o1$estimate,
          pNA=o1$p.value,
          corEU_lat=o2$estimate,
          pEU_lat=o2$p.value,
          corEU_long=o3$estimate,
          pEU_long=o3$p.value)
        }

save(cline_test, file = "cline_test.X.Rdata")
####
weather <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/Weather_NASA_power/WeatherData.mgarvin.csv")
weather %>%
  mutate(Date = as.Date(paste(YEAR,MO,DY, sep = "-"))) -> weather.dat

weather_samps =
  foreach(WinSlide = seq(from=15, to = 90, by = 15),
          .combine = "rbind")%do%{
            
            ######
            foreach(i=1:dim(seasonal.samps)[1], 
                    .combine = "rbind")%do%{
                      
                      print(paste(WinSlide, i, sep = "|"))
                      
                      samp.i <- seasonal.samps[i,]
                      
                      bdate=samp.i$DATE_spec-WinSlide
                      fdate=samp.i$DATE_spec
                      weather.dat %>%
                        filter(city == samp.i$city) %>%
                        filter(Date >= bdate & Date <= fdate) %>%
                        summarize(
                          RH2M_mean =  mean(RH2M),
                          T2M_mean = mean(T2M),
                          PRECTOTCORR_mean = mean(PRECTOTCORR),
                          
                          RH2M_sd =  sd(RH2M),
                          T2M_sd = sd(T2M),
                          PRECTOTCORR_sd = sd(PRECTOTCORR),
                          
                          RH2M_min =  min(RH2M),
                          T2M_min = min(T2M),
                          PRECTOTCORR_min = min(PRECTOTCORR),
                          
                          RH2M_max =  max(RH2M),
                          T2M_max = max(T2M),
                          PRECTOTCORR_max = max(PRECTOTCORR),
                          
                          T2M_5c = sum(T2M < 5),
                          T2M_32c = sum(T2M > 32)
                          
                        ) -> weather.slice
                      data.frame(
                        samp.i,
                        WinSlide = WinSlide,
                        weather.slice)
                      
                    } ### inner
            ######
          } ### outer

save(weather_samps, file = "weather_samps.Rdata")
####
#####
#### note ... USE ONLY CVILLE!!!
X_exploration = 
foreach(POS.i = top_muts_X$POS,
        .combine = "rbind")%do%{
        
          SAMPStoUSE = Cville.seasonal.samps$sampleId
          
          data <- getData(chr="X",
                          start=POS.i,
                          end=POS.i,
                          samplesUse=SAMPStoUSE)
          setkey(data, sampleId)
          data <- data[!duplicated(data)]
          data <- data %>% filter(!is.na(af_nEff))
          
          mod_o =
            foreach(Slide_K = seq(from=15, to = 90, by = 15),
                    .combine = "rbind", .errorhandling = "remove")%do%{
                      foreach(vars = c("RH2M_mean", "T2M_mean", "PRECTOTCORR_mean",
                                       "RH2M_sd", "T2M_sd", "PRECTOTCORR_sd",
                                       "RH2M_min", "T2M_min", "PRECTOTCORR_min",
                                       "RH2M_max", "T2M_max", "PRECTOTCORR_max",
                                       "T2M_5c", "T2M_32c"
                      ),
                      .combine = "rbind")%do%{
                        
                        #print(paste(perm, Slide_K, vars, sep = "|"))      
                        
                        perm = 0
                        
                        filter_samps_Dat = weather_samps %>%
                          filter(WinSlide == Slide_K)
                        
                        data.in = left_join(data, samps, by="sampleId") %>%
                          left_join(filter_samps_Dat, 
                                    by="sampleId")
                        
                        data.in %>%
                          select(year.x, city.x, nEff, af_nEff, WinSlide, 
                                 RH2M_mean, T2M_mean, PRECTOTCORR_mean,
                                 RH2M_sd, T2M_sd, PRECTOTCORR_sd,
                                 RH2M_min, T2M_min, PRECTOTCORR_min,
                                 RH2M_max, T2M_max, PRECTOTCORR_max,
                                 T2M_5c, T2M_32c
                          ) %>%
                          melt(id=c("year.x", "city.x", "nEff", "WinSlide","af_nEff")) %>%
                          mutate(year_pop = paste(year.x, city.x, sep = ":")) ->
                          data.in
                        
                        data.in %>% filter(variable == vars) -> tmp
                        y <- tmp$af_nEff
                        glm.method <- 0		
                        
                        if(perm == 0){
                          X.year.noperm<- model.matrix(~as.factor(year_pop), tmp)
                          t1.year.noperm <- fastglm(x=X.year.noperm, y=y, 
                                                    family=binomial(), weights=tmp$nEff, method=glm.method)
                          
                          X.year.var <- model.matrix(~as.factor(year_pop)+value, tmp)
                          t1.year.var <- fastglm(x=X.year.var, y=y, family=binomial(),
                                                 weights=tmp$nEff, method=glm.method)
                          
                          o =
                            data.table(
                                       POS=POS.i,
                                       mod=vars,
                                       AIC=c(AIC(t1.year.var)),
                                       b_temp=last(t1.year.var$coef), 
                                       se_temp=last(t1.year.var$se),
                                       p_lrt=anovaFun(t1.year.noperm, t1.year.var),
                                       Timewindow = Slide_K,
                                       perm=perm)
                          
                          
                        }
                        print(paste(POS.i, perm, Slide_K, vars, o$p_lrt, sep = "|"))
                        return(o)  
                      }
                    }
          
        }

save(X_exploration, file = "X_exploration.seas.Rdata")
###
###
###
load("X_exploration.seas.Rdata")
load("weather_samps.Rdata")
load("cline_test.X.Rdata")
####
###
###
###

left_join(cline_test, X_exploration) %>%
  filter(b_temp != 0) %>%
  group_by(POS) %>%
  slice_min(p_lrt, n = 1) ->
  dat.for.plot

dat.for.plot %>%
  ggplot(
   aes( x=-log10(pNA), 
    y=-log10(p_lrt))
  ) + geom_point(size = 3.0) + 
  geom_vline(xintercept = -log10(0.05)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_text(data = filter(dat.for.plot, pNA < 0.05 & p_lrt < 0.05),
            aes(label = POS)) + theme_bw()->
  coloc_x

ggsave(coloc_x, file = "coloc_x.pdf",
       h= 5, w = 5)
  
###
dat.for.plot %>% 
  filter(pNA < 0.05 & p_lrt < 0.05) ->
  topXsnps
#PRECTOTCORR_max -- 90 days
#X pos = 15602941 
#X pos = 15648874  (Odz03; sog regulator)
dat.for.plot %>%
  filter(POS == 15648874) %>% as.data.frame()

X_exploration %>%
  filter(POS == 15602941) %>%
  filter(p_lrt < 0.05)


#### Explore Cross
genofile_cross <- seqOpen("/netfiles02/lockwood_lab/IntrogressionProject/SNPcalling_output/BLockIntro.PoolSeq.PoolSNP.001.5.test.ann.gds")
####seqResetFilter(genofile_cross)

snps.CROSS <- data.table(chr=seqGetData(genofile_cross, "chromosome"),
                         pos=seqGetData(genofile_cross, "position"),
                         variant.id=seqGetData(genofile_cross, "variant.id"),
                         alleles=seqGetData(genofile_cross, "allele"),
                         nAlleles=seqNumAllele(genofile_cross))

#15648874
## NOTE that top SNP is --> 15602941
## NOTE that Odz03 is --> 15648874
snps.CROSS %>%
  filter(chr == "X") %>%
  filter(pos %in% 15602941)
### check the haplotype

snps.CROSS %>%
  filter(chr == "X") %>%
  filter(pos %in% 15602941) %>%
  filter(nAlleles == 2) ->
  haplo.CROSS

seqSetFilter(genofile_cross, variant.id=haplo.CROSS$variant.id, 
             sample.id = c("SKF2_2_1",
                           "SKF3_3_1",
                           "SKF_1_1",
                           "SK_0_1",
                           "VT8F2_2_1",
                           "VT8F3_3_1",
                           "VT8F_1_1",
                           "VT8_0_1"))


message("Allele Freqs")
ad.cross <- seqGetData(genofile_cross, "annotation/format/AD")$data
dp.cross <- seqGetData(genofile_cross, "annotation/format/DP")

dat.cross <- ad.cross/dp.cross
colnames(dat.cross) <- paste(seqGetData(genofile_cross, "chromosome"), seqGetData(genofile_cross, "position"), sep="_")
rownames(dat.cross) <- seqGetData(genofile_cross, "sample.id")

dat.cross %>% 
  as.data.frame() %>%
  mutate(sampleid = rownames(.)) %>%
  separate(sampleid,
           into = c("pop","rep1","rep2"),
           sep = "_") %>% 
  melt(id = c("pop", "rep1", "rep2")) %>%
  separate(variable, into = c("chr","pos"), sep = "_") %>%
  mutate(pos = as.numeric(pos)) ->
  dat.cross.o

dat.cross.o %>%
  separate(pop, remove =F, into = c("pop1","pop_etc"), sep = 2) %>%
  group_by(pop1, pos, parent=rep1==0) %>%
  summarize(m.AF = mean(value)) %>%
  ggplot(aes(
    x=pos,
    y=m.AF,
    shape=parent, color = pop1
  )) + geom_point(size = 4) + geom_hline(yintercept = 0.5) + theme_bw() ->
  delta_crossx
ggsave(delta_crossx, file = "delta_crossx.pdf",
       w=5, h = 4)

dat.cross.o %>%
  filter(pop == "SK" & rep1 == 0 & value == 1)

#pop rep1 rep2 chr      pos value
#1  SK    0    1   X 15602941     1 ## dm3-> 15496974
#2  SK    0    1   X 15607604     1 ## dm3-> 15501637
#3  SK    0    1   X 15847814     1 ## dm3-> 15741847

topXsnps %>%
  filter(POS %in% c(15602941, 15607604, 15847814 ))
######
######
## NOTE that top SNP is --> 15602941 .. with PRECTOTCORR_max 90
## NOTE that Odz03 is --> 15648874 ... with RH2M_sd 15

snppp = 15648874
win = 15
var = "RH2M_sd"
mod_o_perms =
  foreach(Slide_K = win, #seq(from=15, to = 90, by = 15),
          .combine = "rbind", .errorhandling = "remove")%do%{
            foreach(vars = c(var
            ),
            .combine = "rbind")%do%{
              
              #print(paste(perm, Slide_K, vars, sep = "|"))      
              
              #perm = 0
              
              filter_samps_Dat = weather_samps %>%
                filter(WinSlide == Slide_K)
              SAMPStoUSE = seasonal.samps$sampleId
              
              data <- getData(chr="X",
                              start=snppp,
                              end=snppp,
                              samplesUse=SAMPStoUSE)
              
              data.in = left_join(data, samps, by="sampleId") %>%
                left_join(filter_samps_Dat, 
                          by="sampleId")
              
              data.in %>%
                select(year.x, city.x, nEff, af_nEff, WinSlide, 
                       RH2M_mean, T2M_mean, PRECTOTCORR_mean,
                       RH2M_sd, T2M_sd, PRECTOTCORR_sd,
                       RH2M_min, T2M_min, PRECTOTCORR_min,
                       RH2M_max, T2M_max, PRECTOTCORR_max,
                       T2M_5c, T2M_32c
                ) %>%
                melt(id=c("year.x", "city.x", "nEff", "WinSlide","af_nEff")) %>%
                mutate(year_pop = paste(year.x, city.x, sep = ":")) ->
                data.in
              
              data.in %>% filter(variable == vars) -> tmp
              tmp = tmp[complete.cases(af_nEff),]
              
              y <- tmp$af_nEff
              glm.method <- 0		
              
              
              foreach(perm = 1:100, .combine = "rbind")%do%{
                
                value_perm = tmp[sample(1:dim(tmp)[1]), value]
                tmp_p = cbind(tmp, value_perm = value_perm)

                X.year.noperm<- model.matrix(~as.factor(year_pop), tmp)
                t1.year.noperm <- fastglm(x=X.year.noperm, y=y, 
                                          family=binomial(), weights=tmp$nEff, method=glm.method)
                
                X.year.var <- model.matrix(~as.factor(year_pop)+value_perm, tmp_p)
                t1.year.var <- fastglm(x=X.year.var, y=y, family=binomial(),
                                       weights=tmp_p$nEff, method=glm.method)
                
                o =
                  data.table(mod=vars,
                             AIC=c(AIC(t1.year.var)),
                             b_temp=last(t1.year.var$coef), 
                             se_temp=last(t1.year.var$se),
                             p_lrt=anovaFun(t1.year.noperm, t1.year.var),
                             Timewindow = Slide_K,
                             perm=perm)
                
                print(paste(perm, Slide_K, vars, perm, o$p_lrt, sep = "|"))
                return(o)  
              }
            }
          }

sum(mod_o_perms$p_lrt < 0.0153053937)/length(mod_o_perms$p_lrt)
sum(mod_o_perms$p_lrt < 0.1189294)/length(mod_o_perms$p_lrt)

##  0.03
## beat 97% of permutaions!

#######
####### SPATIAL ANALYSIS
####### SPATIAL ANALYSIS
####### SPATIAL ANALYSIS
####### SPATIAL ANALYSIS
####### SPATIAL ANALYSIS
####### SPATIAL ANALYSIS
####### SPATIAL ANALYSIS

### in dm3 chrX:15637404-15637404


####
  weather_Select = weather_samps %>%
    filter(WinSlide == win) %>%
    mutate(year_pop = paste(city,year, sep = ":"))
  
  
  Cville.dataX <- getData(chr="X",
                         start=snppp,
                         end=snppp,
                         samplesUse=seasonal.samps$sampleId) %>%
    left_join(seasonal.samps) %>%
    left_join(weather_Select) %>%
    mutate(year_pop = paste(year, city, sep = ":"))
  
  ggplot() +
    geom_smooth(
      data=mutate(Cville.dataX, type = "z.AFS"),
      aes(x=DATE_spec, y = af),  color = "black" ,
      se = F,   linewidth = 0.5, linetype = "dashed") + 
    facet_grid(city~year, scales = "free") + theme_bw()->
    perm_real_plotX
  ggsave(perm_real_plotX, 
         w = 7, h = 3,
         file = "perm_real_plotX.pdf")

  ggplot() +
    geom_smooth(
      data=mutate(Cville.dataX, type = "z.AFS"),
      aes(x=DATE_spec, y = PRECTOTCORR_max),
      se = F,   linewidth = 0.5, linetype = "solid") + 
    facet_grid(city~year, scales = "free") + theme_bw()->
    precipt_plotX
  ggsave(precipt_plotX, 
         w = 7, h = 3,
         file = "precipt_plotX.pdf")
  
  
  Cville.dataX %>%
  ggplot(aes(
    x=PRECTOTCORR_max,
    y=af,
  )) +
    geom_point( shape = 21, size = 4) +
    geom_smooth(method = "lm") + theme_bw() +
    facet_grid(~city)->
    weaplot
      ggsave(weaplot, 
         w = 9, h = 4,
         file = "PRECTOTCORR_max.pdf")

cor.test(~PRECTOTCORR_max+af, 
         data = filter(Cville.dataX, city == "Charlottesville"))
cor.test(~PRECTOTCORR_max+af, 
         data = filter(Cville.dataX, city == "Odesa"))
cor.test(~PRECTOTCORR_max+af, 
         data = filter(Cville.dataX, city == "Yesiloz"))

      ### bring in clinal weather      
weather <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/Weather_NASA_power/WeatherData.mgarvin.csv")
      weather %>%
        group_by(city) %>%
        summarize(
          m_RH2M=mean(RH2M, na.rm = T),
          m_T2M=mean(T2M, na.rm = T),
          m_PRECTOTCORR=mean(PRECTOTCORR, na.rm = T),
          
          v_RH2M=sd(RH2M, na.rm = T),
          v_T2M=sd(T2M, na.rm = T),
          v_PRECTOTCORR=sd(PRECTOTCORR, na.rm = T),
          
          max_RH2M=max(RH2M, na.rm = T),
          max_T2M=max(T2M, na.rm = T),
          max_PRECTOTCORR=max(PRECTOTCORR, na.rm = T),
          
          min_RH2M=min(RH2M, na.rm = T),
          min_T2M=min(T2M, na.rm = T),
          min_PRECTOTCORR=min(PRECTOTCORR, na.rm = T)
          
        ) -> weather.ag
      
      lnadmix <- get(load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/linear_admix/linear.Admix.Data.S2.Rdata"))
      lnadmix %>%
        group_by(city, source_pop, filter,admix.set) %>%
        summarise(mean.est = mean(Estimate),
                  sd.est = sd(Estimate)) %>%
        filter(filter == "silent") %>%
        filter(admix.set == "N.America" & source_pop == "AFRICA") ->
        o_dat

      o_dat %>%
        left_join(X_lat_NAm, 
                  by = c("city")) %>% 
        filter(!is.na(m.AF))-> data_merged
      
      data_merged %>%
        group_by(city) %>%
        left_join(weather.ag, by = "city") ->
        data_merged.wag
      
      lm(m.AF ~ mean.est, data=data_merged.wag) ->
        ancestry_model
      anova(ancestry_model)
      summary(ancestry_model) 
      
      data_merged.wag %>%
        as.data.frame() %>%
        mutate(AF_res = ancestry_model$residuals) ->
        data_merged.wag
      
      cor.test(~AF_res+m.Lat, data = data_merged.wag)
      
      select_cols = c("m_RH2M","m_T2M","m_PRECTOTCORR","v_RH2M",
                      "v_T2M", "v_PRECTOTCORR", "max_RH2M",
                      "max_T2M", "max_PRECTOTCORR", "min_RH2M",
                      "min_T2M", "min_PRECTOTCORR",
                      "AF_res","city")
      
      data_merged.wag[,select_cols] %>%
        .[complete.cases(.),] %>% 
        melt(id = c("AF_res","city")) ->
        AF_res.melt
      
      ####
      mods =
        foreach(k=0:100,
                .combine = "rbind",
                .errorhandling = "remove")%do%{
                  foreach(i=select_cols[1:12],
                          .combine = "rbind",
                          .errorhandling = "remove")%do%{
                            message(paste(k,i, sep = "|"))
                            
                            data=filter(AF_res.melt,
                                        variable == i) ->
                              dat.tmp
                            
                            if(k==0){
                              lm(AF_res ~ value, 
                                 data = dat.tmp) -> mod.o
                              
                              mod.o  %>%
                                summary ->
                                mod_tmp
                              
                              data.frame(
                                mod=i,
                                p_val=mod_tmp$coefficients[2,4],
                                est=mod_tmp$coefficients[2,1],
                                err=mod_tmp$coefficients[2,2],
                                perm=k
                              ) -> o
                              return(o)
                            } # k == 1
                            
                            if(k!=0){
                              
                              filter(AF_res.melt,
                                     variable == i) %>%
                                mutate(value.random = sample(value)) ->
                                tmp.random
                              
                              lm(AF_res ~ value.random, 
                                 data=tmp.random) %>%
                                summary ->
                                mod_tmp
                              
                              data.frame(
                                mod=i,
                                p_val=mod_tmp$coefficients[2,4],
                                est=mod_tmp$coefficients[2,1],
                                err=mod_tmp$coefficients[2,2],
                                perm=k
                              ) -> o
                              return(o)
                            } # if k > 1
                            
                          }#i
                }#k
      
      #save(mods, file = "clinal.mods.Rdata")
      mods$perm %>% table
      
      mods %>%
        filter(perm == 0)
      
      ggplot() +
        geom_boxplot(data = filter(mods,perm !=0 ),
                     aes(
                       x=mod,
                       y=-log10(p_val)
                     ))+
        geom_point(data = filter(mods,perm ==0 ),
                   aes(
                     x=mod,
                     y=-log10(p_val)
                   ), size = 3, color = "red") +
        coord_flip()->
        fig_ests
      
      ggsave(fig_ests, file = "fig_ests.pdf",
             w= 4, h=3.4)
      
      mods %>%
        filter(mod == "v_PRECTOTCORR") %>%
        filter(perm == 0) %>% .$p_val -> realda
      
      mods %>%
        filter(mod == "v_PRECTOTCORR") %>%
        filter(perm != 0) %>% .$p_val -> permsda
      ### Perms - perfirmance
      sum(permsda < realda)

      #####  
      ### Calculate morans' I
      AF_res.melt %>%
        filter(variable == "v_PRECTOTCORR") ->
        eco_var
      
      left_join(data_merged, eco_var) %>%
        filter(!is.na(value)) -> datMorI
      
      snp.dists <- as.matrix(dist(cbind(datMorI$m.Long,
                                        datMorI$m.Lat)))
      
      snp.dists.inv <- 1/snp.dists
      diag(snp.dists.inv) <- 0
      snp.dists.inv[is.infinite(snp.dists.inv)] =0
      
      Moran.I(datMorI$value, snp.dists.inv)
      
      #####      #####      #####
      #####
      #####
      #####
      #####      #####      ##### MAP
      #####
      world <- ne_countries(scale = 110, returnclass = "sf") 
      
      all.dataX <- getData(chr="X",
                           start=snppp,
                           end=snppp,
                           samplesUse=all.samps$sampleId) %>%
        left_join(all.samps) %>%
        mutate(year_pop = paste(year, city, sep = ":"))     
      
      all.dataX %>%
        group_by(city, cluster2.0_k4, continent) %>%
        summarise(
          m.AF = mean(af, na.rm = T),
          m.Lat = mean(lat),
          m.Long = mean(long)
        )-> cline_summaries
      
      ggplot() +
        geom_sf(data = world, 
                color = "black", 
                fill = "lightgray", linewidth = 0.1) +
        xlim(-128,160) + ylim(-40, 68) +
        geom_point(
          data = cline_summaries,
          size = 1.5,
          aes(
            y=m.Lat,
            x=m.Long,
            fill=m.AF,
            shape=as.factor(cluster2.0_k4)
          )) +
        #geom_vline(xintercept = -100) +
        scale_shape_manual(values = 21:24) +  
        scale_fill_gradient2(midpoint = 0.5)->
        lat_clineX
      
      ggsave(lat_clineX, file = "lat_clineX.pdf")
      
####
      ################
      ################
      ################      ################
      ################
      ################


left_join(F16s, VTbgd) %>% 
  filter(pos %in% top_muts_X$POS) %>%
  group_by(pos) %>%
  mutate(delta_bgd = VTbgd-value) %>%
  summarise(m.bgd = mean(delta_bgd)) %>%
  ggplot(aes(
    x=pos,
    y=m.bgd,
  )) + 
  geom_point() +
  geom_vline(xintercept = 15743371)->
  haploX
ggsave(haploX, file = "haploX.pdf")



### Individual SNP

snps.CROSS %>%
  filter(chr == "X") %>%
  filter(pos == 15743371) ->
  focal.CROSS

seqSetFilter(genofile_cross, variant.id=focal.CROSS$variant.id, 
             sample.id = c("SKF2_2_1",
                           "SKF3_3_1",
                           "SKF_1_1",
                           "SK_0_1",
                           "VT8F2_2_1",
                           "VT8F3_3_1",
                           "VT8F_1_1",
                           "VT8_0_1"))

data.table(ann=seqGetData(genofile_cross, "allele")) ->
  alleles.cross

message("Allele Freqs")
ad.cross <- seqGetData(genofile_cross, "annotation/format/AD")$data
dp.cross <- seqGetData(genofile_cross, "annotation/format/DP")

dat.cross <- ad.cross/dp.cross
colnames(dat.cross) <- paste(seqGetData(genofile_cross, "chromosome"), seqGetData(genofile_cross, "position"), sep="_")
rownames(dat.cross) <- seqGetData(genofile_cross, "sample.id")

dat.cross %>%
  as.data.frame() %>%
  mutate(sampleId = row.names(.)) %>%
  mutate(x_g = case_when(
    sampleId %in% c("SK_0_1") ~ "SK",
    sampleId %in% c("VT8_0_1") ~ "VT",
    TRUE ~ "cross"
  )) %>%
  group_by(x_g) %>%
  summarise(me.AF.c = mean(`X_15743371`),
            #med.AF.c = median(`2R_20551633`),
  )-> dat.cross.ag



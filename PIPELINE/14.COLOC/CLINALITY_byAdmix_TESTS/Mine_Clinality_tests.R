#### Latitudinal correlation with R in top SNP

library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)
library(forcats)
library(data.table)
library(gmodels)
library(ggrepel)
library(ape)

library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(rnaturalearth)
library(lmtest)

#library(sf)
#library(tidycensus)
#library(corrr)
#library(tmap)
#library(spdep)
#library(tigris)
##library(rmapshaper)
#library(flextable)
#library(car)
#library(spatialreg)
#library(stargazer)

####
#Get_DEST DATA
dest_loc <-"/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/dest.all.PoolSNP.001.50.24Aug2024.ann.gds"
genofile <- seqOpen(dest_loc)
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.csv")

snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))

snps.dt %>%
  filter(chr == "2R") %>%
  filter(pos == 20551633) ->
  focal

seqSetFilter(genofile,
             focal$variant.id)

data.table(ann=seqGetData(genofile, "annotation/info/ANN")) ->
  focal_ann

data.table(ann=seqGetData(genofile, "allele")) ->
  alleles


ad <- seqGetData(genofile, "annotation/format/AD")
ad <- ad$data
dp <- seqGetData(genofile, "annotation/format/DP")
print("Create dat object")

dat = ad/dp
####
colnames(dat) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , 
                       sep="_")
rownames(dat) <- seqGetData(genofile, "sample.id")


dat %>%
  as.data.frame() %>%
  mutate(sampleId = rownames(.)) %>%
  left_join(samps) ->
  focal_snp_data

### Bring in linear admix

lnadmix <- get(load("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/linear_admix/linear.Admix.Data.S2.Rdata"))
lnadmix %>%
  group_by(sampleId, admix.set, source_pop, lat, long, filter) %>%
  summarise(mean.est = mean(Estimate),
            sd.est = sd(Estimate)) %>%
  filter(filter == "silent") %>%
  filter(admix.set == "N.America" & source_pop == "AFRICA")  ->
  o_dat

o_dat %>%
  left_join(focal_snp_data, 
            by = c("sampleId","lat", "long")) %>% 
  filter(!is.na(`2R_20551633`))  -> data_merged
 
#### Bringin in weather data
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

#####
data_merged %>%
  group_by(city) %>%
  summarise(mean.AF = mean(`2R_20551633`),
            mean.mean.Est=mean(mean.est),
            ) %>%
  left_join(weather.ag, by = "city") ->
  data_merged.wag
 
lm(mean.AF ~ mean.mean.Est, data=data_merged.wag) ->
  ancestry_model
anova(ancestry_model)
summary(ancestry_model)
#ancestry_model =    -1.03824    0.21000  -4.944  2.2e-06 ***

data_merged.wag %>%
  as.data.frame() %>%
  mutate(AF_res = ancestry_model$residuals) ->
  data_merged.wag

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

save(mods, file = "clinal.mods.Rdata")
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

### nominal P value ANOVa
filter(AF_res.melt,
       variable == "v_PRECTOTCORR") ->
  dat.v_PRECTOTCORR

lm(AF_res ~ value, 
   data = dat.v_PRECTOTCORR) -> mod.precip
anova(mod.precip)

### Calculate morans' I
data_merged %>%
  group_by(city) %>%
  summarise(Lat = mean(lat),
            Long=mean(long),
  ) -> city.lat.long

left_join(dat.v_PRECTOTCORR, select(city.lat.long, city, Lat, Long) ) ->
  dat.v_PRECTOTCORR.lat.long

snp.dists <- as.matrix(dist(cbind(dat.v_PRECTOTCORR.lat.long$Long,
                                  dat.v_PRECTOTCORR.lat.long$Lat)))

snp.dists.inv <- 1/snp.dists
diag(snp.dists.inv) <- 0

snp.dists.inv[1:5, 1:5]

Moran.I(dat.v_PRECTOTCORR.lat.long$AF_res, snp.dists.inv)



####
AF_res.melt %>%
  filter(variable == "v_PRECTOTCORR") %>%
  ggplot(aes(
    x=value,
    y=AF_res,
    label=city
  )) +
  geom_smooth(method = "lm", se = F, color = "black") +
    geom_point( shape = 21, size = 4) + 
  geom_text_repel(size = 2.1) +
   theme_bw() ->
    fig_points
  
  ggsave(fig_points, file = "fig_points.pdf",
         w= 4, h=3.4)

  ##
  weather %>%
    filter(city %in% AF_res.melt$city) %>%
    group_by(city, LAT) %>%
    summarise(
      glob_var_PRECTOTCORR = sd(PRECTOTCORR, na.rm = T))->
    weather_glob
      
      
  weather %>%
    filter(city %in% AF_res.melt$city) %>%
    group_by(city, LAT, MO) %>%
    summarise(
      var_PRECTOTCORR = sd(PRECTOTCORR, na.rm = T)
    ) %>%
    left_join(weather_glob) %>%
    filter(MO %in% 5:11) %>%
    ggplot(
      aes(
        x=fct_reorder(city, glob_var_PRECTOTCORR),
        y=as.factor(MO),
        fill=var_PRECTOTCORR
      ) 
    ) +
    geom_point(size = 4, shape = 22, color = "white") +
    coord_flip() + theme_classic() +
    scale_fill_gradient2(
      low = "salmon", 
      mid = "grey", 
      high = "purple4", 
      midpoint = .6
    )->
    cline_weather
  
  ggsave(cline_weather, file = "cline_weather.pdf",
         h=4, w= 4.0)
  
#### CLINALITY ANALYSIS IN EUROPE AND AFRICA
#### CLINALITY ANALYSIS IN EUROPE AND AFRICA
#### CLINALITY ANALYSIS IN EUROPE AND AFRICA
#### CLINALITY ANALYSIS IN EUROPE AND AFRICA
#### CLINALITY ANALYSIS IN EUROPE AND AFRICA
  focal_snp_data %>%
    filter(continent == "Europe") %>%
    group_by(city,cluster1.0) %>%
    summarise(
      m.AF = mean(`2R_20551633`, na.rm = T),
      m.Lat = mean(lat),
      m.Long = mean(long)
    ) %>%
    filter(!is.na(cluster1.0))-> Eu_summaries

cor.test(~m.Lat+m.AF, data = Eu_summaries)
cor.test(~m.Long+m.AF, data = Eu_summaries)

  
lm(m.AF ~ cluster1.0, data=Eu_summaries) ->
  EU_clust_model
anova(EU_clust_model)
summary(EU_clust_model)

Eu_summaries %>%
  as.data.frame() %>%
  mutate(AF_res = EU_clust_model$residuals) ->
  Eu_summaries.wag

select_cols = c("m_RH2M","m_T2M","m_PRECTOTCORR","v_RH2M",
                "v_T2M", "v_PRECTOTCORR", "max_RH2M",
                "max_T2M", "max_PRECTOTCORR", "min_RH2M",
                "min_T2M", "min_PRECTOTCORR",
                "AF_res","city")

Eu_summaries.wag%>%
  left_join(weather.ag, by = "city") %>%
  .[,select_cols] %>%
  .[complete.cases(.),] %>% 
  melt(id = c("AF_res","city")) ->
  Eu_summaries.melt

mods.EU =
  foreach(k=0:100,
          .combine = "rbind",
          .errorhandling = "remove")%do%{
            foreach(i=select_cols[1:12],
                    .combine = "rbind",
                    .errorhandling = "remove")%do%{
                      message(paste(k,i, sep = "|"))
                      
                      data=filter(Eu_summaries.melt,
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

save(mods.EU, file = "clinal.modsEU.Rdata")
mods.EU$perm %>% table

mods.EU %>%
  filter(perm == 0)

ggplot() +
  geom_boxplot(data = filter(mods.EU,perm !=0 ),
               aes(
                 x=mod,
                 y=-log10(p_val)
               ))+
  geom_point(data = filter(mods.EU,perm ==0 ),
             aes(
               x=mod,
               y=-log10(p_val)
             ), size = 3, color = "red") +
  coord_flip()->
  fig_ests.mods.EU

ggsave(fig_ests.mods.EU, file = "fig_ests.mods.EU.pdf",
       w= 4, h=3.4)


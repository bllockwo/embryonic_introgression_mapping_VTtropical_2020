### libraries
  library(data.table)
  library(tidyverse)
  library(magrittr)
  library(lubridate)
  library(foreach)
  library(SeqArray)
  library(doMC)
  registerDoMC(8)
  library(fastglm)

###
weather <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/2023.DEST.2.0._release/Weather_NASA_power/WeatherData.mgarvin.csv")
weather %>%
mutate(Date = as.Date(paste(YEAR,MO,DY, sep = "-"))) -> weather.dat

###
samps <- fread("https://raw.githubusercontent.com/DEST-bio/DESTv2/refs/heads/main/populationInfo/dest_v2.samps_24Aug2024.csv")

samps %>% 
#filter(city %in% 
#c("Charlottesville",
#"Odesa",
#"Yesiloz"
#)) %>%
filter(set != "dgn") %>%
filter(Recommendation == "Pass") %>%
filter(!is.na(year)) %>%
filter(!is.na(min_day)) %>%
group_by(sampleId) %>%
mutate(DATE_spec = as.Date(str_split(sampleId, "_")[[1]][5],
format = "%Y-%m-%d")) -> seasonal.samps

###
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
#2R: 20,551,633 

  #CITY="Charlottesville"
  #CITY = "Odesa"
  #CITY = "Yesiloz"
  
#GLMs_from_mutation = 
#foreach(CITY = c("Charlottesville" , "Odesa"#, #"Yesiloz"
#                 ),
#.combine = "rbind",.errorhandling = "remove")%do%{
#%>% filter(year %in% 2016:2018) 
  
  SAMPStoUSE = seasonal.samps$sampleId
    
        data <- getData(chr="2R",
                        start=20551633,
                        end=20551633,
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
  							 data.table(mod=vars,
                             AIC=c(AIC(t1.year.var)),
                             b_temp=last(t1.year.var$coef), 
                             se_temp=last(t1.year.var$se),
                             p_lrt=anovaFun(t1.year.noperm, t1.year.var),
                             Timewindow = Slide_K,
                             perm=perm)
							
							
}
  print(paste(perm, Slide_K, vars, o$p_lrt, sep = "|"))
  return(o)  
}
}


#######
mod_o %>%
  slice_min(p_lrt) ->
  target_model

#######
mod_o_perms =
  foreach(Slide_K = 45, #seq(from=15, to = 90, by = 15),
          .combine = "rbind", .errorhandling = "remove")%do%{
            foreach(vars = c("T2M_sd"
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
              
              
              foreach(perm = 1:100, .combine = "rbind")%do%{
                
                value_perm = tmp[sample(1:dim(tmp)[1]), value]
                tmp_p = cbind(tmp, value_perm = value_perm)
                
                X.year.noperm<- model.matrix(~as.factor(year_pop), tmp_p)
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
####

sum(mod_o_perms$p_lrt < target_model$p_lrt)

rbind(target_model, mod_o_perms) -> mods_ploting
  ggplot() +
    geom_boxplot(data = filter(mods_ploting, perm != 0),
                aes(x=paste(mod,Timewindow, sep = "|"), y = -log10(p_lrt))) +
    geom_point(data = filter(mods_ploting, perm == 0),
               aes(x=paste(mod,Timewindow, sep = "|"), y = -log10(p_lrt),
               size = 3, color = "red")) -> model_plot
  ggsave(model_plot, file = "model_plot.pdf",
         w= 3, h = 3)

####
  weather_Select = weather_samps %>%
    filter(WinSlide == 45) %>%
    mutate(year_pop = paste(city,year, sep = ":"))
  
  data_downstream = left_join(data, weather_Select, by="sampleId") 
  
  data_downstream %>%
    filter(city %in% c("Charlottesville","Odesa","Yesiloz")) %>%
    ggplot(aes(
      x=sqrt(T2M_sd),
      y=af,
    )) +
    geom_point( shape = 21, size = 4) +
    geom_smooth(method = "lm", se = F, color = "black") +
    #scale_fill_gradient2(midpoint = 0, high = "firebrick") +
    facet_grid(.~city) + theme_bw() ->
    grand_plot
  ggsave(grand_plot, file = "grand_plot.pdf", w= 6, h = 2.5)



####
samps %>%
filter(city %in% "Charlottesville") %>% 
filter(year %in% 2016:2018) %>% .$sampleId -> cville.SAMPStoUSE
samps %>%
  filter(city %in% "Odesa") %>% .$sampleId -> Odesa.SAMPStoUSE

Cville.data <- getData(chr="2R",
                        start=20551633,
                        end=20551633,
                        samplesUse=cville.SAMPStoUSE) %>%
                        left_join(seasonal.samps) %>%
                        left_join(weather_Select) %>%
  filter(year %in% 2016:2018) %>%
  mutate(year_pop = paste(year, city, sep = ":"))

weather.dat %>%
  filter(city == "Charlottesville") %>%
  filter(YEAR %in% 2016:2018) ->
  cville_weather
names(cville_weather)[3] = "year"

cville_weather %>%
  group_by(year, MO) %>%
  summarise(T32c = sum(T2M > 32)) %>%
  mutate(Date = as.Date(paste(year,MO, 1, sep = "-")))->
  heatwave_Counter

cville_weather %>%
  group_by(year, MO) %>%
  summarise(T5c = sum(T2M < 4)) %>%
  mutate(Date = as.Date(paste(year,MO, 1, sep = "-")))->
  frost_Counter


ggplot() +
geom_line(
  data=mutate(cville_weather, type = "weather"),
aes(
  x=Date, y = T2M, color = T2M)) + 
  geom_line(
    data=mutate(heatwave_Counter, type = "weather2"),
    aes(
      x=Date, y = T32c), color = "firebrick",) +
  geom_line(
    data=mutate(frost_Counter, type = "weather3"),
    aes(
      x=Date, y = T5c), color = "steelblue",) +
  geom_smooth(
    data=mutate(Cville.data, type = "z.AFS"),
  aes(x=DATE_spec, y = af), color = "black", se = F,   linewidth = 0.5, linetype = "dashed") + 
  scale_color_gradient2(high = "red", low = "blue", mid = "grey") +
facet_grid(type~year, scales = "free") + theme_bw()->
perm_real_plot
ggsave(perm_real_plot, 
w = 9, h = 4,
file = "weather_Scaled_plot.pdf")


Odesa.data <- getData(chr="2R",
                       start=20551633,
                       end=20551633,
                       samplesUse=Odesa.SAMPStoUSE) %>%
  left_join(seasonal.samps) %>%
  left_join(weather_Select) %>%
  mutate(year_pop = paste(year, city, sep = ":")) %>%
  filter(!is.na(locality))


### #### All samples
samps %>% 
  filter(Recommendation == "Pass") %>%
  filter(!is.na(year)) %>%
  filter(!is.na(min_day)) %>%
  group_by(sampleId) %>%
  mutate(DATE_spec = as.Date(str_split(sampleId, "_")[[1]][5],
                             format = "%Y-%m-%d")) -> all.pass.samps

WinSlide=15

all.precip =    
foreach(i=1:dim(all.pass.samps)[1], 
                    .combine = "rbind",
        .errorhandling = "remove")%do%{
                      print(paste(WinSlide, i, sep = "|"))
                      
                      samp.i <- all.pass.samps[i,]
                      
                      bdate=samp.i$DATE_spec-WinSlide
                      fdate=samp.i$DATE_spec
                      weather.dat %>%
                        filter(city == samp.i$city) %>%
                        filter(Date >= bdate & Date <= fdate) %>%
                        summarize(
                          PRECTOTCORR_max = max(PRECTOTCORR),
                        ) -> weather.slice
                      data.frame(
                        samp.i,
                        WinSlide = WinSlide,
                        weather.slice)
                      
                    } ### inner

all.data <- getData(chr="2R",
                       start=20551633,
                       end=20551633,
                       samplesUse=all.precip$sampleId) %>%
  left_join(all.precip) %>%
  filter(!is.na(PRECTOTCORR_max)) %>%
  filter(!is.infinite(PRECTOTCORR_max)) %>%
  filter(!is.na(af))

ggplot(data = filter(all.data, sr_season == "fall"),
       aes(x=PRECTOTCORR_max, y= af)) +
  geom_point(aes(color = continent)) + 
  geom_smooth(method = "lm", color = "black")  ->
  Seasonal_plot_mut
ggsave(Seasonal_plot_mut, file = "PRECTOTCORR_max_world.pdf", w= 3, h = 3)


cor.test(~m.AF+m.Precip, data = all.data)

data = filter(all.data, continent == "Europe"))
cor.test(~m.AF+m.Precip, data = filter(all.data, continent == "North_America"))
cor.test(~af+PRECTOTCORR_max, data = all.data)

ggplot(data = all.data,
       aes(x=PRECTOTCORR_max, y= af)) +
  geom_point(aes(color = continent)) + 
  geom_smooth(method = "lm", color = "black")  ->
  Seasonal_plot_mut
ggsave(Seasonal_plot_mut, file = "PRECTOTCORR_max_world.pdf", w= 3, h = 3)
